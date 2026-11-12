#include "local_rhs.hpp"
#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>
struct Node {
  int point, axis, left, right;
};
struct Hit {
  double dist;
  int point;
  bool operator<(const Hit &h) const { return dist < h.dist; }
};
struct Field {
  int n, k;
  double damping;
  const double *c, *m, *p, *q, *coef, *scale;
  std::vector<Node> nodes;
  Field(int N, int rows, int K, double d, const double *C, const double *M,
        const double *P, const double *Q, const double *A, const double *S)
      : n(N), k(K), damping(d), c(C), m(M), p(P), q(Q), coef(A), scale(S) {
    std::vector<int> idx(rows);
    std::iota(idx.begin(), idx.end(), 0);
    nodes.reserve(rows);
    build(idx, 0, rows);
  }
  int build(std::vector<int> &ids, int lo, int hi) {
    if (lo == hi)
      return -1;
    int axis = 0;
    double widest = -1;
    for (int j = 0; j < n; j++) {
      double a = 1e300, b = -1e300;
      for (int i = lo; i < hi; i++) {
        double v = c[ids[i] * n + j];
        a = std::min(a, v);
        b = std::max(b, v);
      }
      if (b - a > widest) {
        axis = j;
        widest = b - a;
      }
    }
    int mid = (lo + hi) / 2;
    std::nth_element(
        ids.begin() + lo, ids.begin() + mid, ids.begin() + hi,
        [&](int a, int b) { return c[a * n + axis] < c[b * n + axis]; });
    int id = (int)nodes.size();
    nodes.push_back({ids[mid], axis, -1, -1});
    int left = build(ids, lo, mid), right = build(ids, mid + 1, hi);
    nodes[id].left = left;
    nodes[id].right = right;
    return id;
  }
  void query(int node, const double *x, Hit *h, int &count) const {
    if (node < 0)
      return;
    const auto &t = nodes[node];
    const double *C = c + t.point * n;
    double ds = 0;
    for (int j = 0; j < n; j++) {
      double v = x[j] - C[j];
      ds += v * v;
    }
    if (count < k + 1) {
      h[count++] = {ds, t.point};
      std::push_heap(h, h + count);
    } else if (ds < h[0].dist) {
      std::pop_heap(h, h + count);
      h[count - 1] = {ds, t.point};
      std::push_heap(h, h + count);
    }
    double split = x[t.axis] - C[t.axis];
    int near = split < 0 ? t.left : t.right, far = split < 0 ? t.right : t.left;
    query(near, x, h, count);
    if (count < k + 1 || split * split <= h[0].dist)
      query(far, x, h, count);
  }
  void rhs(const double *b, double *out) const {
    double x[64], dx[64];
    Hit hits[65];
    int count = 0;
    for (int j = 0; j < n; j++) {
      x[j] = b[j] / scale[j];
      out[j] = 0;
    }
    query(0, x, hits, count);
    double radius = std::max(1e-12, std::sqrt(hits[0].dist)), sum = 0;
    const bool tied = std::all_of(hits+1,hits+count,[&](const Hit& h) {
      return std::sqrt(h.dist) >= radius*(1.-1e-6);
    });
    for (int h = 1; h < count; h++) {
      int idx = hits[h].point;
      double ratio = std::min(1., std::sqrt(hits[h].dist) / radius),
             u = 1 - ratio;
      double weight = u * u * u * u * (1 + 4 * ratio), a = 0, v = 0, A = 0,
             B = 0;
      if (tied) weight = 1.;
      sum += weight;
      const double *P = p + idx * n * 2, *Q = q + idx * n * 2, *M = m + idx * n,
                   *C = coef + idx * 6 * n;
      for (int j = 0; j < n; j++) {
        dx[j] = x[j] - M[j];
        a += dx[j] * P[2 * j];
        v += dx[j] * P[2 * j + 1];
        A += dx[j] * Q[2 * j];
        B += dx[j] * Q[2 * j + 1];
      }
      double f[] = {1., a, v, a * a, a * v, v * v};
      for (int j = 0; j < n; j++) {
        double r = 0;
        for (int i = 0; i < 6; i++)
          r += f[i] * C[i * n + j];
        out[j] +=
            weight * (r - damping * (dx[j] - A * Q[2 * j] - B * Q[2 * j + 1]));
      }
    }
    for (int j = 0; j < n; j++)
      out[j] = out[j] * scale[j] / sum;
  }
};
extern "C" {
void *mc_local_field_v1(int n, int rows, int k, double d, const double *c,
                        const double *m, const double *p, const double *q,
                        const double *a, const double *s) {
  if (n < 2 || n > 64 || k < 1 || k > 64 || k >= rows || !std::isfinite(d) ||
      d < 0 || !c || !m || !p || !q || !a || !s)
    return nullptr;
  try { return new Field(n, rows, k, d, c, m, p, q, a, s); }
  catch (...) { return nullptr; }
}
void mc_local_batch_v1(void *handle, int count, const double *x, double *y) {
  auto *f = (Field *)handle;
  for (int i = 0; i < count; i++)
    f->rhs(x + i * f->n, y + i * f->n);
}
void mc_local_drop_v1(void *h) { delete (Field *)h; }
}
