#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>

namespace {
struct Trace {
    const double *x, *y;
    std::ptrdiff_t sx, sy;
    double tolX, tolY;
    std::int64_t* output;
    std::size_t count=0;

    // Every accepted segment has been checked against every intervening point.
    // A midpoint split keeps recursion shallow even across a discontinuity.
    void appendSegment(std::size_t a,std::size_t b) {
        if(b>a+1) {
            const double x0=x[sx*std::ptrdiff_t(a)], y0=y[sy*std::ptrdiff_t(a)];
            const double dx=x[sx*std::ptrdiff_t(b)]-x0, dy=y[sy*std::ptrdiff_t(b)]-y0;
            for(std::size_t i=a+1;i<b;++i) {
                const double f=double(i-a)/double(b-a);
                if(std::abs(x[sx*std::ptrdiff_t(i)]-(x0+f*dx))>tolX ||
                   std::abs(y[sy*std::ptrdiff_t(i)]-(y0+f*dy))>tolY) {
                    const auto middle=a+(b-a)/2;
                    appendSegment(a,middle);
                    appendSegment(middle,b);
                    return;
                }
            }
        }
        output[count++]=static_cast<std::int64_t>(b);
    }
};
}

// Strides are in doubles. Output holds at most n indices in chronological order.
// The caller supplies finite arrays. Returning zero indicates invalid arguments.
extern "C" std::size_t trajectory_vertices(const double* x,std::ptrdiff_t sx,
        const double* y,std::ptrdiff_t sy,std::size_t n,double tolerance,std::int64_t* indices) {
    if(!n || !x || !y || !indices || !(tolerance>=0 && tolerance<=1))return 0;
    std::array<std::size_t,6> anchors{0,0,0,0,0,n-1};
    for(std::ptrdiff_t i=0;i<static_cast<std::ptrdiff_t>(n);++i) {
        const auto ax=x[sx*i], ay=y[sy*i];
        if(!std::isfinite(ax) || !std::isfinite(ay))return 0;
        if(ax<x[sx*std::ptrdiff_t(anchors[1])])anchors[1]=i;
        if(ax>x[sx*std::ptrdiff_t(anchors[2])])anchors[2]=i;
        if(ay<y[sy*std::ptrdiff_t(anchors[3])])anchors[3]=i;
        if(ay>y[sy*std::ptrdiff_t(anchors[4])])anchors[4]=i;
    }
    Trace trace{x,y,sx,sy,(x[sx*std::ptrdiff_t(anchors[2])]-x[sx*std::ptrdiff_t(anchors[1])])*tolerance,
                          (y[sy*std::ptrdiff_t(anchors[4])]-y[sy*std::ptrdiff_t(anchors[3])])*tolerance,indices};
    std::sort(anchors.begin(),anchors.end());
    const auto end=std::unique(anchors.begin(),anchors.end());
    indices[trace.count++]=0;
    for(auto p=anchors.begin()+1;p<end;++p) {
        for(std::size_t start=p[-1];start<*p;) {
            const auto stop=std::min(*p,start+8192);
            trace.appendSegment(start,stop);
            start=stop;
        }
    }
    return trace.count;
}
