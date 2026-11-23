#include <dynsyslab/math/Topology.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
// Graph homology: Hatcher, Algebraic Topology, section 1.A.
// H0 persistence uses union-find and the elder rule; Edelsbrunner and Harer,
// Computational Topology (AMS, 2010), persistent homology chapters.
// Shortest paths: Dijkstra, Numerische Mathematik 1 (1959), 269-271,
// doi:10.1007/BF01386390; nonnegative undirected edge weights.
namespace dynsyslab::math {
namespace {
void validate(const FiniteGraph& g){
    if(g.vertices>1000000||g.edges.size()>10000000)throw std::invalid_argument("Graph exceeds analysis limits");
    for(const auto& e:g.edges)if(e.from>=g.vertices||e.to>=g.vertices||!std::isfinite(e.length)||e.length<0)
        throw std::invalid_argument("Graph edge endpoint or nonnegative length invalid");
}
struct UnionFind {
    std::vector<std::size_t> parent,size;
    explicit UnionFind(std::size_t n):parent(n),size(n,1){std::iota(parent.begin(),parent.end(),0);}
    std::size_t root(std::size_t a){while(parent[a]!=a){parent[a]=parent[parent[a]];a=parent[a];}return a;}
    bool join(std::size_t a,std::size_t b){a=root(a);b=root(b);if(a==b)return false;
        if(size[a]<size[b])std::swap(a,b);parent[b]=a;size[a]+=size[b];return true;}
};
}
GraphTopology graphTopology(const FiniteGraph& g){
    validate(g);UnionFind u(g.vertices);GraphTopology r;
    r.degree.resize(g.vertices);r.component.resize(g.vertices);r.connectedComponents=g.vertices;
    for(std::size_t i=0;i<g.edges.size();++i){auto e=g.edges[i];++r.degree[e.from];++r.degree[e.to];
        if(u.join(e.from,e.to)){--r.connectedComponents;r.spanningForestEdges.push_back(i);}
        else r.cycleClosingEdges.push_back(i);}
    r.firstBettiNumber=r.cycleClosingEdges.size();
    std::vector<std::size_t> labels(g.vertices,g.vertices);std::size_t next=0;
    for(std::size_t i=0;i<g.vertices;++i){auto root=u.root(i);if(labels[root]==g.vertices)labels[root]=next++;r.component[i]=labels[root];}
    return r;
}
std::vector<std::optional<double>> graphDistances(const FiniteGraph& g,std::size_t source){
    validate(g);if(source>=g.vertices)throw std::invalid_argument("Graph source out of range");
    using Item=std::pair<double,std::size_t>;
    std::vector<std::vector<Item>> adjacent(g.vertices);
    for(const auto& e:g.edges){adjacent[e.from].push_back({e.length,e.to});adjacent[e.to].push_back({e.length,e.from});}
    std::vector<std::optional<double>> distance(g.vertices);distance[source]=0;
    std::priority_queue<Item,std::vector<Item>,std::greater<Item>> queue;queue.push({0,source});
    while(!queue.empty()){auto [d,v]=queue.top();queue.pop();if(d>*distance[v])continue;
        for(auto [length,w]:adjacent[v]){const double candidate=d+length;
            if(!std::isfinite(candidate))throw std::overflow_error("Graph path length overflow");
            if(!distance[w]||candidate<*distance[w]){distance[w]=candidate;queue.push({candidate,w});}}}
    return distance;
}
std::vector<PersistenceInterval> lowerStarPersistence(const FiniteGraph& g,const std::vector<double>& values){
    validate(g);if(values.size()!=g.vertices)throw std::invalid_argument("One value required per graph vertex");
    for(auto v:values)if(!std::isfinite(v))throw std::invalid_argument("Filtration values must be finite");
    std::vector<std::size_t> edges(g.edges.size()),birth(g.vertices);std::iota(edges.begin(),edges.end(),0);std::iota(birth.begin(),birth.end(),0);
    auto level=[&](std::size_t i){return std::max(values[g.edges[i].from],values[g.edges[i].to]);};
    std::stable_sort(edges.begin(),edges.end(),[&](auto a,auto b){return level(a)<level(b);});
    UnionFind u(g.vertices);std::vector<PersistenceInterval> result;
    for(auto i:edges){auto a=u.root(g.edges[i].from),b=u.root(g.edges[i].to);if(a==b)continue;
        auto older=[&](std::size_t x,std::size_t y){return std::pair{values[birth[x]],birth[x]}<std::pair{values[birth[y]],birth[y]};};
        if(!older(a,b))std::swap(a,b);
        result.push_back({birth[b],values[birth[b]],level(i)});
        auto elder=birth[a];u.join(a,b);birth[u.root(a)]=elder;
    }
    for(std::size_t i=0;i<g.vertices;++i)if(u.root(i)==i)result.push_back({birth[i],values[birth[i]],{}});
    return result;
}
}
