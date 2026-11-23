#include <dynsyslab/math/ReebGraph.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
// Lower-star persistence: Edelsbrunner, Letscher & Zomorodian, Discrete Comput.
// Geom. 28 (2002) 511-533. H0 uses union-find with the elder rule; (H1,H2)
// pairs use the same rule on the dual graph in reverse filtration order, with
// one outside node for a boundary (Lefschetz duality).
// Join-tree sweep: Carr, Snoeyink & Axen, Comput. Geom. 24 (2003) 75-94.
// Link-based critical points: Banchoff, Amer. Math. Monthly 77 (1970) 475-485.
namespace dynsyslab::math {
namespace {
using Id=std::uint32_t;
constexpr Id none=std::numeric_limits<Id>::max();
constexpr std::size_t pieceLimit=20000000;
struct UnionFind {
    std::vector<Id> parent;
    explicit UnionFind(std::size_t n):parent(n){std::iota(parent.begin(),parent.end(),Id{0});}
    Id root(Id a){while(parent[a]!=a){parent[a]=parent[parent[a]];a=parent[a];}return a;}
};
struct Edge { Id a{},b{},t0{none},t1{none},hi{},lo{}; };
struct Triangle { std::array<Id,3> v{},e{}; Id hi{},mid{},lo{}; };
struct Complex {
    std::vector<Edge> edges;
    std::vector<Triangle> triangles;
    std::vector<Id> rank,byRank;
};
Id wrap(std::size_t x,std::size_t n){return Id(x==n?0:x);}
Complex buildComplex(const ScalarGrid2D& g){
    const std::size_t count=g.nx*g.ny;
    Complex c;c.rank.resize(count);c.byRank.resize(count);
    std::iota(c.byRank.begin(),c.byRank.end(),Id{0});
    std::sort(c.byRank.begin(),c.byRank.end(),[&](Id a,Id b){
        return g.values[a]!=g.values[b]?g.values[a]<g.values[b]:a<b;});
    for(std::size_t r=0;r<count;++r)c.rank[c.byRank[r]]=Id(r);
    std::vector<Id> slot(3*count,none);
    auto edge=[&](Id from,int direction,Id to,Id triangle){
        Id& id=slot[3*std::size_t(from)+direction];
        if(id==none){
            id=Id(c.edges.size());Edge e;e.a=from;e.b=to;
            e.hi=std::max(c.rank[from],c.rank[to]);e.lo=std::min(c.rank[from],c.rank[to]);
            c.edges.push_back(e);
        }
        auto& e=c.edges[id];(e.t0==none?e.t0:e.t1)=triangle;
        return id;
    };
    const std::size_t cx=g.periodicX?g.nx:g.nx-1,cy=g.periodicY?g.ny:g.ny-1;
    c.triangles.reserve(2*cx*cy);c.edges.reserve(3*count);
    auto add=[&](Id x,Id y,Id z,std::array<std::pair<Id,int>,3> from,std::array<Id,3> to){
        Triangle t;t.v={x,y,z};const Id id=Id(c.triangles.size());
        for(int k=0;k<3;++k)t.e[k]=edge(from[k].first,from[k].second,to[k],id);
        std::array<Id,3> r{c.rank[x],c.rank[y],c.rank[z]};std::sort(r.begin(),r.end());
        t.lo=r[0];t.mid=r[1];t.hi=r[2];c.triangles.push_back(t);
    };
    for(std::size_t j=0;j<cy;++j)for(std::size_t i=0;i<cx;++i){
        const Id i1=wrap(i+1,g.nx),j1=wrap(j+1,g.ny);
        const Id a=Id(i+g.nx*j),b=Id(i1+g.nx*j),d=Id(i+g.nx*j1),e=Id(i1+g.nx*j1);
        add(a,b,e,{{{a,0},{a,2},{b,1}}},{b,e,e}); // Directions: 0 x, 1 y, 2 diagonal.
        add(a,d,e,{{{a,1},{a,2},{d,0}}},{d,e,e});
    }
    return c;
}
std::vector<GridCriticalPoint> classify(const ScalarGrid2D& g,const std::vector<Id>& rank){
    // Cyclic link order of the Freudenthal star.
    static constexpr std::array<std::array<int,2>,6> ring{{{1,0},{1,1},{0,1},{-1,0},{-1,-1},{0,-1}}};
    std::vector<GridCriticalPoint> out(g.nx*g.ny);
    auto neighbour=[&](std::size_t i,std::size_t j,int dx,int dy)->std::optional<std::size_t>{
        long long x=(long long)i+dx,y=(long long)j+dy;const auto nx=(long long)g.nx,ny=(long long)g.ny;
        if(g.periodicX)x=(x+nx)%nx;else if(x<0||x>=nx)return std::nullopt;
        if(g.periodicY)y=(y+ny)%ny;else if(y<0||y>=ny)return std::nullopt;
        return std::size_t(x)+g.nx*std::size_t(y);
    };
    for(std::size_t j=0;j<g.ny;++j)for(std::size_t i=0;i<g.nx;++i){
        const auto v=i+g.nx*j;std::array<int,6> state{}; // 0 missing, 1 lower, 2 upper.
        int present=0,lower=0;
        for(int k=0;k<6;++k)if(auto w=neighbour(i,j,ring[k][0],ring[k][1])){
            ++present;state[k]=rank[*w]<rank[v]?1:2;lower+=state[k]==1;
        }
        unsigned lowerRuns=0,upperRuns=0;
        for(int k=0;k<6;++k){
            const int here=state[k],next=state[(k+1)%6];
            if(here==1&&next!=1)++lowerRuns;
            if(here==2&&next!=2)++upperRuns;
        }
        const bool boundary=present<6;
        if(!boundary&&lower==6)lowerRuns=1; // Whole circle has no run transition.
        if(!boundary&&lower==0)upperRuns=1;
        auto& p=out[v];p.vertex=v;p.i=i;p.j=j;p.value=g.values[v];p.onBoundary=boundary;
        p.lowerLinkComponents=lowerRuns;p.upperLinkComponents=upperRuns;
        p.eulerContribution=!boundary&&lower==6?1:1-int(lowerRuns);
        if(lowerRuns==0)p.kind=GridCriticalKind::Minimum;
        else if(upperRuns==0)p.kind=GridCriticalKind::Maximum;
        else if(lowerRuns>1||upperRuns>1){p.kind=GridCriticalKind::Saddle;p.multiplicity=std::max(lowerRuns,upperRuns)-1;}
    }
    return out;
}
std::vector<GridPersistencePair> persistence(const ScalarGrid2D& g,const Complex& c,double minimum){
    std::vector<GridPersistencePair> out;
    auto finitePair=[&](unsigned dimension,Id birth,Id death){
        if(birth==death)return; // Paired inside one lower star.
        if(g.values[death]-g.values[birth]<minimum)return;
        out.push_back({dimension,birth,g.values[birth],std::size_t(death),g.values[death]});
    };
    std::vector<Id> edgeOrder(c.edges.size());std::iota(edgeOrder.begin(),edgeOrder.end(),Id{0});
    std::sort(edgeOrder.begin(),edgeOrder.end(),[&](Id x,Id y){
        const auto& a=c.edges[x];const auto& b=c.edges[y];
        return a.hi!=b.hi?a.hi<b.hi:a.lo!=b.lo?a.lo<b.lo:x<y;});
    UnionFind vertices(g.nx*g.ny);std::vector<Id> oldest(g.nx*g.ny);std::iota(oldest.begin(),oldest.end(),Id{0});
    std::vector<char> negative(c.edges.size());
    for(Id id:edgeOrder){
        const auto& e=c.edges[id];Id a=vertices.root(e.a),b=vertices.root(e.b);
        if(a==b)continue;
        negative[id]=1;
        if(c.rank[oldest[a]]>c.rank[oldest[b]])std::swap(a,b); // a holds the elder minimum.
        finitePair(0,oldest[b],c.byRank[e.hi]);
        vertices.parent[b]=a;
    }
    const Id globalMinimum=c.byRank.front();
    out.push_back({0,globalMinimum,g.values[globalMinimum],std::nullopt,std::nullopt});
    // Dual sweep: triangles are dual vertices, born in reverse filtration order.
    const Id triangleCount=Id(c.triangles.size());
    std::vector<Id> triangleOrder(triangleCount);std::iota(triangleOrder.begin(),triangleOrder.end(),Id{0});
    std::sort(triangleOrder.begin(),triangleOrder.end(),[&](Id x,Id y){
        const auto& a=c.triangles[x];const auto& b=c.triangles[y];
        if(a.hi!=b.hi)return a.hi<b.hi;
        if(a.mid!=b.mid)return a.mid<b.mid;
        return a.lo!=b.lo?a.lo<b.lo:x<y;});
    std::vector<Id> position(triangleCount+1);
    for(Id k=0;k<triangleCount;++k)position[triangleOrder[k]]=k;
    position[triangleCount]=triangleCount; // The outside node is elder to every triangle.
    UnionFind dual(triangleCount+1);std::vector<Id> elder(triangleCount+1);std::iota(elder.begin(),elder.end(),Id{0});
    for(auto it=edgeOrder.rbegin();it!=edgeOrder.rend();++it){
        if(negative[*it])continue;
        const auto& e=c.edges[*it];
        Id a=dual.root(e.t0),b=dual.root(e.t1==none?triangleCount:e.t1);
        const Id birth=c.byRank[e.hi];
        if(a==b){out.push_back({1,birth,g.values[birth],std::nullopt,std::nullopt});continue;}
        if(position[elder[a]]<position[elder[b]])std::swap(a,b); // a holds the elder dual class.
        finitePair(1,birth,c.byRank[c.triangles[elder[b]].hi]);
        dual.parent[b]=a;
    }
    const Id top=elder[dual.root(triangleOrder.back())];
    if(top!=triangleCount){ // Closed surface: the fundamental class survives.
        const Id vertex=c.byRank[c.triangles[top].hi];
        out.push_back({2,vertex,g.values[vertex],std::nullopt,std::nullopt});
    }
    return out;
}
LevelGraphNode nodeFor(const GridCriticalPoint& p){return {p.vertex,p.value,p.kind,p.onBoundary,0,0};}
void countDegrees(ScalarLevelGraph& g){
    for(auto& n:g.nodes)n.lowerDegree=n.upperDegree=0;
    for(const auto& a:g.arcs){++g.nodes[a.lowerNode].upperDegree;++g.nodes[a.upperNode].lowerDegree;}
}
// Join tree of sublevel components (Carr, Snoeyink & Axen 2003).
ScalarLevelGraph mergeTree(const Complex& c,const std::vector<GridCriticalPoint>& critical){
    const std::size_t count=c.rank.size();
    std::vector<std::vector<Id>> lower(count);
    for(const auto& e:c.edges)lower[e.hi].push_back(c.byRank[e.lo]);
    UnionFind sets(count);std::vector<Id> head(count,none);
    ScalarLevelGraph out;
    auto newNode=[&](Id vertex){out.nodes.push_back(nodeFor(critical[vertex]));return Id(out.nodes.size()-1);};
    for(std::size_t r=0;r<count;++r){
        const Id v=c.byRank[r];std::vector<Id> roots;
        for(Id w:lower[r]){const Id root=sets.root(w);if(std::find(roots.begin(),roots.end(),root)==roots.end())roots.push_back(root);}
        Id current=none;
        if(roots.empty())current=newNode(v);
        else if(roots.size()==1)current=head[roots[0]];
        else{current=newNode(v);for(Id root:roots)out.arcs.push_back({head[root],current});}
        for(Id root:roots)sets.parent[root]=v;
        head[v]=current;
    }
    const Id top=c.byRank.back(),last=head[sets.root(top)];
    if(out.nodes[last].vertex!=top){const Id root=newNode(top);out.arcs.push_back({last,root});}
    countDegrees(out);
    return out;
}
// Reeb graph by slab components. Between consecutive critical ranks every
// contour is a regular cylinder, and each triangle meets a slab in one convex
// piece. Pieces join across shared edges; a component is attached to a critical
// node when it touches that vertex's star, otherwise it continues through the level.
ScalarLevelGraph reebGraph(const Complex& c,const std::vector<GridCriticalPoint>& critical){
    std::vector<Id> ranks;
    for(std::size_t r=0;r<c.byRank.size();++r)if(critical[c.byRank[r]].kind!=GridCriticalKind::Regular)ranks.push_back(Id(r));
    ScalarLevelGraph out;
    for(Id r:ranks)out.nodes.push_back(nodeFor(critical[c.byRank[r]]));
    auto slabRange=[&](Id low,Id high){
        const auto first=std::size_t(std::upper_bound(ranks.begin(),ranks.end(),low)-ranks.begin())-1;
        const auto last=std::size_t(std::lower_bound(ranks.begin(),ranks.end(),high)-ranks.begin())-1;
        return std::pair<std::size_t,std::size_t>{first,last};
    };
    const std::size_t triangleCount=c.triangles.size();
    std::vector<std::size_t> start(triangleCount+1);std::vector<Id> firstSlab(triangleCount);
    for(std::size_t t=0;t<triangleCount;++t){
        const auto [a,b]=slabRange(c.triangles[t].lo,c.triangles[t].hi);
        firstSlab[t]=Id(a);start[t+1]=start[t]+(b-a+1);
        if(start[t+1]>pieceLimit)throw std::invalid_argument("Reeb extraction exceeds work limit; coarsen or smooth the grid");
    }
    const std::size_t pieces=start[triangleCount];
    auto piece=[&](Id t,std::size_t slab){return Id(start[t]+slab-firstSlab[t]);};
    auto slabCount=[&](Id t){return start[t+1]-start[t];};
    UnionFind sets(pieces);
    for(const auto& e:c.edges){
        if(e.t1==none)continue;
        const auto [a,b]=slabRange(e.lo,e.hi);
        for(std::size_t s=a;s<=b;++s){
            const Id x=sets.root(piece(e.t0,s)),y=sets.root(piece(e.t1,s));
            if(x!=y)sets.parent[y]=x;
        }
    }
    std::vector<char> hasLower(pieces),hasUpper(pieces);
    std::vector<Id> criticalIndex(c.rank.size(),none);
    for(std::size_t k=0;k<ranks.size();++k)criticalIndex[c.byRank[ranks[k]]]=Id(k);
    for(Id t=0;t<triangleCount;++t){
        const auto& tri=c.triangles[t];
        for(Id v:tri.v){
            const Id k=criticalIndex[v];if(k==none)continue;
            if(c.rank[v]<tri.hi)hasLower[sets.root(piece(t,k))]=1;
            if(c.rank[v]>tri.lo)hasUpper[sets.root(piece(t,k-1))]=1;
        }
    }
    // Regular contours continue across a critical level as one arc.
    UnionFind chains(pieces);
    for(Id t=0;t<triangleCount;++t)for(std::size_t k=0;k+1<slabCount(t);++k){
        const std::size_t s=firstSlab[t]+k;
        const Id a=sets.root(piece(t,s)),b=sets.root(piece(t,s+1));
        if(hasUpper[a]||hasLower[b])continue;
        const Id x=chains.root(a),y=chains.root(b);
        if(x!=y)chains.parent[y]=x;
    }
    std::vector<Id> slabOf(pieces,none);
    for(Id t=0;t<triangleCount;++t)for(std::size_t k=0;k<slabCount(t);++k)
        slabOf[sets.root(Id(start[t]+k))]=Id(firstSlab[t]+k);
    std::vector<Id> from(pieces,none),to(pieces,none);
    for(Id p=0;p<pieces;++p){
        if(sets.root(p)!=p)continue;
        const Id chain=chains.root(p);
        if(hasLower[p]){
            if(from[chain]!=none)throw std::logic_error("Reeb arc has two lower nodes");
            from[chain]=slabOf[p];
        }
        if(hasUpper[p]){
            if(to[chain]!=none)throw std::logic_error("Reeb arc has two upper nodes");
            to[chain]=slabOf[p]+1;
        }
    }
    for(Id p=0;p<pieces;++p){
        if(sets.root(p)!=p||chains.root(p)!=p)continue;
        if(from[p]==none||to[p]==none)throw std::logic_error("Reeb arc is missing an endpoint");
        out.arcs.push_back({from[p],to[p]});
    }
    std::sort(out.arcs.begin(),out.arcs.end(),[](const auto& a,const auto& b){
        return a.lowerNode!=b.lowerNode?a.lowerNode<b.lowerNode:a.upperNode<b.upperNode;});
    countDegrees(out);
    return out;
}
void validate(const ScalarGrid2D& g,double minimum){
    if(g.nx<2||g.ny<2||g.nx>1000000||g.ny>1000000||g.nx*g.ny>1000000)
        throw std::invalid_argument("Scalar grid needs nx,ny>=2 and at most 1000000 vertices");
    if((g.periodicX&&g.nx<3)||(g.periodicY&&g.ny<3))throw std::invalid_argument("Periodic grid axes need at least 3 samples");
    if(g.values.size()!=g.nx*g.ny)throw std::invalid_argument("Scalar grid needs nx*ny values");
    for(double v:g.values)if(!std::isfinite(v))throw std::invalid_argument("Scalar grid values must be finite");
    if(!std::isfinite(minimum)||minimum<0)throw std::invalid_argument("Minimum persistence must be finite and nonnegative");
}
}
const char* gridCriticalKindName(GridCriticalKind kind){
    switch(kind){
        case GridCriticalKind::Minimum:return "minimum";
        case GridCriticalKind::Saddle:return "saddle";
        case GridCriticalKind::Maximum:return "maximum";
        default:return "regular";
    }
}
FiniteGraph toFiniteGraph(const ScalarLevelGraph& graph){
    FiniteGraph out;out.vertices=graph.nodes.size();out.edges.reserve(graph.arcs.size());
    for(const auto& a:graph.arcs)out.edges.push_back({a.lowerNode,a.upperNode,
        std::abs(graph.nodes[a.upperNode].value-graph.nodes[a.lowerNode].value)});
    return out;
}
ScalarTopology analyzeScalarGrid(const ScalarGrid2D& grid,double minimumPersistence){
    validate(grid,minimumPersistence);
    const auto complex=buildComplex(grid);
    const auto all=classify(grid,complex.rank);
    ScalarTopology out;
    for(const auto& p:all){
        out.eulerCharacteristic+=p.eulerContribution;
        if(p.kind!=GridCriticalKind::Regular)out.criticalPoints.push_back(p);
    }
    out.persistence=persistence(grid,complex,minimumPersistence);
    out.mergeTree=mergeTree(complex,all);
    out.reebGraph=reebGraph(complex,all);
    out.reebFirstBetti=graphTopology(toFiniteGraph(out.reebGraph)).firstBettiNumber;
    for(const auto& n:out.reebGraph.nodes){
        if(n.kind!=GridCriticalKind::Saddle)continue;
        if(n.onBoundary&&n.lowerDegree+n.upperDegree<3)continue; // Contour reshaped by the boundary only.
        out.separatrices.push_back({n.vertex,n.value,n.onBoundary,n.lowerDegree,n.upperDegree});
    }
    return out;
}
}
