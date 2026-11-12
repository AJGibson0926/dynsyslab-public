#include "crossing_kernel.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
namespace {
struct Segment { std::size_t i; double xmin,xmax,ymin,ymax; };
double det(double x,double y,double u,double v){return x*v-y*u;}
}
extern "C" int mc_crossings_v1(const double* p,std::size_t n,mc_crossing_v1* out,
                              std::size_t cap,std::size_t* found){
    if(!found)return 1;
    *found=0;
    if(!p||n<2||n>std::numeric_limits<std::size_t>::max()/2)return 1;
    for(std::size_t j=0;j<2*n;++j)if(!std::isfinite(p[j]))return 1;
    try {
        std::vector<Segment> segments,active;
        for(std::size_t i=0;i+1<n;++i){
            auto a=p+2*i,b=a+2;
            if(a[0]!=b[0]||a[1]!=b[1])
                segments.push_back({i,std::min(a[0],b[0]),std::max(a[0],b[0]),std::min(a[1],b[1]),std::max(a[1],b[1])});
        }
        std::sort(segments.begin(),segments.end(),[](const auto& a,const auto& b){
            return a.xmin==b.xmin ? a.i<b.i : a.xmin<b.xmin;
        });
        for(const auto& s:segments){
            active.erase(std::remove_if(active.begin(),active.end(),[&](const auto& a){return a.xmax<s.xmin;}),active.end());
            for(const auto& q:active){
                auto i=std::min(s.i,q.i),j=std::max(s.i,q.i);
                if(j-i<=1||s.ymax<q.ymin||q.ymax<s.ymin)continue;
                const auto a=p+2*i,b=p+2*j;
                double x=a[2]-a[0],y=a[3]-a[1],u=b[2]-b[0],v=b[3]-b[1],d=det(x,y,u,v);
                if(std::abs(d)<=32*std::numeric_limits<double>::epsilon()*std::hypot(x,y)*std::hypot(u,v))continue;
                double dx=b[0]-a[0],dy=b[1]-a[1],f=det(dx,dy,u,v)/d,g=det(dx,dy,x,y)/d;
                // Assign shared vertices once using half-open segments; the final endpoint is included.
                if(f<0||g<0||f>=1||g>1||(g==1&&j+2<n))continue;
                if(out && *found>=cap)return 2;
                if(out)out[*found]={i,j,f,g};
                ++*found;
            }
            active.push_back(s);
        }
    }catch(...){return 3;}
    return 0;
}
