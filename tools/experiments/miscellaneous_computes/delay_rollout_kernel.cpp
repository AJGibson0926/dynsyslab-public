#include "delay_rollout_kernel.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

extern "C" int mc_delay_rk4_v1(mc_delay_rhs_v1 rhs, void* context, std::size_t dim,
    const double* lags, std::size_t taps, const double* history, std::size_t history_rows,
    double h, std::size_t output_rows, const double* limits, double* output, std::size_t* written) {
    if(!written)return 3;
    *written=0;
    if(!rhs||!dim||!taps||!history||!lags||!output||history_rows<2||!output_rows||
       !std::isfinite(h)||h<=0.||lags[0]!=0.)return 3;
    const std::size_t max=std::numeric_limits<std::size_t>::max();
    if(std::max({taps,history_rows,output_rows})>max/dim || output_rows>max-history_rows)return 3;
    for(std::size_t j=1;j<taps;++j)
        if(!std::isfinite(lags[j])||lags[j]<h||lags[j]<=lags[j-1]||lags[j]/h>double(history_rows-1))return 3;
    for(std::size_t j=0;j<history_rows*dim;++j)if(!std::isfinite(history[j]))return 3;
    for(std::size_t j=0;j<dim;++j)if(limits&&(!(limits[j]>0.)||std::isnan(limits[j])))return 3;
    try {
        std::vector<double> window(taps*dim), stage(dim), k1(dim),k2(dim),k3(dim),k4(dim);
        struct Delay {std::size_t behind;double frac;};
        std::vector<Delay> delays[3];
        for(int phase=0;phase<3;++phase){
            for(std::size_t j=1;j<taps;++j){
                double distance=lags[j]/h-.5*phase;
                double whole=std::floor(distance);
                delays[phase].push_back({static_cast<std::size_t>(whole),distance-whole});
            }
        }
        auto row=[&](std::size_t index){return index<history_rows?history+index*dim:static_cast<const double*>(output)+(index-history_rows+1)*dim;};
        auto check=[&](const double* x){
            for(std::size_t j=0;j<dim;++j){if(!std::isfinite(x[j]))return 2;if(limits&&std::abs(x[j])>=limits[j])return 1;}
            return 0;
        };
        const double* first=history+(history_rows-1)*dim;
        if(int code=check(first))return code;
        std::copy_n(first,dim,output);*written=1;
        for(std::size_t i=0;i+1<output_rows;++i){
            const double* x=output+i*dim;
            const std::size_t last=history_rows-1+i;
            auto evaluate=[&](const double* current,int phase,double* velocity){
                if(int code=check(current))return code;
                std::copy_n(current,dim,window.data());
                for(std::size_t j=1;j<taps;++j){
                    auto d=delays[phase][j-1];std::size_t late=last-d.behind;
                    const double* b=row(late);const double* a=d.frac==0.?b:row(late-1);
                    for(std::size_t v=0;v<dim;++v)window[j*dim+v]=b[v]+d.frac*(a[v]-b[v]);
                }
                if(rhs(context,window.data(),velocity))return 2;
                for(std::size_t v=0;v<dim;++v)if(!std::isfinite(velocity[v]))return 2;
                return 0;
            };
            if(int c=evaluate(x,0,k1.data()))return c;
            for(std::size_t j=0;j<dim;++j)stage[j]=x[j]+.5*h*k1[j];
            if(int c=evaluate(stage.data(),1,k2.data()))return c;
            for(std::size_t j=0;j<dim;++j)stage[j]=x[j]+.5*h*k2[j];
            if(int c=evaluate(stage.data(),1,k3.data()))return c;
            for(std::size_t j=0;j<dim;++j)stage[j]=x[j]+h*k3[j];
            if(int c=evaluate(stage.data(),2,k4.data()))return c;
            for(std::size_t j=0;j<dim;++j)stage[j]=x[j]+(h/6.)*(k1[j]+2.*k2[j]+2.*k3[j]+k4[j]);
            if(int c=check(stage.data()))return c;
            std::copy(stage.begin(),stage.end(),output+(i+1)*dim);++*written;
        }
        return 0;
    } catch(...){return 3;}
}
