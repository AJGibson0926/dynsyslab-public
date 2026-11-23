#include <dynsyslab/numerics/NonlinearLeastSquares.hpp>
#include <dynsyslab/numerics/LeastSquares.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
namespace dynsyslab::numerics {
NonlinearLeastSquaresResult boundedNonlinearLeastSquares(
    const std::function<std::vector<double>(const std::vector<double>&)>& function,
    const std::vector<double>& initial,const std::vector<double>& lower,
    const std::vector<double>& upper,const NonlinearLeastSquaresOptions& o){
    const auto n=initial.size();if(!n||n>1024||lower.size()!=n||upper.size()!=n||!o.maximumIterations||
        !std::isfinite(o.finiteDifferenceStep)||o.finiteDifferenceStep<=0||o.finiteDifferenceStep>=1||
        !std::isfinite(o.stepTolerance)||o.stepTolerance<=0||!std::isfinite(o.gradientTolerance)||o.gradientTolerance<=0)
        throw std::invalid_argument("Invalid nonlinear least-squares settings");
    std::vector<double> z(n),scale(n);NonlinearLeastSquaresResult out;
    for(std::size_t j=0;j<n;++j){scale[j]=upper[j]-lower[j];
        if(!std::isfinite(initial[j])||!std::isfinite(lower[j])||!std::isfinite(upper[j])||!std::isfinite(scale[j])||scale[j]<=0||initial[j]<lower[j]||initial[j]>upper[j])throw std::invalid_argument("Invalid least-squares box/initial point");
        z[j]=(initial[j]-lower[j])/scale[j];}
    auto physical=[&](const auto& q){auto x=q;for(std::size_t j=0;j<n;++j)x[j]=lower[j]+scale[j]*q[j];return x;};
    std::size_t m=0;
    auto evaluate=[&](const auto& q){auto r=function(physical(q));++out.evaluations;
        if(r.empty()||(m&&r.size()!=m))throw std::invalid_argument("Nonlinear residual changed dimension");
        for(double x:r)if(!std::isfinite(x))throw std::domain_error("Nonfinite nonlinear residual");return r;};
    auto norm=[](const auto& v){double s=0;for(double x:v)s+=x*x;if(!std::isfinite(s))throw std::overflow_error("Nonlinear objective overflow");return s;};
    auto residual=evaluate(z);m=residual.size();double cost=norm(residual),lambda=1e-3;
    for(std::size_t iteration=0;iteration<o.maximumIterations;++iteration){out.iterations=iteration+1;
        std::vector<double> jac(m*n),gradient(n);
        for(std::size_t j=0;j<n;++j){auto a=z,b=z;a[j]=std::max(0.,z[j]-o.finiteDifferenceStep);b[j]=std::min(1.,z[j]+o.finiteDifferenceStep);
            const auto ra=evaluate(a),rb=evaluate(b);for(std::size_t i=0;i<m;++i){jac[i*n+j]=(rb[i]-ra[i])/(b[j]-a[j]);gradient[j]+=jac[i*n+j]*residual[i];}}
        double pg=0;for(std::size_t j=0;j<n;++j){double g=gradient[j];if((z[j]==0&&g>0)||(z[j]==1&&g<0))g=0;pg=std::max(pg,std::abs(g));}
        if(pg<=o.gradientTolerance){out.converged=true;out.reason="projected_gradient";break;}
        bool accepted=false;
        for(int attempt=0;attempt<20;++attempt){std::vector<double> matrix((m+n)*n),rhs(m+n);
            std::copy(jac.begin(),jac.end(),matrix.begin());for(std::size_t i=0;i<m;++i)rhs[i]=-residual[i];
            for(std::size_t j=0;j<n;++j){double d=0;for(std::size_t i=0;i<m;++i)d+=jac[i*n+j]*jac[i*n+j];matrix[(m+j)*n+j]=std::sqrt(lambda*std::max(1.,d));}
            const auto step=solveLeastSquares(matrix,m+n,n,rhs,1);auto candidate=z;double length=0;
            for(std::size_t j=0;j<n;++j){candidate[j]=std::clamp(z[j]+step[j],0.,1.);length=std::max(length,std::abs(candidate[j]-z[j]));}
            auto next=evaluate(candidate);const double nextCost=norm(next);
            if(nextCost<cost){z=std::move(candidate);residual=std::move(next);cost=nextCost;lambda=std::max(1e-12,lambda/3);accepted=true;
                if(length<o.stepTolerance){out.converged=true;out.reason="accepted_step";}break;}
            lambda*=10;
        }
        if(out.converged)break;
        if(!accepted){out.reason="no_descent";break;}
    }
    if(out.reason.empty())out.reason="iteration_limit";
    out.parameters=physical(z);out.squaredResidual=cost;return out;
}
}
