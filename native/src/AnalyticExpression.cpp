#include <dynsyslab/core/AnalyticExpression.hpp>
#include <dynsyslab/core/EllipticIntegrals.hpp>
#include <dynsyslab/core/PolynomialRoots.hpp>
#include "ExpressionFunctions.hpp"
#include <cstdint>
#include <limits>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <utility>

namespace dynsyslab::core {
namespace {
bool unaryName(const std::string& op) {
    return op=="negative"||op=="exp"||op=="log"||op=="sin"||op=="cos"||op=="tanh"||op=="sech"||op=="ellipk"||op=="ellipe";
}
bool functionName(const std::string& op) { return op=="atan2"||op=="ellipf"||op=="ellipeinc"; }
bool callName(const std::string& op) { return op=="ellippi"||op=="polyroot"; }
AnalyticExpression power_(const AnalyticExpression& x,double p) { return x.power(p); }
}  // namespace
struct AnalyticExpression::Node {
    std::string operation;
    std::size_t dimension{}, index{};
    double scalar{};
    std::vector<AnalyticExpression> arguments;
    std::vector<double> parameters;  // a special function's constants; a sum's bounds
    std::uint64_t free{};            // the summation variables this node depends on, one bit per slot
    bool stateful{};                 // depends on the state
};
AnalyticExpression::AnalyticExpression(std::shared_ptr<const Node> node) : node_(std::move(node)) {}
std::shared_ptr<const AnalyticExpression::Node> AnalyticExpression::make(Node n) {
    for (const auto& a : n.arguments) { n.free |= a.node_->free; n.stateful = n.stateful || a.node_->stateful; }
    if (n.operation == "state") n.stateful = true;
    if (n.operation == "index") n.free |= std::uint64_t(1) << n.index;
    if (n.operation == "sum") n.free &= ~(std::uint64_t(1) << n.index);
    return std::make_shared<const Node>(std::move(n));
}
AnalyticExpression AnalyticExpression::constant(double value, std::size_t dimension) {
    if (!dimension || !std::isfinite(value)) throw std::invalid_argument("invalid expression constant");
    return AnalyticExpression(make(Node{"constant",dimension,0,value,{}}));
}
AnalyticExpression AnalyticExpression::state(std::size_t dimension, std::size_t index) {
    if (!dimension || index>=dimension) throw std::invalid_argument("invalid expression state axis");
    return AnalyticExpression(make(Node{"state",dimension,index,0,{}}));
}
std::size_t AnalyticExpression::stateDimension() const noexcept { return node_->dimension; }
double AnalyticExpression::evaluate(const double* x) const {
    std::unordered_map<const Node*, double> memo;
    std::vector<double> indices(64, std::numeric_limits<double>::quiet_NaN());
    return evaluate(x, memo, indices);
}
double AnalyticExpression::evaluate(const double* x, std::unordered_map<const Node*, double>& memo, std::vector<double>& indices) const {
    const auto& n=*node_; const auto& op=n.operation;
    if(op=="constant")return n.scalar;
    if(op=="state")return x[n.index];
    if(op=="index")return indices[n.index];
    const bool cached=n.free==0;  // a node that depends on a summation variable changes with it
    if(cached)if(const auto found=memo.find(&n);found!=memo.end())return found->second;
    const auto value=[&]() -> double {
    if(op=="sum") {
        const double saved=indices[n.index], from=n.parameters[0], to=n.parameters[1];
        double total=0, carry=0;
        int quiet=0;
        double result=std::numeric_limits<double>::quiet_NaN();
        const double limit=std::isinf(to)?from+1e6:to;
        for(double k=from;k<=limit;++k) {
            indices[n.index]=k;
            const double term=n.arguments[0].evaluate(x,memo,indices);
            if(!std::isfinite(term)){result=std::numeric_limits<double>::quiet_NaN();break;}
            const double y=term-carry,t=total+y;carry=(t-total)-y;total=t;  // compensated summation
            if(!std::isinf(to)){result=total;continue;}
            quiet=std::abs(term)<=1e-17*std::abs(total)+1e-300?quiet+1:0;
            if(quiet>=5){result=total;break;}
        }
        if(!std::isinf(to)&&from>to)result=0;
        indices[n.index]=saved;
        return result;
    }
    if(const auto* f=detail::specialFunction(op)) {  // the leading arguments are its parameters
        std::vector<double> p,a;
        for(std::size_t i=0;i<n.arguments.size();++i)(i<f->parameters?p:a).push_back(n.arguments[i].evaluate(x,memo,indices));
        return f->value(p,a);
    }
    if(op=="ellippi")return ellipticPi(n.arguments[0].evaluate(x,memo,indices),n.arguments[1].evaluate(x,memo,indices),n.arguments[2].evaluate(x,memo,indices));
    if(op=="polyroot") {
        std::vector<double> c;
        for(const auto& item:n.arguments)c.push_back(item.evaluate(x,memo,indices));
        return realPolynomialRoot(n.index,c);
    }
    const double a=n.arguments[0].evaluate(x,memo,indices);
    if(op=="negative")return -a;
    if(op=="power")return std::pow(a,n.scalar);
    if(op=="exp")return std::exp(a);
    if(op=="log")return std::log(a);
    if(op=="sin")return std::sin(a);
    if(op=="cos")return std::cos(a);
    if(op=="tanh")return std::tanh(a);
    if(op=="sech") {const double e=std::exp(-std::abs(a));return 2*e/(1+e*e);}
    if(op=="ellipk")return ellipticK(a);
    if(op=="ellipe")return ellipticE(a);
    const double b=n.arguments[1].evaluate(x,memo,indices);
    if(op=="add")return a+b;
    if(op=="subtract")return a-b;
    if(op=="multiply")return a*b;
    if(op=="divide")return a/b;
    if(op=="atan2")return std::atan2(a,b);
    if(op=="ellipf")return ellipticF(a,b);
    if(op=="ellipeinc")return ellipticEIncomplete(a,b);
    throw std::logic_error("unknown expression operation");
    }();
    if(cached)memo.emplace(&n,value);
    return value;
}
AnalyticExpression AnalyticExpression::binary(const std::string& op,const AnalyticExpression& b) const {
    if(stateDimension()!=b.stateDimension())throw std::invalid_argument("expression dimensions differ");
    const bool ac=node_->operation=="constant",bc=b.node_->operation=="constant";
    if(op=="add") { if(ac&&node_->scalar==0)return b; if(bc&&b.node_->scalar==0)return *this; }
    if(op=="subtract"&&bc&&b.node_->scalar==0)return *this;
    if(op=="multiply") {
        if((ac&&node_->scalar==0)||(bc&&b.node_->scalar==0))return constant(0,stateDimension());
        if(ac&&node_->scalar==1)return b;
        if(bc&&b.node_->scalar==1)return *this;
    }
    if(op=="divide"&&bc&&b.node_->scalar==1)return *this;
    auto value=AnalyticExpression(make(Node{op,stateDimension(),0,0,{*this,b}}));
    if(ac&&bc){const double v=value.evaluate(nullptr);if(std::isfinite(v))return constant(v,stateDimension());}  // folding is an optimization: a non-finite value stays a node
    return value;
}
AnalyticExpression AnalyticExpression::power(double value) const {
    if(!std::isfinite(value))throw std::invalid_argument("expression exponent must be finite");
    if(value==0)return constant(1,stateDimension());
    if(value==1)return *this;
    if(node_->operation=="constant"&&std::isfinite(std::pow(node_->scalar,value)))return constant(std::pow(node_->scalar,value),stateDimension());
    return AnalyticExpression(make(Node{"power",stateDimension(),0,value,{*this}}));
}
AnalyticExpression AnalyticExpression::unary(const std::string& op) const {
    if(!unaryName(op))
        throw std::invalid_argument("unsupported unary expression operation");
    auto result=AnalyticExpression(make(Node{op,stateDimension(),0,0,{*this}}));
    if(node_->operation=="constant"&&std::isfinite(result.evaluate(nullptr)))return constant(result.evaluate(nullptr),stateDimension());
    return result;
}
AnalyticExpression AnalyticExpression::function(const std::string& op,const AnalyticExpression& b) const {
    if(!functionName(op))throw std::invalid_argument("unsupported two-argument expression function");
    return binary(op,b);
}
AnalyticExpression AnalyticExpression::call(const std::string& op,std::vector<AnalyticExpression> args) {
    if(!callName(op)||args.empty())throw std::invalid_argument("unsupported expression function "+op);
    const auto dimension=args.front().stateDimension();
    for(const auto& item:args)if(item.stateDimension()!=dimension)throw std::invalid_argument("expression dimensions differ");
    std::size_t index=0;
    if(op=="ellippi"&&args.size()!=3)throw std::invalid_argument("ellippi(phi, n, m) takes three arguments");
    if(op=="polyroot") {
        const auto& k=*args.front().node_;
        if(k.operation!="constant"||k.scalar<0||k.scalar!=std::floor(k.scalar)||args.size()<3)
            throw std::invalid_argument("polyroot(k, c_0, ..., c_d) needs a constant root index and a polynomial of degree >= 1");
        index=std::size_t(k.scalar);
        args.erase(args.begin());
    }
    auto value=AnalyticExpression(make(Node{op,dimension,index,0,std::move(args)}));
    bool constantArguments=true;
    for(const auto& item:value.node_->arguments)constantArguments=constantArguments&&item.node_->operation=="constant";
    if(constantArguments&&std::isfinite(value.evaluate(nullptr)))return constant(value.evaluate(nullptr),dimension);
    return value;
}
bool AnalyticExpression::dependsOnState() const noexcept { return node_->stateful; }
AnalyticExpression AnalyticExpression::special(const std::string& name,std::vector<double> parameters,std::vector<AnalyticExpression> args) {
    if(args.empty())throw std::invalid_argument(name+" needs an argument");
    std::vector<AnalyticExpression> p;
    for(const double v:parameters)p.push_back(constant(v,args.front().stateDimension()));
    return special(name,std::move(p),std::move(args));
}
AnalyticExpression AnalyticExpression::special(const std::string& name,std::vector<AnalyticExpression> parameters,std::vector<AnalyticExpression> args) {
    const auto* f=detail::specialFunction(name);
    if(!f)throw std::invalid_argument("unknown special function "+name);
    if(parameters.size()!=f->parameters||args.size()!=f->arguments)
        throw std::invalid_argument(name+" takes "+std::to_string(f->parameters)+" parameter(s) and "+std::to_string(f->arguments)+" argument(s)");
    for(const auto& p:parameters)if(p.dependsOnState())throw std::invalid_argument(name+" parameters may not depend on the state");
    const auto dimension=args.front().stateDimension();
    parameters.insert(parameters.end(),args.begin(),args.end());
    for(const auto& item:parameters)if(item.stateDimension()!=dimension)throw std::invalid_argument("expression dimensions differ");
    bool constantArguments=true;
    for(const auto& item:parameters)constantArguments=constantArguments&&item.node_->operation=="constant";
    auto value=AnalyticExpression(make(Node{name,dimension,0,0,std::move(parameters)}));
    if(constantArguments) {
        const double v=value.evaluate(nullptr);
        if(std::isfinite(v))return constant(v,dimension);
    }
    return value;
}
AnalyticExpression AnalyticExpression::index(std::size_t dimension,std::size_t slot) {
    if(!dimension||slot>=64)throw std::invalid_argument("a summation variable needs a slot from 0 to 63");
    return AnalyticExpression(make(Node{"index",dimension,slot,0,{}}));
}
AnalyticExpression AnalyticExpression::sum(std::size_t slot,double from,double to,const AnalyticExpression& term) {
    if(slot>=64||!std::isfinite(from)||from!=std::floor(from)||std::isnan(to)||(std::isfinite(to)&&to!=std::floor(to))||to==-INFINITY)
        throw std::invalid_argument("a sum runs over integers from a finite start to an integer end or infinity");
    return AnalyticExpression(make(Node{"sum",term.stateDimension(),slot,0,{term},{from,to}}));
}
AnalyticExpression AnalyticExpression::derivative(std::size_t axis) const {
    if(axis>=stateDimension())throw std::out_of_range("expression derivative axis");
    std::unordered_map<const Node*, AnalyticExpression> memo;
    return derivative(axis, memo);
}
AnalyticExpression AnalyticExpression::derivative(std::size_t axis, std::unordered_map<const Node*, AnalyticExpression>& memo) const {
    const auto& n=*node_; const auto& op=n.operation;
    if(op=="constant")return constant(0,stateDimension());
    if(op=="state")return constant(axis==n.index?1:0,stateDimension());
    if(op=="index")return constant(0,stateDimension());
    if(const auto found=memo.find(&n);found!=memo.end())return found->second;
    const auto result=[&]() -> AnalyticExpression {
    if(op=="sum")return sum(n.index,n.parameters[0],n.parameters[1],n.arguments[0].derivative(axis,memo));
    if(const auto* f=detail::specialFunction(op)) {  // parameters do not depend on the state
        const std::vector<AnalyticExpression> p(n.arguments.begin(),n.arguments.begin()+std::ptrdiff_t(f->parameters));
        const std::vector<AnalyticExpression> a(n.arguments.begin()+std::ptrdiff_t(f->parameters),n.arguments.end());
        auto total=constant(0,stateDimension());
        for(std::size_t i=0;i<a.size();++i) {
            const auto d=a[i].derivative(axis,memo);
            if(d.node_->operation=="constant"&&d.node_->scalar==0)continue;
            total=total+f->partial(p,a,*this,i)*d;
        }
        return total;
    }
    if(op=="ellippi") {
        // DLMF 19.4 in the parameter m: dPi/dphi = 1/(Delta (1 - n sin^2 phi)), and
        // dPi/dn = [E + (m-n)F/n + (n^2-m)Pi/n - n Delta sin cos/(1 - n sin^2)] / (2(m-n)(n-1)),
        // dPi/dm = [E/(m-1) + Pi - m sin cos/((m-1) Delta)] / (2(n-m)).
        const auto& phi=n.arguments[0]; const auto& nn=n.arguments[1]; const auto& m=n.arguments[2];
        const auto s=phi.unary("sin"), c=phi.unary("cos"), delta=(1-m*s.power(2)).power(0.5), d=1-nn*s.power(2);
        const auto dphi=phi.derivative(axis,memo), dn=nn.derivative(axis,memo), dm=m.derivative(axis,memo);
        auto result=dphi/(delta*d);
        const auto F=phi.function("ellipf",m), E=phi.function("ellipeinc",m);
        if(dn.node_->operation!="constant"||dn.node_->scalar!=0)
            result=result+(E+(m-nn)*F/nn+(nn*nn-m)**this/nn-nn*delta*s*c/d)/(2*(m-nn)*(nn-1))*dn;
        if(dm.node_->operation!="constant"||dm.node_->scalar!=0)
            result=result+(E/(m-1)+*this-m*s*c/((m-1)*delta))/(2*(nn-m))*dm;
        return result;
    }
    if(op=="polyroot") {
        // x = root of P(x) = sum_i c_i x^(d-i): dx = -sum_i x^(d-i) dc_i / P'(x).
        const std::size_t degree=n.arguments.size()-1;
        auto slope=constant(0,stateDimension()), change=constant(0,stateDimension());
        for(std::size_t i=0;i<=degree;++i) {
            const double power=double(degree-i);
            if(i<degree)slope=slope+power*n.arguments[i]*power_(*this,power-1);
            change=change+power_(*this,power)*n.arguments[i].derivative(axis,memo);
        }
        return -change/slope;
    }
    const auto& a=n.arguments[0]; const auto da=a.derivative(axis,memo);
    if(op=="negative")return -da;
    if(op=="power")return n.scalar*a.power(n.scalar-1)*da;
    if(op=="exp")return a.unary("exp")*da;
    if(op=="log")return da/a;
    if(op=="sin")return a.unary("cos")*da;
    if(op=="cos")return -a.unary("sin")*da;
    if(op=="tanh")return (1-a.unary("tanh").power(2))*da;
    if(op=="sech")return -a.unary("sech")*a.unary("tanh")*da;
    // dK/dm = (E - (1-m)K)/(2m(1-m)), dE/dm = (E - K)/(2m) (DLMF 19.4.1, with m = k^2).
    if(op=="ellipk")return (a.unary("ellipe")-(1-a)*a.unary("ellipk"))/(2*a*(1-a))*da;
    if(op=="ellipe")return (a.unary("ellipe")-a.unary("ellipk"))/(2*a)*da;
    const auto& b=n.arguments[1]; const auto db=b.derivative(axis,memo);
    if(op=="add")return da+db;
    if(op=="subtract")return da-db;
    if(op=="multiply")return da*b+a*db;
    if(op=="divide")return (da*b-a*db)/b.power(2);
    if(op=="atan2")return (b*da-a*db)/(a.power(2)+b.power(2));
    // DLMF 19.4.2 in the parameter m: dF/dphi = 1/Delta, dE/dphi = Delta, Delta = sqrt(1 - m sin^2 phi).
    const auto delta=(1-b*a.unary("sin").power(2)).power(0.5);
    if(op=="ellipf")
        return da/delta+(a.function("ellipeinc",b)/(2*b*(1-b))-a.function("ellipf",b)/(2*b)
                         -(2*a).unary("sin")/(4*(1-b)*delta))*db;
    if(op=="ellipeinc")return delta*da+(a.function("ellipeinc",b)-a.function("ellipf",b))/(2*b)*db;
    throw std::logic_error("unknown expression derivative");
    }();
    memo.emplace(&n,result);
    return result;
}
AnalyticExpression AnalyticExpression::flowDerivative(
    const std::vector<AnalyticExpression>& f,std::size_t order) const {
    if(f.size()!=stateDimension())throw std::invalid_argument("expression vector field dimension mismatch");
    for(const auto& item:f)if(item.stateDimension()!=stateDimension())
        throw std::invalid_argument("vector field expression dimension mismatch");
    auto result=*this;
    for(std::size_t k=0;k<order;++k) {
        auto next=constant(0,stateDimension());
        for(std::size_t axis=0;axis<stateDimension();++axis)next=next+result.derivative(axis)*f[axis];
        result=std::move(next);
    }
    return result;
}
json::Value AnalyticExpression::toJson() const {
    // Nodes in post order, each once, with the number of references to each.
    std::unordered_map<const Node*,std::size_t> uses;std::vector<const AnalyticExpression*> order;
    std::function<void(const AnalyticExpression&)> visit=[&](const AnalyticExpression& e) {
        if(uses[e.node_.get()]++)return;
        for(const auto& item:e.node_->arguments)visit(item);
        order.push_back(&e);
    };
    visit(*this);
    bool shared=false;
    for(const auto& [node,count]:uses)shared=shared||(count>1&&node->operation!="constant"&&node->operation!="state");
    const auto object=[](const Node& n) {
        json::Value::Object value{{"operation",n.operation},{"state_dimension",double(n.dimension)}};
        if(n.operation=="state")value["state_index"]=double(n.index);
        if(n.operation=="constant"||n.operation=="power")value["scalar"]=n.scalar;
        if(n.operation=="polyroot")value["root_index"]=double(n.index);
        if(n.operation=="index"||n.operation=="sum")value["slot"]=double(n.index);
        if(n.operation=="sum") {
            value["from"]=n.parameters[0];
            value["to"]=std::isinf(n.parameters[1])?json::Value():json::Value(n.parameters[1]);  // null: infinity
        }
        return value;
    };
    if(!shared) {  // a tree: the nested v1 encoding
        std::unordered_map<const Node*,json::Value> done;
        for(const auto* e:order) {
            auto value=object(*e->node_);
            value["encoding"]="analytic-state-expression/v1";
            if(!e->node_->arguments.empty()) {
                json::Value::Array args;
                for(const auto& item:e->node_->arguments)args.push_back(done.at(item.node_.get()));
                value["arguments"]=std::move(args);
            }
            done.emplace(e->node_.get(),json::Value(std::move(value)));
        }
        return done.at(node_.get());
    }
    std::unordered_map<const Node*,std::size_t> position;json::Value::Array nodes;
    for(const auto* e:order) {
        auto value=object(*e->node_);
        if(!e->node_->arguments.empty()) {
            json::Value::Array args;
            for(const auto& item:e->node_->arguments)args.emplace_back(double(position.at(item.node_.get())));
            value["arguments"]=std::move(args);
        }
        position.emplace(e->node_.get(),nodes.size());nodes.emplace_back(std::move(value));
    }
    const double root=double(nodes.size()-1);
    return json::Value::Object{{"encoding","analytic-state-expression/v2"},{"state_dimension",double(stateDimension())},
                               {"nodes",json::Value(std::move(nodes))},{"root",root}};
}
AnalyticExpression AnalyticExpression::rebuild(const json::Value& entry,std::size_t dimension,std::vector<AnalyticExpression> args) {
    const auto op=entry.at("operation").string();
    for(const auto& item:args)if(item.stateDimension()!=dimension)throw std::invalid_argument("expression dimensions differ");
    if(op=="constant")return constant(entry.at("scalar").number(),dimension);
    if(op=="state") {
        const double index=entry.at("state_index").number(-1);
        if(index<0||index!=std::floor(index)||index>=double(dimension))throw std::invalid_argument("invalid expression axis");
        return state(dimension,std::size_t(index));
    }
    if(op=="index")return index(dimension,std::size_t(entry.at("slot").number()));
    if(op=="sum") {
        if(args.size()!=1)throw std::invalid_argument("expression arity mismatch");
        return sum(std::size_t(entry.at("slot").number()),entry.at("from").number(),
                   entry.at("to").isNull()?INFINITY:entry.at("to").number(),args[0]);
    }
    if(const auto* f=detail::specialFunction(op)) {  // parameters are its leading arguments
        if(args.size()!=f->parameters+f->arguments)throw std::invalid_argument("expression arity mismatch");
        std::vector<AnalyticExpression> p(args.begin(),args.begin()+std::ptrdiff_t(f->parameters));
        args.erase(args.begin(),args.begin()+std::ptrdiff_t(f->parameters));
        return special(op,std::move(p),std::move(args));
    }
    if(callName(op)) {
        if(op=="polyroot")args.insert(args.begin(),constant(entry.at("root_index").number(),dimension));
        return call(op,std::move(args));
    }
    const bool binaryOp=op=="add"||op=="subtract"||op=="multiply"||op=="divide"||functionName(op);
    if(args.size()!=(binaryOp?2u:1u))throw std::invalid_argument("expression arity mismatch");
    if(functionName(op))return args[0].function(op,args[1]);
    if(binaryOp)return args[0].binary(op,args[1]);
    if(op=="power")return args[0].power(entry.at("scalar").number());
    return args[0].unary(op);
}
AnalyticExpression AnalyticExpression::fromJson(const json::Value& value) {
    const auto encoding=value.at("encoding").string();
    if(encoding!="analytic-state-expression/v1"&&encoding!="analytic-state-expression/v2")
        throw std::invalid_argument("unsupported analytic expression encoding");
    const double d=value.at("state_dimension").number();
    if(d<1||d>100000||d!=std::floor(d))throw std::invalid_argument("invalid expression dimension");
    const auto dimension=std::size_t(d);
    if(encoding=="analytic-state-expression/v1") {  // a tree
        std::vector<AnalyticExpression> args;
        if(!value.at("arguments").isNull())for(const auto& item:value.at("arguments").array())args.push_back(fromJson(item));
        return rebuild(value,dimension,std::move(args));
    }
    // A node table: each node's arguments are indices of earlier nodes, so each is built once.
    std::vector<AnalyticExpression> built;
    for(const auto& entry:value.at("nodes").array()) {
        std::vector<AnalyticExpression> args;
        if(!entry.at("arguments").isNull())
            for(const auto& index:entry.at("arguments").array()) {
                const double k=index.number(-1);
                if(k<0||k!=std::floor(k)||k>=double(built.size()))throw std::invalid_argument("expression node table refers forward");
                args.push_back(built[std::size_t(k)]);
            }
        built.push_back(rebuild(entry,dimension,std::move(args)));
    }
    const double root=value.at("root").number(-1);
    if(root<0||root!=std::floor(root)||root>=double(built.size()))throw std::invalid_argument("invalid expression root");
    return built[std::size_t(root)];
}
AnalyticExpression AnalyticExpression::operator+(const AnalyticExpression& b) const{return binary("add",b);}
AnalyticExpression AnalyticExpression::operator-(const AnalyticExpression& b) const{return binary("subtract",b);}
AnalyticExpression AnalyticExpression::operator*(const AnalyticExpression& b) const{return binary("multiply",b);}
AnalyticExpression AnalyticExpression::operator/(const AnalyticExpression& b) const{return binary("divide",b);}
AnalyticExpression AnalyticExpression::operator-() const{return unary("negative");}
AnalyticExpression AnalyticExpression::operator+(double b) const{return *this+constant(b,stateDimension());}
AnalyticExpression AnalyticExpression::operator-(double b) const{return *this-constant(b,stateDimension());}
AnalyticExpression AnalyticExpression::operator*(double b) const{return *this*constant(b,stateDimension());}
AnalyticExpression AnalyticExpression::operator/(double b) const{return *this/constant(b,stateDimension());}
AnalyticExpression operator+(double a,const AnalyticExpression& b){return b+a;}
AnalyticExpression operator-(double a,const AnalyticExpression& b){return AnalyticExpression::constant(a,b.stateDimension())-b;}
AnalyticExpression operator*(double a,const AnalyticExpression& b){return b*a;}
AnalyticExpression operator/(double a,const AnalyticExpression& b){return AnalyticExpression::constant(a,b.stateDimension())/b;}
}  // namespace dynsyslab::core
