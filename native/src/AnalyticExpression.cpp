#include <dynsyslab/core/AnalyticExpression.hpp>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace dynsyslab::core {
struct AnalyticExpression::Node {
    std::string operation;
    std::size_t dimension{}, index{};
    double scalar{};
    std::vector<AnalyticExpression> arguments;
};
AnalyticExpression::AnalyticExpression(std::shared_ptr<const Node> node) : node_(std::move(node)) {}
AnalyticExpression AnalyticExpression::constant(double value, std::size_t dimension) {
    if (!dimension || !std::isfinite(value)) throw std::invalid_argument("invalid expression constant");
    return AnalyticExpression(std::make_shared<Node>(Node{"constant",dimension,0,value,{}}));
}
AnalyticExpression AnalyticExpression::state(std::size_t dimension, std::size_t index) {
    if (!dimension || index>=dimension) throw std::invalid_argument("invalid expression state axis");
    return AnalyticExpression(std::make_shared<Node>(Node{"state",dimension,index,0,{}}));
}
std::size_t AnalyticExpression::stateDimension() const noexcept { return node_->dimension; }
double AnalyticExpression::evaluate(const double* x) const {
    const auto& n=*node_; const auto& op=n.operation;
    if(op=="constant")return n.scalar;
    if(op=="state")return x[n.index];
    const double a=n.arguments[0].evaluate(x);
    if(op=="negative")return -a;
    if(op=="power")return std::pow(a,n.scalar);
    if(op=="exp")return std::exp(a);
    if(op=="log")return std::log(a);
    if(op=="sin")return std::sin(a);
    if(op=="cos")return std::cos(a);
    if(op=="tanh")return std::tanh(a);
    if(op=="sech") {const double e=std::exp(-std::abs(a));return 2*e/(1+e*e);}
    const double b=n.arguments[1].evaluate(x);
    if(op=="add")return a+b;
    if(op=="subtract")return a-b;
    if(op=="multiply")return a*b;
    if(op=="divide")return a/b;
    throw std::logic_error("unknown expression operation");
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
    auto value=AnalyticExpression(std::make_shared<Node>(Node{op,stateDimension(),0,0,{*this,b}}));
    if(ac&&bc)return constant(value.evaluate(nullptr),stateDimension());
    return value;
}
AnalyticExpression AnalyticExpression::power(double value) const {
    if(!std::isfinite(value))throw std::invalid_argument("expression exponent must be finite");
    if(value==0)return constant(1,stateDimension());
    if(value==1)return *this;
    if(node_->operation=="constant")return constant(std::pow(node_->scalar,value),stateDimension());
    return AnalyticExpression(std::make_shared<Node>(Node{"power",stateDimension(),0,value,{*this}}));
}
AnalyticExpression AnalyticExpression::unary(const std::string& op) const {
    if(op!="negative"&&op!="exp"&&op!="log"&&op!="sin"&&op!="cos"&&op!="tanh"&&op!="sech")
        throw std::invalid_argument("unsupported unary expression operation");
    auto result=AnalyticExpression(std::make_shared<Node>(Node{op,stateDimension(),0,0,{*this}}));
    if(node_->operation=="constant")return constant(result.evaluate(nullptr),stateDimension());
    return result;
}
AnalyticExpression AnalyticExpression::derivative(std::size_t axis) const {
    if(axis>=stateDimension())throw std::out_of_range("expression derivative axis");
    const auto& n=*node_; const auto& op=n.operation;
    if(op=="constant")return constant(0,stateDimension());
    if(op=="state")return constant(axis==n.index?1:0,stateDimension());
    const auto& a=n.arguments[0]; const auto da=a.derivative(axis);
    if(op=="negative")return -da;
    if(op=="power")return n.scalar*a.power(n.scalar-1)*da;
    if(op=="exp")return a.unary("exp")*da;
    if(op=="log")return da/a;
    if(op=="sin")return a.unary("cos")*da;
    if(op=="cos")return -a.unary("sin")*da;
    if(op=="tanh")return (1-a.unary("tanh").power(2))*da;
    if(op=="sech")return -a.unary("sech")*a.unary("tanh")*da;
    const auto& b=n.arguments[1]; const auto db=b.derivative(axis);
    if(op=="add")return da+db;
    if(op=="subtract")return da-db;
    if(op=="multiply")return da*b+a*db;
    if(op=="divide")return (da*b-a*db)/b.power(2);
    throw std::logic_error("unknown expression derivative");
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
    const auto& n=*node_;
    json::Value::Object value{{"encoding","analytic-state-expression/v1"},
        {"operation",n.operation},{"state_dimension",double(n.dimension)}};
    if(n.operation=="state")value["state_index"]=double(n.index);
    if(n.operation=="constant"||n.operation=="power")value["scalar"]=n.scalar;
    if(!n.arguments.empty()) {
        json::Value::Array args;
        for(const auto& item:n.arguments)args.push_back(item.toJson());
        value["arguments"]=std::move(args);
    }
    return value;
}
AnalyticExpression AnalyticExpression::fromJson(const json::Value& value) {
    if(value.at("encoding").string()!="analytic-state-expression/v1")
        throw std::invalid_argument("unsupported analytic expression encoding");
    const double d=value.at("state_dimension").number();
    if(d<1||d>100000||d!=std::floor(d))throw std::invalid_argument("invalid expression dimension");
    const auto dimension=std::size_t(d); const auto op=value.at("operation").string();
    if(op=="constant")return constant(value.at("scalar").number(),dimension);
    if(op=="state") {
        const double index=value.at("state_index").number(-1);
        if(index<0||index!=std::floor(index)||index>=d)throw std::invalid_argument("invalid expression axis");
        return state(dimension,std::size_t(index));
    }
    const auto& args=value.at("arguments").array();
    const bool binaryOp=op=="add"||op=="subtract"||op=="multiply"||op=="divide";
    if(args.size()!=(binaryOp?2:1))throw std::invalid_argument("expression arity mismatch");
    const auto a=fromJson(args[0]);
    if(a.stateDimension()!=dimension)throw std::invalid_argument("expression dimensions differ");
    if(binaryOp)return a.binary(op,fromJson(args[1]));
    if(op=="power")return a.power(value.at("scalar").number());
    return a.unary(op);
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
