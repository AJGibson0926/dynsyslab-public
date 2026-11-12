#include "JetDecision.hpp"

#include <cstdlib>
#include <filesystem>
#include <memory>
#include <sstream>
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

// Opaque types of the optional Z3 C ABI. The library location comes from
// process configuration only. A query can never select code to load.
struct _Z3_config;
struct _Z3_context;
namespace dynsyslab::core::jet {
namespace {
std::string product(const std::vector<std::string>& factors) {
    if (factors.empty()) return "1";
    if (factors.size()==1) return factors.front();
    std::string result="(*";
    for (const auto& f:factors) result+=" "+f;
    return result+")";
}
std::string monomial(const Powers& powers,char state) {
    std::vector<std::string> factors;
    for (std::size_t i=0;i<powers.size();++i)
        for (int j=0;j<powers[i];++j) factors.push_back(state+std::to_string(i));
    return product(factors);
}
// Compare numerators after multiplying by both nonzero source denominators.
// All obligations are quantifier-free real polynomial feasibility problems.
std::string numeratorText(const Polynomial& p,const Powers& denominatorPowers,char state) {
    if (p.empty()) return "0";
    std::vector<std::string> terms;
    for (const auto& [powers,c]:p) {
        const Rational positive=abs(c);
        const std::string scalar="(/ "+numerator(positive).str()+" "+denominator(positive).str()+")";
        Powers exponents=denominatorPowers;
        for (std::size_t i=0;i<powers.size();++i) exponents[i]+=powers[i];
        auto term=product({scalar,monomial(exponents,state)});
        terms.push_back(c<0?"(- "+term+")":term);
    }
    if (terms.size()==1) return terms[0];
    std::string sum="(+";
    for (const auto& term:terms) sum+=" "+term;
    return sum+")";
}
std::string equalValues(const Polynomial& p,std::size_t dimension) {
    if (p.empty()) return "true";
    Powers denominatorPowers(dimension);
    for (const auto& term:p) for (std::size_t i=0;i<dimension;++i)
        denominatorPowers[i]=std::max(denominatorPowers[i],-term.first[i]);
    return "(= "+product({numeratorText(p,denominatorPowers,'a'),monomial(denominatorPowers,'b')})+" "+
        product({numeratorText(p,denominatorPowers,'b'),monomial(denominatorPowers,'a')})+")";
}
class Library {
public:
    explicit Library(const std::filesystem::path& path) {
        if (!path.is_absolute()) throw std::runtime_error("DYNSYSLAB_Z3_LIBRARY must be absolute");
#ifdef _WIN32
        handle=LoadLibraryExW(path.c_str(),nullptr,LOAD_LIBRARY_SEARCH_DLL_LOAD_DIR|LOAD_LIBRARY_SEARCH_DEFAULT_DIRS);
#else
        handle=dlopen(path.c_str(),RTLD_NOW|RTLD_LOCAL);
#endif
        if (!handle) throw std::runtime_error("cannot load configured Z3 library");
    }
    ~Library() {
#ifdef _WIN32
        FreeLibrary(handle);
#else
        dlclose(handle);
#endif
    }
    Library(const Library&)=delete;
    Library& operator=(const Library&)=delete;
    template<class Function> Function symbol(const char* name) const {
#ifdef _WIN32
        auto address=GetProcAddress(handle,name);
#else
        auto address=dlsym(handle,name);
#endif
        if (!address) throw std::runtime_error(std::string("missing Z3 API: ")+name);
        return reinterpret_cast<Function>(address);
    }
private:
#ifdef _WIN32
    HMODULE handle{};
#else
    void* handle{};
#endif
};
} // namespace

std::string fiberProblem(const std::vector<Polynomial>& jet,const std::vector<Polynomial>& different,
    std::size_t n,const std::set<std::size_t>& poles,const std::vector<std::size_t>& positiveAxes,
    const std::vector<Polynomial>& graph) {
    std::ostringstream out;
    out<<"(set-logic QF_NRA)\n";
    for (char prefix : {'a','b'}) {
        for (std::size_t i=0;i<n;++i) out<<"(declare-const "<<prefix<<i<<" Real)\n";
        for (auto i:poles) out<<"(assert (distinct "<<prefix<<i<<" 0))\n";
        for (auto i:positiveAxes) out<<"(assert (> "<<prefix<<i<<" 0))\n";
        for (const auto& p:graph) out<<"(assert (= "<<numeratorText(p,Powers(n),prefix)<<" 0))\n";
    }
    for (const auto& p:jet) out<<"(assert "<<equalValues(p,n)<<")\n";
    out<<"(assert (not (and";
    for (const auto& p:different) out<<" "<<equalValues(p,n);
    out<<")))\n(check-sat)\n";
    return out.str();
}
json::Value decideFiber(const std::string& problem,std::size_t timeoutMs) {
    std::filesystem::path configured;
#ifdef _WIN32
    wchar_t* path=nullptr;
    std::size_t size=0;
    if (_wdupenv_s(&path,&size,L"DYNSYSLAB_Z3_LIBRARY")==0 && path) {
        std::unique_ptr<wchar_t,decltype(&std::free)> owner(path,&std::free);
        configured=path;
    }
#else
    if (const char* path=std::getenv("DYNSYSLAB_Z3_LIBRARY")) configured=path;
#endif
    if (configured.empty()) return json::Value::Object{{"status","backend_unavailable"}};
    try {
        Library lib(configured);
        auto makeConfig=lib.symbol<_Z3_config*(*)()>("Z3_mk_config");
        auto delConfig=lib.symbol<void(*)(_Z3_config*)>("Z3_del_config");
        auto makeContext=lib.symbol<_Z3_context*(*)(_Z3_config*)>("Z3_mk_context_rc");
        auto delContext=lib.symbol<void(*)(_Z3_context*)>("Z3_del_context");
        auto run=lib.symbol<const char*(*)(_Z3_context*,const char*)>("Z3_eval_smtlib2_string");
        std::unique_ptr<_Z3_config,decltype(delConfig)> config(makeConfig(),delConfig);
        std::unique_ptr<_Z3_context,decltype(delContext)> context(makeContext(config.get()),delContext);
        const auto command="(set-option :timeout "+std::to_string(timeoutMs)+")\n"+problem;
        const char* reply=run(context.get(),command.c_str());
        std::string answer=reply?reply:"";
        answer.erase(answer.find_last_not_of(" \n\r\t")+1);
        if (answer=="sat") {
            const char* model=run(context.get(),"(get-model)");
            return json::Value::Object{{"status","sat"},{"model",model?model:""}};
        }
        if (answer=="unsat" || answer=="unknown") return json::Value::Object{{"status",answer}};
        return json::Value::Object{{"status","solver_error"},{"reason",answer}};
    } catch (const std::exception& e) {
        return json::Value::Object{{"status","backend_unavailable"},{"reason",e.what()}};
    }
}
} // namespace dynsyslab::core::jet
