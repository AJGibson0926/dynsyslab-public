#pragma once
#include <functional>
#include <vector>
namespace dynsyslab::core::detail {
double analyticIntegral(const std::function<double(double)> &f, double a, double b);
// Periodic integrand, unwrapped primitive. A positive integrand additionally gives
// a monotone clock which can be inverted without stepping an ODE.
class PeriodicPrimitive {
  public:
    PeriodicPrimitive(std::function<double(double)> rate, double period, bool positive = false);
    double value(double parameter) const;
    double inverse(double value) const;
    double increment() const { return table_.back(); }

  private:
    std::function<double(double)> rate_;
    double period_;
    bool positive_;
    std::vector<double> table_;
};
} // namespace dynsyslab::core::detail
