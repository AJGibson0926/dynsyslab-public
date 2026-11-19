// Viewer views of a linear_systems archive.
#include "dynsyslab_viewer/NativePayload.hpp"

#include <dynsyslab/core/Archive.hpp>

#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) { if (!condition) throw std::runtime_error(message); }
}  // namespace

int main() {
    const auto path = std::filesystem::temp_directory_path() / "dynsyslab-linear-systems-view-test.npz";
    try {
        using dynsyslab::core::NpyArray;
        dynsyslab::core::NpzRecord a;
        a.arrays["system_poles"] = NpyArray::fromDoubles({3, 2}, {-1, 0, -5, 0, 0.2, 0});
        a.arrays["reduced_poles"] = NpyArray::fromDoubles({2, 2}, {-1.01, 0, 0.2, 0});
        a.arrays["transmission_zeros"] = NpyArray::fromDoubles({1, 2}, {-3, 0});
        a.arrays["hankel_singular_values"] = NpyArray::fromDoubles({2}, {0.5, 0.01});
        a.arrays["hankel_error_bound"] = NpyArray::fromDoubles({3}, {1.02, 0.02, 0});
        a.arrays["controllability_gramian_eigenvalues"] = NpyArray::fromDoubles({2}, {0.6, 0.05});
        a.arrays["observability_gramian_eigenvalues"] = NpyArray::fromDoubles({2}, {0.4, 0.01});
        a.arrays["frequencies"] = NpyArray::fromDoubles({3}, {0.1, 1, 10});
        a.arrays["full_singular_values"] = NpyArray::fromDoubles({3, 1}, {1, 0.7, 0.1});
        a.arrays["reduced_singular_values"] = NpyArray::fromDoubles({3, 1}, {1, 0.69, 0.1});
        a.arrays["full_phase"] = NpyArray::fromDoubles({3}, {0, -0.7, -1.5});
        a.arrays["lqr_closed_loop_poles"] = NpyArray::fromDoubles({1, 2}, {-2, 0});
        a.arrays["lqg_closed_loop_poles"] = NpyArray::fromDoubles({2, 2}, {-2, 0, -4, 0});
        a.arrays["reduced_lqg_full_closed_loop_poles"] = NpyArray::fromDoubles({2, 2}, {-2.1, 0, -3.9, 0});
        a.arrays["bode_frequencies"] = NpyArray::fromDoubles({3}, {0.1, 1, 10});
        a.arrays["bode_magnitude_db"] = NpyArray::fromDoubles({3, 1}, {0, -3, -20});
        a.arrays["bode_phase_deg"] = NpyArray::fromDoubles({3, 1}, {-5, -45, -200});
        for (const char* loop : {"lqr_input", "lqg_input", "lqg_output", "ltr_0_input"}) {
            const std::string prefix = std::string("loop__") + loop + "__";
            a.arrays[prefix + "frequencies"] = NpyArray::fromDoubles({4}, {0.1, 1, 10, 100});
            a.arrays[prefix + "response"] = NpyArray::fromDoubles({4, 2}, {20, -5, 1, -2, -0.5, -0.3, 0.01, -0.001});
        }
        a.arrays["metadata_json"] = NpyArray::fromString(
            R"({"label":"linear systems view","method":"linear-systems","schema":"dynsyslab.linear-systems/v1",)"
            R"("time_domain":"continuous","gramian_coordinates":"model coordinates","reduction_error_hinf":0.011,)"
            R"("bode_channels":["echo / u0"],"transfer_functions":[{"channel":"y0_u0","output":"echo","input":"u0","gain":2,"relative_degree":1}],)"
            R"("loop_names":["lqr_input","lqg_input","lqg_output","ltr_0_input"],)"
            R"("loop_robustness":{"lqr_input":{"minimum_return_difference":1},"lqg_input":{"minimum_return_difference":0.6},)"
            R"("lqg_output":{"minimum_return_difference":0.6},"ltr_0_input":{"minimum_return_difference":0.7,"ltr_weight":100}},)"
            R"("reduction":{"method":"truncation","order":2,"stable_order":1,"unstable_order":1,"error_bound":0.02}})");
        dynsyslab::core::writeNpzAtomic(path, a);
        const auto hankel = dsl::buildNativePayload(path, "hankel_spectrum");
        require(hankel.figures.size() == 1 && hankel.figures[0].traces.size() == 3, "Hankel figure traces");
        require(hankel.figures[0].valueScale == "log", "Hankel log scale");
        require(hankel.figures[0].traces[2].x.size() == 2 && hankel.figures[0].traces[2].x[0] == 1.0, "chosen order line");
        const auto gramian = dsl::buildNativePayload(path, "gramian_spectrum");
        require(gramian.figures.size() == 2 && gramian.figures[1].traces.size() == 2, "gramian spectra and ellipses");
        const auto frequency = dsl::buildNativePayload(path, "frequency_response");
        require(frequency.figures.size() == 2 && frequency.figures[0].traces.size() == 2, "frequency response figures");
        require(frequency.figures[0].traces[0].x[1] == 0.0, "log10 frequency axis");
        const auto poles = dsl::buildNativePayload(path, "pole_zero");
        require(poles.figures.size() == 3 && poles.figures[0].traces.size() == 3, "pole-zero figures");
        require(poles.figures[2].traces.size() == 2, "output-feedback pole figure");
        require(poles.figures[0].traces[0].lineStyle == "none", "poles drawn as a scatter");
        const auto bode = dsl::buildNativePayload(path, "bode");
        require(bode.figures.size() == 2 && bode.figures[0].traces.size() == 1 && bode.figures[1].traces[0].y[2] == -200,
                "plant Bode magnitude and phase");
        const auto nyquist = dsl::buildNativePayload(path, "nyquist");
        require(nyquist.figures.size() == 2, "input and output Nyquist figures");
        const auto& in = nyquist.figures[0];  // disk, LQR, mirror, LQG, mirror, LTR, mirror, -1
        require(in.traces.size() == 8 && in.square, "input Nyquist traces");
        require(std::isnan(in.traces[1].x[0]) && in.traces[1].x[1] == 1.0, "|L| > 10 not drawn");
        require(in.traces[1].color == "#004DB8" && in.traces[3].color == "#D8242A", "LQR blue, LQG red");
        require(in.traces.back().x[0] == -1.0, "critical point");
        const auto sens = dsl::buildNativePayload(path, "sensitivity");
        require(sens.figures.size() == 2 && sens.figures[0].traces.size() == 4, "sensitivity per loop");
        require(std::abs(sens.figures[0].traces[0].y[3] - 20 * std::log10(1 / std::abs(std::complex<double>(1.01, -0.001)))) < 1e-9,
                "|S| in dB");
        const auto lb = dsl::buildNativePayload(path, "loop_bode");
        require(lb.figures.size() == 2 && lb.figures[1].traces[0].y.size() == 4, "loop Bode");
        std::filesystem::remove(path);
        std::cout << "dynsyslab_linear_systems_payload_tests: PASS\n";
        return 0;
    } catch (const std::exception& error) {
        std::error_code ignored; std::filesystem::remove(path, ignored);
        std::cerr << "dynsyslab_linear_systems_payload_tests: FAIL: " << error.what() << '\n';
        return 1;
    }
}
