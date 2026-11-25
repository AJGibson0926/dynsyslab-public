#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/ControlRegression.hpp>
#include <dynsyslab/core/Json.hpp>
#include <dynsyslab/core/KoopmanArchive.hpp>
#include <dynsyslab/core/KoopmanValidation.hpp>
#include <dynsyslab/core/KromScreening.hpp>
#include <dynsyslab/core/LQR.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <memory>
#include <stdexcept>

namespace {

dynsyslab::core::KoopmanModel velocityModel() {
    dynsyslab::core::KoopmanSINDyResult fit;
    fit.continuousEigenvalues = {{-1.0, 0.0}};
    fit.coefficients = {{1.0, 0.0}};
    fit.exponents = {0, 1};
    fit.center = {0.0, 0.0}; fit.scale = {1.0, 1.0};
    fit.stateDimension = 2; fit.featureCount = 1;
    dynsyslab::core::SystemConfiguration system;
    system.identifier = "duffing"; system.preset = "double_well_reference";
    return {std::move(fit), std::move(system)};
}

dynsyslab::core::KoopmanModel energyCoordinateModel() {
    dynsyslab::core::KoopmanSINDyResult fit;
    fit.continuousEigenvalues={{-1.0,0.0}}; fit.coefficients={{1.0,0.0}};
    fit.exponents={1}; fit.center={0.0}; fit.scale={1.0};
    fit.stateDimension=1; fit.dynamicalStateDimension=2; fit.featureCount=1;
    fit.coordinateObservableNames={"total_energy"};
    dynsyslab::core::SystemConfiguration system;
    system.identifier="duffing"; system.preset="double_well_reference";
    return {std::move(fit),std::move(system)};
}

dynsyslab::core::KoopmanModel screeningModel() {
    dynsyslab::core::KoopmanSINDyResult fit;
    fit.continuousEigenvalues={{-1.0,0.0},{-2.0,0.0}};
    fit.coefficients={{0.0,0.0},{0.0,0.0},{1.0,0.0},
                      {0.0,0.0},{1.0,0.0},{0.0,0.0}};
    fit.exponents={0,0,1,0,0,1};fit.center={0.0,0.0};fit.scale={1.0,1.0};
    fit.stateDimension=2;fit.featureCount=3;
    dynsyslab::core::SystemConfiguration system;
    system.identifier="duffing";system.preset="double_well_reference";
    return {std::move(fit),std::move(system)};
}

void close(double actual, double expected, double tolerance, const char* message) {
    if (std::abs(actual - expected) > tolerance) throw std::runtime_error(message);
}

void legacyArchiveCompatibility() {
    namespace core=dynsyslab::core;
    const auto root=std::filesystem::temp_directory_path()
        / "dynsyslab-native-legacy-krom-test";
    std::filesystem::remove_all(root); std::filesystem::create_directories(root);
    const auto kromPath=root/"legacy-krom.npz";
    core::NpzRecord krom;
    krom.arrays["eigenvalues"]=core::NpyArray::fromComplex(
        {2},{std::exp(std::complex<double>{-0.5,1.0}*0.01),
             std::exp(std::complex<double>{-0.5,-1.0}*0.01)});
    krom.arrays["coefficients"]=core::NpyArray::fromComplex(
        {2,3},{{0.0,0.0},{1.0,0.0},{0.0,1.0},
               {0.0,0.0},{1.0,0.0},{0.0,-1.0}});
    krom.arrays["training_states"]=core::NpyArray::fromDoubles({1,2},{0.0,0.0});
    krom.arrays["metadata_json"]=core::NpyArray::fromString(
        R"({"analysis_metadata":{"dictionary":{"degree":1,"include_constant":true},"state_preprocessing":{"center":[0,0],"scale":[1,1],"multiplier":[1,1]},"krom":{"delta_t":0.01,"eigenvalue_domain":"continuous","retained_indices":[0]}}})");
    core::writeNpzAtomic(kromPath,krom);
    core::SystemConfiguration configuration;
    configuration.identifier="duffing";
    configuration.preset="gibson_damped_stabilization_2022";
    auto imported=core::loadPolynomialKoopmanModel(kromPath,configuration);
    if (imported.liftedDimension()!=1 || imported.splitDimension()!=2)
        throw std::runtime_error("legacy explicit conjugate dual was not canonicalized");
    const double point[]{0.25,-0.5};
    const auto psi=imported.evaluate(point);
    close(psi[0].real(),0.25,1e-12,"legacy KROM real evaluation");
    close(psi[0].imag(),-0.5,1e-12,"legacy KROM imaginary evaluation");

    const auto mapsPath=root/"legacy-control-maps.npz";
    core::NpzRecord maps;
    maps.arrays["evaluation_0_matrix"]=core::NpyArray::fromComplex(
        {1,1},{{0.2,0.3}});
    maps.arrays["evaluation_1_operators"]=core::NpyArray::fromComplex(
        {1,1,1},{{0.1,0.0}});
    maps.arrays["evaluation_1_offsets"]=core::NpyArray::fromComplex(
        {1,1},{{0.2,0.3}});
    maps.arrays["evaluation_1_conjugate_operators"]=core::NpyArray::fromComplex(
        {1,1,1},{{0.0,0.1}});
    core::writeNpzAtomic(mapsPath,maps);
    imported.setControlMap("constant",core::loadApproximateControlMap(
        mapsPath,"evaluation_0",1,1));
    imported.setControlMap("bilinear",core::loadApproximateControlMap(
        mapsPath,"evaluation_1",1,1));
    const auto constant=imported.controlMapSplit(point,"constant");
    close(constant[0],0.2,1e-12,"legacy constant control map real");
    close(constant[1],0.3,1e-12,"legacy constant control map imaginary");
    const auto bilinear=imported.controlMapSplit(point,"bilinear");
    if (bilinear.size()!=2 || !std::isfinite(bilinear[0])
        || !std::isfinite(bilinear[1]))
        throw std::runtime_error("legacy bilinear control map evaluation");
    std::filesystem::remove_all(root);
}

void selfContainedArchiveRoundTrip() {
    namespace core=dynsyslab::core;
    const auto root=std::filesystem::temp_directory_path()
        / "dynsyslab-native-self-contained-krom-test";
    std::filesystem::remove_all(root); std::filesystem::create_directories(root);
    const auto path=root/"independent-krom.npz";
    auto model=velocityModel();
    core::ConstantControlMap constant{{{2.0,0.5}},1,1,"round-trip test"};
    constant.discreteMatrix={{0.2,0.05}};constant.sampleTime=0.1;
    model.setControlMap("constant",std::move(constant));
    model.setReconstruction("state",{{{1.0,0.0},{0.0,1.0}},2,1,0.0,0.0,0.0});
    core::savePolynomialKoopmanModel(path,model,R"({"label":"Independent KROM"})");
    const auto metadata=core::json::parse(core::readNpz(path).metadataJson());
    if (metadata.at("system_configuration").at("identifier").string()!="duffing"
        || metadata.at("krom_mode_count").number()!=1.0)
        throw std::runtime_error("self-contained KROM system metadata");
    const auto loaded=core::loadPolynomialKoopmanModel(path);
    if (loaded.systemConfiguration().identifier!="duffing"
        || loaded.systemConfiguration().parameters.empty()
        || !loaded.hasControlMap("constant") || !loaded.hasReconstruction("state"))
        throw std::runtime_error("self-contained KROM round trip");
    const auto& restored=std::get<core::ConstantControlMap>(
        loaded.storedControlMap("constant"));
    if(restored.discreteMatrix.size()!=1||restored.sampleTime!=0.1)
        throw std::runtime_error("constant control map time-domain evidence round trip");
    std::filesystem::remove_all(root);
}

void scatteredArchiveStoresOneValueArray() {
    namespace core=dynsyslab::core;
    const auto root=std::filesystem::temp_directory_path()
        / "dynsyslab-native-scattered-krom-deduplication-test";
    std::filesystem::remove_all(root);std::filesystem::create_directories(root);
    const auto path=root/"scattered-krom.npz";
    core::InterpolationData fit;
    fit.trainingStates={-1.0,0.0,0.0,0.5,1.0,0.0};
    fit.trainingValues={{-1.0,0.0},{0.0,0.5},{1.0,0.0}};
    fit.samples=3;fit.stateDimension=2;fit.neighbors=2;
    core::SystemConfiguration system;
    system.identifier="duffing";system.preset="double_well_reference";
    core::saveKoopmanModel(path,core::KoopmanModel(std::move(fit),core::Spectrum::fromContinuousRates({{-0.5,0.25}}),std::move(system)));
    const auto archive=core::readNpz(path);
    if(!archive.contains("eigenfunction_values")
        ||archive.contains("training_eigenfunction_values"))
        throw std::runtime_error("scattered KROM duplicated its sampled value array");
    const auto loaded=core::loadKoopmanModel(path);
    const double state[]{0.0,0.5};
    const auto value=loaded.evaluate(state);
    if(value.size()!=1||std::abs(value.front()-std::complex<double>{0.0,0.5})>1e-12)
        throw std::runtime_error("deduplicated scattered KROM did not round trip");
    std::filesystem::remove_all(root);
}

void semigroupValidation() {
    dynsyslab::core::KoopmanSINDyResult fit;
    fit.continuousEigenvalues={{0.0,0.0}}; fit.coefficients={{1.0,0.0}};
    fit.exponents={0,0}; fit.center={0.0,0.0}; fit.scale={1.0,1.0};
    fit.stateDimension=2; fit.featureCount=1;
    dynsyslab::core::SystemConfiguration configuration;
    configuration.identifier="duffing"; configuration.preset="double_well_reference";
    dynsyslab::core::KoopmanModel model(std::move(fit),std::move(configuration));
    const std::vector<double> time{0.0,0.1,0.25,0.5};
    const std::vector<double> states{-1.0,0.0,-0.8,0.2,0.1,-0.4,1.2,0.0};
    const auto result=dynsyslab::core::validateKoopmanSemigroup(
        model,time,states,4,2,{0.1,0.25},{-1.0,-0.5},{1.0,0.5});
    if(result.modes.size()!=1||result.modes[0].horizons.size()!=2
        ||result.outsideTrainingSamples!=1)
        throw std::runtime_error("Koopman semigroup validation shape");
    close(result.modes[0].oneStep.relativeL2,0.0,1e-14,
          "constant eigenfunction one-step error");
    close(result.modes[0].anchoredLongTerm.relativeL2,0.0,1e-14,
          "constant eigenfunction anchored error");
    close(result.modes[0].generator.relativeL2,0.0,1e-14,
          "constant eigenfunction generator error");
    if(result.modes[0].oneStepChangeDefined)
        throw std::runtime_error("constant observable has no nonzero change normalization");
    const auto wrong=dynsyslab::core::KoopmanModel::fromPolynomialCoefficients(
        model.systemConfiguration(),2,{1,0},{{1.,0.}},{{0.,0.}});
    const auto changes=dynsyslab::core::validateKoopmanSemigroup(wrong,
        {0.,.1,.2},{1.,0.,2.,0.,4.,0.},3,2,{.1});
    close(changes.modes[0].oneStepChangeRelativeL2,1.,1e-14,
        "identity prediction error equals actual change for wrong constant-rate mode");
    const auto exact=dynsyslab::core::KoopmanModel::fromPolynomialCoefficients(
        model.systemConfiguration(),2,{1,0},{{1.,0.}},{{std::log(2.)/.1,0.}});
    const auto perfect=dynsyslab::core::validateKoopmanSemigroup(exact,
        {0.,.1,.2},{1.,0.,2.,0.,4.,0.},3,2,{.1});
    close(perfect.modes[0].oneStepChangeRelativeL2,0.,1e-14,"exact exponential change-relative error");
}

// KLQR adapters on a scalar KROM psi = x0 (lambda = -1), constant B = 1, readout "state".
void klqrIntegralAndPidAdapters() {
    namespace core = dynsyslab::core;
    core::SystemConfiguration system;
    system.identifier = "duffing"; system.preset = "double_well_reference";
    auto build = [&] {
        auto model = std::make_shared<core::KoopmanModel>(core::KoopmanModel::fromPolynomialCoefficients(
            system, 2, {1, 0}, {{1.0, 0.0}}, {{-1.0, 0.0}}));
        core::ReconstructionResult state{{{1.0, 0.0}, {0.0, 0.0}}, 2, 1};
        state.offset = {0.0, 0.0};
        model->setReconstruction("state", state);
        model->setControlMap("constant", core::ConstantControlMap{{{1.0, 0.0}}, 1, 1, "adapter test"});
        return model;
    };
    const auto target = [](double) { return std::vector<double>{1.0, 0.0}; };
    core::LiftedOutputChannel channel{{1.0}, {0.0}, {0}};
    core::KoopmanLQRController plain(build(), {10.0}, {1.0}, target, "constant");
    core::KoopmanLQRController integral(build(), {10.0}, {1.0}, target, "constant");
    integral.setIntegralAction({channel, {10.0}, {}});
    core::KoopmanPIDCorrection pid{channel, {0.0}, {2.0}, {0.0}};
    core::KoopmanLQRController corrected(build(), {10.0}, {1.0}, target, "constant");
    corrected.setPIDCorrection(pid);
    if (!integral.stateful() || !corrected.stateful() || plain.stateful())
        throw std::runtime_error("KLQR augmentation flags");
    // Plant psi' = -psi + u + d (an input disturbance the KROM does not know).
    const double d = 0.5, dt = 1e-3;
    double a = 0.0, b = 0.0, c = 0.0;
    for (int k = 0; k < 40000; ++k) {
        const double t = k * dt;
        const double ua = plain.valueFromLifted(t, {{a, 0.0}}).front();
        const double ub = integral.valueFromLifted(t, {{b, 0.0}}).front();
        const double uc = corrected.valueFromLifted(t, {{c, 0.0}}).front();
        a += dt * (-a + ua + d); b += dt * (-b + ub + d); c += dt * (-c + uc + d);
    }
    if (std::abs(a - 1.0) < 0.05) throw std::runtime_error("plain KLQR should keep an offset");
    close(b, 1.0, 1e-5, "KLQR integral action removes the offset");
    close(c, 1.0, 1e-5, "KLQR PID correction removes the offset");
    if (integral.lastSolution().gain.size() != 2) throw std::runtime_error("LQI gain on the augmented KROM");
}

}  // namespace

int main() {
    try {
        klqrIntegralAndPidAdapters();
        auto model = velocityModel(); const double point[]{0.25, -0.5};
        legacyArchiveCompatibility();
        selfContainedArchiveRoundTrip();
        scatteredArchiveStoresOneValueArray();
        semigroupValidation();
        close(model.evaluate(point)[0].real(), -0.5, 1e-12, "eigenfunction evaluation");
        close(model.exactControlMap(point)[0].real(), 1.0, 1e-12, "exact B_psi");
        auto energy=energyCoordinateModel();
        close(energy.evaluate(point)[0].real(),
              0.5*point[1]*point[1]-0.5*point[0]*point[0]+0.25*std::pow(point[0],4),
              1e-12,"observable-coordinate eigenfunction");
        close(energy.exactControlMap(point)[0].real(),point[1],1e-12,
              "observable-coordinate B_psi chain rule");
        energy.setReconstruction("energy",{{{1.0,0.0}},1,1,0.0,0.0,0.0});
        close(energy.reconstruct("energy",point)[0].real(),energy.evaluate(point)[0].real(),
              1e-12,"KROM observable reconstruction");
        const std::vector<double> samples{-1,-1, 0,0, 1,1, 0.5,-0.5};
        const auto constant = dynsyslab::core::fitConstantExactControlMap(model, samples, 4);
        close(constant.representation.matrix[0].real(), 1.0, 1e-12, "constant exact fit");
        dynsyslab::core::ControlledData controlled;
        controlled.samples=8;controlled.stateDimension=2;controlled.inputDimension=1;
        controlled.inputs.assign(controlled.samples-1,2.0);
        const double step=0.1,transition=std::exp(-step),integral=1.0-transition;
        double velocity=0.3;
        for(std::size_t index=0;index<controlled.samples;++index){
            controlled.time.push_back(index*step);controlled.states.push_back(0.0);
            controlled.states.push_back(velocity);velocity=transition*velocity+integral*2.0;
        }
        auto system=dynsyslab::core::createSystem(model.systemConfiguration());
        const auto multistep=dynsyslab::core::fitConstantControlMultistepOutput(
            model,controlled,system->observableCatalog().at("state_1"),
            {{{1.0,0.0}},1,1,0.0,0.0,0.0});
        close(multistep.representation.matrix[0].real(),1.0,1e-10,
              "multistep continuous control fit");
        close(multistep.representation.discreteMatrix[0].real(),integral,1e-10,
              "multistep discrete control fit");
        close(multistep.representation.sampleTime,step,1e-14,
              "multistep sample time");
        model.setControlMap("constant", constant.representation);
        const auto bilinear = dynsyslab::core::fitBilinearExactControlMap(model, samples, 4);
        model.setControlMap("bilinear", bilinear.representation);
        const auto errors = dynsyslab::core::evaluateControlMapError2D(
            model, "bilinear", -1, 1, -1, 1, 7, 7);
        if (errors.absoluteMaximum > 1e-10) throw std::runtime_error("bilinear error grid");
        const auto candidates=dynsyslab::core::rankKromControlCandidates(
            screeningModel(),{0.0,0.0},{0.0,1.0},{2,20,3,2});
        if(candidates.size()!=3||candidates.front().modeIndices!=std::vector<std::size_t>{0}
            ||candidates.front().minimumExactInputAuthority<=0.0)
            throw std::runtime_error("native KROM candidate screening");
        const auto reduced = model.extract({0});
        if (!reduced.hasControlMap("constant") || !reduced.hasControlMap("bilinear"))
            throw std::runtime_error("sub-KROM control map inheritance");
        const auto care = dynsyslab::core::solveContinuousLQR(
            {0.0}, {1.0}, {1.0}, {1.0}, 1, 1);
        close(care.riccati[0], 1.0, 1e-10, "scalar CARE");
        close(care.gain[0], 1.0, 1e-10, "scalar LQR gain");
        auto shared = std::make_shared<dynsyslab::core::KoopmanModel>(std::move(model));
        dynsyslab::core::KoopmanLQRController controller(shared, {1.0}, {1.0},
            [](double) { return std::vector<double>{0.0, 0.0}; }, "exact");
        const auto input = controller.value(0.0, point);
        if (input.size() != 1 || !std::isfinite(input[0])) throw std::runtime_error("KLQR value");
        dynsyslab::core::KoopmanMPCController mpc(shared,{1.0},{1.0},
            [](double) { return std::vector<double>{0.0,0.0}; },12,0.01,"exact");
        const auto mpcInput=mpc.value(0.0,point);
        if (mpcInput.size()!=1||!std::isfinite(mpcInput[0])||mpc.lastGain().size()!=1)
            throw std::runtime_error("Koopman MPC value");
        std::cout << "dynsyslab_koopman_control_tests: PASS\n"; return 0;
    } catch (const std::exception& error) {
        std::cerr << "dynsyslab_koopman_control_tests: FAIL: " << error.what() << '\n'; return 1;
    }
}
