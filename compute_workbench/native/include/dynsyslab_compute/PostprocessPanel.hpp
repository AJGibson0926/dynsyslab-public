#pragma once

#include "dynsyslab_compute/Registry.hpp"
#include "dynsyslab_viewer/Theme.hpp"

#include <windows.h>
#include <filesystem>
#include <string>
#include <vector>

namespace dsl::compute {

constexpr int kPostPlanId = 5301;
constexpr int kPostRunId = 5302;
constexpr int kKromExtractId = 5303;
constexpr int kKoopmanPipelineId = 5304;
constexpr int kKromModesId = 5305;
constexpr int kControlRegressionSettingsId = 5306;
constexpr int kEstimationSettingsId = 5307;
constexpr int kDecompositionSettingsId = 5308;
constexpr int kKromValidationSettingsId = 5309;
constexpr int kBoundaryEigenfunctionSettingsId = 5310;
constexpr int kEigenvalueOptimizationSettingsId = 5311;
constexpr int kKromTrainingDomainId = 5312;
constexpr int kKromControlSearchId = 5313;
constexpr int kDmdControlSettingsId = 5314;
constexpr int kDmdControlStageId = 5315;
constexpr int kEdmdControlSettingsId = 5316;
constexpr int kEdmdControlStageId = 5317;
constexpr int kSpectralKoopmanSettingsId = 5318;

class PostprocessPanel {
public:
    bool create(HWND parent, const dsl::theme::Fonts& fonts, const Registry* registry);
    void show(bool visible);
    void layout(int x, int y, int width, int height);
    void setSource(const std::filesystem::path& source);
    std::filesystem::path source() const;
    void command(int id, int notification);
    std::string jobJson() const;
private:
    HWND parent_ = nullptr;
    const dsl::theme::Fonts* fonts_ = nullptr;
    const Registry* registry_ = nullptr;
    std::vector<HWND> controls_;
    std::vector<HWND> labels_;
    HWND source_ = nullptr; HWND label_ = nullptr;
    HWND hdmd_ = nullptr; HWND sindy_ = nullptr; HWND havok_ = nullptr;
    HWND dmd_ = nullptr; HWND edmd_ = nullptr; HWND pod_ = nullptr; HWND pca_ = nullptr; HWND svd_ = nullptr;
    HWND dmdControl_ = nullptr; HWND dmdControlButton_ = nullptr;
    HWND edmdControl_ = nullptr; HWND edmdControlButton_ = nullptr;
    HWND spectralKoopman_ = nullptr; HWND spectralKoopmanButton_ = nullptr;
    std::string spectralKoopmanSettingsJson_{"{\"method\":\"resdmd\",\"dictionary\":{\"kind\":\"linear\"},\"rank\":0}"};
    HWND baggedPca_ = nullptr; HWND robustPca_ = nullptr; HWND hankelPod_ = nullptr;
    HWND controlRegression_ = nullptr; HWND kromExtract_ = nullptr;
    HWND kromValidation_ = nullptr; HWND boundaryEigenfunctions_ = nullptr;
    HWND eigenvalueOptimization_ = nullptr;
    HWND observability_ = nullptr; HWND klqe_ = nullptr;
    HWND hdmdSettings_ = nullptr; HWND sindySettings_ = nullptr; HWND dmdSettings_ = nullptr; HWND havokSettings_ = nullptr;
    HWND edmdSettings_ = nullptr; HWND pcaSettings_ = nullptr;
    HWND controlRegressionSettings_ = nullptr; HWND kromExtractSettings_ = nullptr;
    HWND estimationSettings_ = nullptr;
    HWND archive_ = nullptr; HWND pipelineButton_ = nullptr; HWND kromModesButton_ = nullptr;
    HWND controlRegressionButton_ = nullptr;
    HWND kromValidationButton_ = nullptr; HWND boundaryEigenfunctionButton_ = nullptr;
    HWND eigenvalueOptimizationButton_ = nullptr;
    HWND trainingDomainButton_ = nullptr;
    HWND kromSearchButton_ = nullptr;
    HWND estimationButton_ = nullptr;
    HWND decompositionButton_ = nullptr;
    std::string svdSettingsJson_{"{\"rank\":0,\"center\":false}"};
    std::string dmdControlSettingsJson_{"{\"snapshot_array\":\"values\",\"rank\":0,\"augmented_rank\":0,\"svd_tolerance\":0,\"save_eigenfunction_values\":true}"};
    std::string edmdControlSettingsJson_{"{\"snapshot_array\":\"values\",\"dictionary\":\"polynomial\",\"polynomial_degree\":4,\"include_constant\":true,\"normalize\":true,\"rank\":0,\"augmented_rank\":0,\"save_eigenfunction_values\":true}"};
    std::string podSettingsJson_{"{\"rank\":0}"};
    std::string baggedPcaSettingsJson_{"{\"rank\":0,\"bagged_center\":true,\"bags\":32,\"seed\":0}"};
    std::string robustPcaSettingsJson_{"{\"lambda\":0,\"tolerance\":1e-7,\"max_iterations\":1000}"};
    std::string hankelPodSettingsJson_{"{\"rank\":0,\"matrix_size\":10,\"delay_stride\":1,\"center\":false}"};
    std::string kromValidationSettingsJson_{"{\"krom\":\"\",\"dataset_role\":\"test\",\"horizons\":[]}"};
    std::string boundaryEigenfunctionSettingsJson_{"{\"spectrum\":\"\",\"observable\":\"state_0\",\"additional_sources\":[],\"mode_groups\":[],\"interpolation_neighbors\":8,\"gradient_relative_step\":1e-6,\"svd_tolerance\":1e-12}"};
    std::string eigenvalueOptimizationSettingsJson_{"{\"initial_spectrum\":\"\",\"observable\":\"state_0\",\"additional_sources\":[],\"mode_groups\":[],\"lower_bound\":-1.05,\"upper_bound\":1.05,\"max_iterations\":500,\"memory\":8,\"tolerance\":1e-10,\"preserve_conjugates\":true}"};
    void editKoopmanPipeline();
    void editKromModes();
    void editControlRegression();
    void editEstimation();
    void editDecompositions();
    void editDmdControl(bool extended = false);
    void editKromValidation();
    void editBoundaryEigenfunctions();
    void editEigenvalueOptimization();
    void add(HWND value) { controls_.push_back(value); }
};

}  // namespace dsl::compute
