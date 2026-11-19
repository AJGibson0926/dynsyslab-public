#pragma once
#include "PostprocessInternal.hpp"
#include <dynsyslab/core/KoopmanModel.hpp>

namespace dynsyslab::workflow::postprocess_detail {
std::string datasetDynamics(const SourceData& source);

// Training-data role of a post-processing method (docs/architecture/TRAINING_DATA_ROLES.md).
//   autonomous      fits a Koopman spectrum / eigenfunctions / free dynamics: unforced data only
//   forced_allowed  DMDc / EDMDc: identify eigenfunctions jointly with the input from forced data
//   forced_required regresses an input map B_psi: forced data only
//   any             data-agnostic (decompositions, coordinate maps, reconstruction, analysis)
std::string trainingDataRole(const std::string& method, const json::Value& settings);
// Throws unless the source's dynamics suit the method's role: autonomous methods refuse
// forced sources; input-map regressions require forced ones. Sources without forcing
// provenance ("unknown": legacy, synthetic or measured records) are accepted.
void requireTrainingDataRole(const std::string& method, const SourceData& source, const json::Value& settings);
// Same check from an archive path, reading only its metadata and recorded control history
// (so measurement-only archives without a state array can be classified).
void requireTrainingDataRole(const std::string& method, const std::filesystem::path& source, const json::Value& settings);
void requireMatchingUnforcedDataset(const dynsyslab::core::KoopmanModel& model,
                                    const SourceData& source);
// Select existing samples only; never interpolates or edits a source archive.
SourceData validationSamples(const SourceData& source,
    const dynsyslab::core::KoopmanModel& model,const json::Value& settings);
}
