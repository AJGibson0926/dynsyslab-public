#pragma once
// Receiver-array post-processing methods (ReceiverPostprocess.cpp): receiver_signals,
// state_readout and ring_errors, each run per source.
#include <dynsyslab/core/Json.hpp>

#include <filesystem>
#include <string>

namespace dynsyslab::workflow::postprocess_detail {

std::filesystem::path runReceiverSignals(const std::filesystem::path& source, const dynsyslab::core::json::Value& settings,
                                         const std::filesystem::path& root, const std::string& label);
std::filesystem::path runStateReadout(const std::filesystem::path& source, const dynsyslab::core::json::Value& settings,
                                      const std::filesystem::path& root, const std::string& label);
std::filesystem::path runRingErrors(const std::filesystem::path& source, const dynsyslab::core::json::Value& settings,
                                    const std::filesystem::path& root, const std::string& label);

}  // namespace dynsyslab::workflow::postprocess_detail
