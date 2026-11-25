#pragma once

namespace dsl {

// Loads the workflow library once so readNpz can replay recorded archives (RecordedArchives.cpp).
void loadRecordedArchiveReaders() noexcept;

}  // namespace dsl
