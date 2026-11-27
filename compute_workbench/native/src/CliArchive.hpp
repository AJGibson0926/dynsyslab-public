#pragma once
// dynsyslab-cli archive ... and storage ... (archive management; docs/ARCHIVE_MANAGEMENT.md).
//   archive index | audit [--report DIR] | deduplicate [--apply]
//   archive recompress [FILE...] [--compression C] [--apply]
//   archive share-time-grids [--apply]
//   archive edit-metadata <source.npz> <output.npz> --metadata-json JSON
//   archive campaigns [--output FILE]
//   storage list | register <name> <folder> [--read-only] [--purpose P]
//   storage copy <file> <folder> [--remove-source]
//   storage relocate <campaign-id> <location> | storage append <plan.json> <file>...
// Scope options (audit, deduplicate, recompress, share-time-grids): --days N | --all-ages,
// --section NAME, --system NAME, --campaign ID, --include-derived. Each prints one JSON line.
#include "CliArguments.hpp"

namespace dynsyslab::cli {

void archiveCommand(const Arguments& args);
void storageCommand(const Arguments& args);

}  // namespace dynsyslab::cli
