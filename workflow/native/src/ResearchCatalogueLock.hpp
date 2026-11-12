#pragma once
#include "MetadataPaths.hpp"
#include <filesystem>
#include <stdexcept>
#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>
#endif
namespace dynsyslab::workflow::detail {
// Cross-process writer exclusion. The OS releases the lock if a writer crashes.
class ResearchCatalogueLock {
  public:
    explicit ResearchCatalogueLock(const std::filesystem::path &directory) {
        const auto ioDirectory = metadataIoPath(directory);
        std::filesystem::create_directories(ioDirectory);
        const auto file = ioDirectory / ".writer-lock";
#ifdef _WIN32
        handle_ = CreateFileW(file.c_str(), GENERIC_READ | GENERIC_WRITE,
                              FILE_SHARE_READ | FILE_SHARE_WRITE, nullptr, OPEN_ALWAYS,
                              FILE_ATTRIBUTE_NORMAL, nullptr);
        if (handle_ == INVALID_HANDLE_VALUE)
            throw std::runtime_error("Cannot open bibliography writer lock");
        if (!LockFileEx(handle_, LOCKFILE_EXCLUSIVE_LOCK | LOCKFILE_FAIL_IMMEDIATELY, 0, 1, 0,
                        &overlap_)) {
            CloseHandle(handle_);
            handle_ = INVALID_HANDLE_VALUE;
            throw std::runtime_error("Bibliography has another active writer; retry the addition");
        }
#else
        handle_ = ::open(file.c_str(), O_RDWR | O_CREAT | O_CLOEXEC, 0600);
        if (handle_ < 0)
            throw std::runtime_error("Cannot open bibliography writer lock");
        if (::flock(handle_, LOCK_EX | LOCK_NB) != 0) {
            ::close(handle_);
            handle_ = -1;
            throw std::runtime_error("Bibliography has another active writer; retry the addition");
        }
#endif
    }
    ResearchCatalogueLock(const ResearchCatalogueLock &) = delete;
    ResearchCatalogueLock &operator=(const ResearchCatalogueLock &) = delete;
    ~ResearchCatalogueLock() {
#ifdef _WIN32
        if (handle_ != INVALID_HANDLE_VALUE) {
            UnlockFileEx(handle_, 0, 1, 0, &overlap_);
            CloseHandle(handle_);
        }
#else
        if (handle_ >= 0) {
            ::flock(handle_, LOCK_UN);
            ::close(handle_);
        }
#endif
    }

  private:
#ifdef _WIN32
    HANDLE handle_{INVALID_HANDLE_VALUE};
    OVERLAPPED overlap_{};
#else
    int handle_{-1};
#endif
};
} // namespace dynsyslab::workflow::detail
