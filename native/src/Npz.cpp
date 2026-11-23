#include <dynsyslab/core/LongPath.hpp>
#include <dynsyslab/core/Archive.hpp>
#include <dynsyslab/core/BuildInfo.hpp>
#include <dynsyslab/core/Json.hpp>

#include <archive.h>
#include <archive_entry.h>

#include <algorithm>
#include <chrono>
#include <cstring>
#include <sstream>
#include <stdexcept>

#ifdef _WIN32
#include <windows.h>
#endif

namespace dynsyslab::core {
namespace {

[[noreturn]] void archiveError(archive* handle, const char* context) {
    const char* detail = archive_error_string(handle);
    throw std::runtime_error(std::string(context) + (detail ? ": " + std::string(detail) : ""));
}

ArrayType arrayType(const std::string& descriptor) {
    if (descriptor == "<f8" || descriptor == "=f8" || descriptor == "|f8") return ArrayType::float64;
    if (descriptor == "<f4" || descriptor == "=f4" || descriptor == "|f4") return ArrayType::float32;
    if (descriptor == "<i8" || descriptor == "=i8" || descriptor == "|i8") return ArrayType::int64;
    if (descriptor == "<i4" || descriptor == "=i4" || descriptor == "|i4") return ArrayType::int32;
    if (descriptor == "<u8" || descriptor == "=u8" || descriptor == "|u8") return ArrayType::uint64;
    if (descriptor == "<u4" || descriptor == "=u4" || descriptor == "|u4") return ArrayType::uint32;
    if (descriptor == "|u1") return ArrayType::uint8;
    if (descriptor == "|b1") return ArrayType::boolean;
    if (descriptor == "<c16" || descriptor == "=c16") return ArrayType::complex128;
    if (descriptor.rfind("<U", 0) == 0 || descriptor.rfind("=U", 0) == 0) return ArrayType::unicode;
    return ArrayType::unknown;
}

std::string field(const std::string& header, const std::string& name) {
    const auto key = header.find("'" + name + "'");
    const auto colon = key == std::string::npos ? key : header.find(':', key);
    const auto quote = colon == std::string::npos ? colon : header.find_first_of("'\"", colon);
    const auto end = quote == std::string::npos ? quote : header.find(header[quote], quote + 1);
    return end == std::string::npos ? std::string{} : header.substr(quote + 1, end - quote - 1);
}

std::vector<std::size_t> shape(const std::string& header) {
    const auto key = header.find("'shape'");
    const auto open = key == std::string::npos ? key : header.find('(', key);
    const auto close = open == std::string::npos ? open : header.find(')', open);
    if (close == std::string::npos) throw std::runtime_error("invalid NPY shape header");
    std::vector<std::size_t> result;
    std::stringstream values(header.substr(open + 1, close - open - 1));
    std::string token;
    while (std::getline(values, token, ',')) {
        const auto first = token.find_first_not_of(" \t");
        if (first != std::string::npos) result.push_back(static_cast<std::size_t>(std::stoull(token.substr(first))));
    }
    return result;
}

void makeRowMajor(NpyArray& array) {
    if (array.shape.size()<2) return;
    const std::size_t count=array.elementCount();
    if (!count || array.data.size()%count)
        throw std::runtime_error("Fortran-ordered NPY byte count is invalid");
    const std::size_t width=array.data.size()/count;
    const auto source=array.data;
    std::vector<std::size_t> cStride(array.shape.size(),1);
    for (std::size_t axis=array.shape.size()-1;axis>0;--axis)
        cStride[axis-1]=cStride[axis]*array.shape[axis];
    for (std::size_t destination=0;destination<count;++destination) {
        std::size_t remainder=destination,origin=0,fortranStride=1;
        for (std::size_t axis=0;axis<array.shape.size();++axis) {
            const std::size_t index=remainder/cStride[axis];
            remainder%=cStride[axis]; origin+=index*fortranStride;
            fortranStride*=array.shape[axis];
        }
        std::memcpy(array.data.data()+destination*width,
                    source.data()+origin*width,width);
    }
}

NpyArray parseNpy(const std::vector<std::uint8_t>& bytes) {
    constexpr unsigned char magic[]{0x93, 'N', 'U', 'M', 'P', 'Y'};
    if (bytes.size() < 10 || !std::equal(std::begin(magic), std::end(magic), bytes.begin()))
        throw std::runtime_error("invalid NPY magic");
    const unsigned major = bytes[6];
    std::size_t headerLength = 0, offset = 0;
    if (major == 1) {
        headerLength = bytes[8] | (static_cast<std::size_t>(bytes[9]) << 8); offset = 10;
    } else if (major == 2 || major == 3) {
        if (bytes.size() < 12) throw std::runtime_error("truncated NPY header");
        headerLength = bytes[8] | (static_cast<std::size_t>(bytes[9]) << 8)
            | (static_cast<std::size_t>(bytes[10]) << 16) | (static_cast<std::size_t>(bytes[11]) << 24);
        offset = 12;
    } else throw std::runtime_error("unsupported NPY version");
    if (offset + headerLength > bytes.size()) throw std::runtime_error("truncated NPY header");
    const std::string header(reinterpret_cast<const char*>(bytes.data() + offset), headerLength);
    const bool fortran=header.find("'fortran_order': True") != std::string::npos;
    NpyArray result;
    result.descriptor = field(header, "descr");
    result.type = arrayType(result.descriptor);
    if (result.type == ArrayType::unknown)
        throw std::runtime_error("unsupported NPY dtype: " + result.descriptor);
    result.shape = shape(header);
    result.data.assign(bytes.begin() + static_cast<std::ptrdiff_t>(offset + headerLength), bytes.end());
    if (fortran) makeRowMajor(result);
    return result;
}

std::vector<std::uint8_t> encodeNpy(const NpyArray& array) {
    std::ostringstream header;
    header << "{'descr': '" << array.descriptor << "', 'fortran_order': False, 'shape': (";
    for (std::size_t index = 0; index < array.shape.size(); ++index) {
        if (index) header << ", ";
        header << array.shape[index];
    }
    if (array.shape.size() == 1) header << ',';
    header << "), }";
    std::string text = header.str();
    const std::size_t prefix = 10;
    const std::size_t padding = (64 - ((prefix + text.size() + 1) % 64)) % 64;
    text.append(padding, ' '); text.push_back('\n');
    if (text.size() > 65535) throw std::runtime_error("NPY v1 header is too large");
    std::vector<std::uint8_t> result;
    result.reserve(prefix + text.size() + array.data.size());
    const std::uint8_t magic[]{0x93, 'N', 'U', 'M', 'P', 'Y', 1, 0,
        static_cast<std::uint8_t>(text.size() & 0xff), static_cast<std::uint8_t>((text.size() >> 8) & 0xff)};
    result.insert(result.end(), std::begin(magic), std::end(magic));
    result.insert(result.end(), text.begin(), text.end());
    result.insert(result.end(), array.data.begin(), array.data.end());
    return result;
}

std::string arrayName(const char* path) {
    std::filesystem::path value = path ? std::filesystem::u8path(path) : std::filesystem::path{};
    return value.extension() == L".npy" ? value.stem().u8string() : std::string{};
}

}  // namespace

NpzRecord readNpz(const std::filesystem::path& path) {
    archive* handle = archive_read_new();
    if (!handle) throw std::runtime_error("unable to allocate archive reader");
    archive_read_support_format_zip(handle);
    archive_read_support_filter_all(handle);
#ifdef _WIN32
    const int opened = archive_read_open_filename_w(handle, path.c_str(), 10240);
#else
    const int opened = archive_read_open_filename(handle, path.c_str(), 10240);
#endif
    if (opened != ARCHIVE_OK) { const std::string message = archive_error_string(handle); archive_read_free(handle); throw std::runtime_error(message); }
    NpzRecord result;
    archive_entry* entry = nullptr;
    while (archive_read_next_header(handle, &entry) == ARCHIVE_OK) {
        const std::string name = arrayName(archive_entry_pathname(entry));
        const auto size = archive_entry_size(entry);
        if (name.empty() || size < 0) { archive_read_data_skip(handle); continue; }
        std::vector<std::uint8_t> bytes(static_cast<std::size_t>(size));
        std::size_t offset = 0;
        while (offset < bytes.size()) {
            const auto count = archive_read_data(handle, bytes.data() + offset, bytes.size() - offset);
            if (count < 0) { archive_read_free(handle); throw std::runtime_error("unable to read NPZ entry"); }
            if (count == 0) break;
            offset += static_cast<std::size_t>(count);
        }
        bytes.resize(offset);
        result.arrays.emplace(name, parseNpy(bytes));
    }
    archive_read_free(handle);
    return result;
}

bool NpzIndex::contains(const std::string& name) const {
    return std::find(arrays.begin(),arrays.end(),name)!=arrays.end();
}

NpzIndex readNpzIndex(const std::filesystem::path& path) {
    archive* handle=archive_read_new();
    if (!handle) throw std::runtime_error("unable to allocate archive reader");
    archive_read_support_format_zip(handle); archive_read_support_filter_all(handle);
#ifdef _WIN32
    const int opened=archive_read_open_filename_w(handle,path.c_str(),10240);
#else
    const int opened=archive_read_open_filename(handle,path.c_str(),10240);
#endif
    if (opened!=ARCHIVE_OK) {
        const std::string message=archive_error_string(handle);
        archive_read_free(handle); throw std::runtime_error(message);
    }
    NpzIndex result; archive_entry* entry=nullptr;
    while (archive_read_next_header(handle,&entry)==ARCHIVE_OK) {
        const std::string name=arrayName(archive_entry_pathname(entry));
        const auto size=archive_entry_size(entry);
        if (name.empty() || size<0) { archive_read_data_skip(handle); continue; }
        result.arrays.push_back(name);
        if (name!="metadata_json") { archive_read_data_skip(handle); continue; }
        std::vector<std::uint8_t> bytes(static_cast<std::size_t>(size));
        std::size_t offset=0;
        while (offset<bytes.size()) {
            const auto count=archive_read_data(handle,bytes.data()+offset,bytes.size()-offset);
            if (count<0) { archive_read_free(handle);
                throw std::runtime_error("unable to read NPZ metadata entry"); }
            if (count==0) break;
            offset+=static_cast<std::size_t>(count);
        }
        bytes.resize(offset); result.metadataJson=parseNpy(bytes).string();
    }
    archive_read_free(handle); return result;
}

// Adds the active run provenance to metadata_json; numerical arrays are never touched.
static NpyArray withRunProvenance(const NpyArray& array) {
    const auto& provenance = activeRunProvenance();
    if (provenance.empty()) return array;
    auto text = array.string();
    const auto end = text.find_last_not_of(" \r\n\t");
    if (end == std::string::npos || text[end] != '}' || text.front() != '{') return array;
    if (text.find("\"run_provenance\"") != std::string::npos) {  // copied from a parent archive
        auto metadata = json::parse(text);
        metadata.at("run_provenance") = json::parse(provenance);
        return NpyArray::fromString(json::dump(metadata));
    }
    text.resize(end);
    if (text.find_first_not_of(" \r\n\t", 1) != std::string::npos) text += ",";
    return NpyArray::fromString(text + "\"run_provenance\":" + provenance + "}");
}

void writeTextZip(const std::filesystem::path& path, const std::vector<std::pair<std::string, std::string>>& entries) {
    std::filesystem::create_directories(extendedLengthPath(path.parent_path()));
    archive* handle = archive_write_new();
    if (!handle) throw std::runtime_error("unable to allocate archive writer");
    archive_write_set_format_zip(handle);
    archive_write_set_options(handle, "zip:compression=deflate");
#ifdef _WIN32
    const int opened = archive_write_open_filename_w(handle, path.c_str());
#else
    const int opened = archive_write_open_filename(handle, path.c_str());
#endif
    if (opened != ARCHIVE_OK) archiveError(handle, "unable to create ZIP");
    for (const auto& [name, text] : entries) {
        archive_entry* entry = archive_entry_new();
        archive_entry_set_pathname(entry, name.c_str());
        archive_entry_set_filetype(entry, AE_IFREG);
        archive_entry_set_perm(entry, 0644);
        archive_entry_set_size(entry, static_cast<la_int64_t>(text.size()));
        archive_entry_set_mtime(entry, 315532800, 0);  // 1980-01-01, the ZIP epoch
        if (archive_write_header(handle, entry) != ARCHIVE_OK
            || archive_write_data(handle, text.data(), text.size()) < 0) {
            archive_entry_free(entry); archiveError(handle, "unable to write ZIP entry");
        }
        archive_entry_free(entry);
    }
    archive_write_close(handle);
    archive_write_free(handle);
}

static void writeNpzAtomicImpl(const std::filesystem::path& path, const NpzRecord& record, bool replaceExisting) {
    std::filesystem::create_directories(extendedLengthPath(path.parent_path()));
    const auto token = std::chrono::high_resolution_clock::now().time_since_epoch().count();
    const auto temporary = path.parent_path() / (path.filename().wstring() + L"."
        + std::to_wstring(token) + L".partial");
    archive* handle = archive_write_new();
    if (!handle) throw std::runtime_error("unable to allocate archive writer");
    archive_write_set_format_zip(handle);
    archive_write_set_options(handle, "zip:compression=deflate");
#ifdef _WIN32
    const int opened = archive_write_open_filename_w(handle, temporary.c_str());
#else
    const int opened = archive_write_open_filename(handle, temporary.c_str());
#endif
    if (opened != ARCHIVE_OK) archiveError(handle, "unable to create NPZ archive");
    try {
        for (const auto& [name, array] : record.arrays) {
            const auto bytes = encodeNpy(name == "metadata_json" ? withRunProvenance(array) : array);
            archive_entry* entry = archive_entry_new();
            archive_entry_set_pathname(entry, (name + ".npy").c_str());
            archive_entry_set_filetype(entry, AE_IFREG);
            archive_entry_set_perm(entry, 0644);
            archive_entry_set_size(entry, static_cast<la_int64_t>(bytes.size()));
            if (archive_write_header(handle, entry) != ARCHIVE_OK) {
                archive_entry_free(entry); archiveError(handle, "unable to write NPZ header");
            }
            if (archive_write_data(handle, bytes.data(), bytes.size()) < 0) {
                archive_entry_free(entry); archiveError(handle, "unable to write NPZ entry");
            }
            archive_entry_free(entry);
        }
        archive_write_close(handle);
        archive_write_free(handle);
        handle = nullptr;
#ifdef _WIN32
        // libarchive supports long paths; the final Win32 rename must use them too.
        const auto extended = [](const std::filesystem::path& source) { return extendedLengthPath(source).wstring(); };
        if (!MoveFileExW(extended(temporary).c_str(), extended(path).c_str(), (replaceExisting ? MOVEFILE_REPLACE_EXISTING : 0) | MOVEFILE_WRITE_THROUGH))
            throw std::runtime_error("unable to atomically publish NPZ archive");
#else
        if (replaceExisting) std::filesystem::rename(temporary, path);
        else {
            std::filesystem::create_hard_link(temporary, path);
            std::error_code ignored; std::filesystem::remove(temporary, ignored);
        }
#endif
    } catch (...) {
        if (handle) { archive_write_close(handle); archive_write_free(handle); }
        std::error_code ignored; std::filesystem::remove(temporary, ignored);
        throw;
    }
}

void writeNpzAtomic(const std::filesystem::path& path, const NpzRecord& record) {
    writeNpzAtomicImpl(path,record,true);
}
void writeNpzAtomicNew(const std::filesystem::path& path, const NpzRecord& record) {
    writeNpzAtomicImpl(path,record,false);
}

}  // namespace dynsyslab::core
