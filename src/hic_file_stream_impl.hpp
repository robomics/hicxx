// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef STRAW_USE_ZLIBNG
#include <zlib-ng.h>
#else
#include <zlib.h>
#endif

#include <fmt/format.h>

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <filestream/filestream.hpp>
#include <fstream>
#include <ios>
#include <stdexcept>
#include <string>
#include <utility>

#include "straw/internal/common.h"

namespace internal {

constexpr HiCHeader::operator bool() const noexcept { return masterIndexOffset >= 0; }

inline bool HiCHeader::operator==(const HiCHeader &other) const noexcept {
    return url == other.url && masterIndexOffset == other.masterIndexOffset;
}

inline bool HiCHeader::operator!=(const HiCHeader &other) const noexcept {
    return !(*this == other);
}

inline std::int32_t HiCHeader::nChromosomes() const noexcept {
    return static_cast<std::int32_t>(chromosomes.size());
}

inline std::int32_t HiCHeader::nResolutions() const noexcept {
    return static_cast<std::int32_t>(resolutions.size());
}

inline const chromosome &HiCHeader::getChromosome(std::int32_t id) const noexcept {
    // Chromosomes are sorted by id, so we can use simple arithmetic on iterators to find the
    // chromosome with the given id
    assert(id < chromosomes.size());

    const auto it = std::next(chromosomes.begin(), id);
    assert(it->second.index == id);
    return it->second;
}

constexpr HiCFooterMetadata::operator bool() const noexcept { return fileOffset >= 0; }

inline bool HiCFooterMetadata::operator==(const HiCFooterMetadata &other) const noexcept {
    return url == other.url && matrixType == other.matrixType &&
           normalization == other.normalization && unit == other.unit &&
           resolution == other.resolution && chrom1 == other.chrom1 && chrom2 == other.chrom2;
}

inline bool HiCFooterMetadata::operator!=(const HiCFooterMetadata &other) const noexcept {
    return !(*this == other);
}

inline HiCFooter::HiCFooter(HiCFooterMetadata metadata_) noexcept
    : _metadata(std::move(metadata_)) {}

constexpr HiCFooter::operator bool() const noexcept { return !metadata(); }
inline bool HiCFooter::operator==(const HiCFooter &other) const noexcept {
    return metadata() == other.metadata();
}
inline bool HiCFooter::operator!=(const HiCFooter &other) const noexcept {
    return !(*this == other);
}
constexpr const HiCFooterMetadata &HiCFooter::metadata() const noexcept { return _metadata; }
constexpr HiCFooterMetadata &HiCFooter::metadata() noexcept { return _metadata; }
constexpr const std::string &HiCFooter::url() const noexcept { return metadata().url; }
constexpr MatrixType HiCFooter::matrixType() const noexcept { return metadata().matrixType; }
constexpr NormalizationMethod HiCFooter::normalization() const noexcept {
    return metadata().normalization;
}
constexpr MatrixUnit HiCFooter::unit() const noexcept { return metadata().unit; }
constexpr std::int32_t HiCFooter::resolution() const noexcept { return metadata().resolution; }
constexpr const chromosome &HiCFooter::chrom1() const noexcept { return metadata().chrom1; }
constexpr const chromosome &HiCFooter::chrom2() const noexcept { return metadata().chrom2; }
constexpr std::int64_t HiCFooter::fileOffset() const noexcept { return metadata().fileOffset; }

constexpr const std::vector<double> &HiCFooter::expectedValues() const noexcept {
    return _expectedValues;
}

constexpr const std::vector<double> &HiCFooter::c1Norm() const noexcept { return _c1Norm; }

constexpr const std::vector<double> &HiCFooter::c2Norm() const noexcept {
    if (chrom1() == chrom2()) {
        return _c1Norm;
    }
    return _c2Norm;
}

constexpr std::vector<double> &HiCFooter::expectedValues() noexcept { return _expectedValues; }

constexpr std::vector<double> &HiCFooter::c1Norm() noexcept { return _c1Norm; }

constexpr std::vector<double> &HiCFooter::c2Norm() noexcept {
    if (chrom1() == chrom2()) {
        return _c1Norm;
    }
    return _c2Norm;
}

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline T BinaryBuffer::read() {
    static_assert(sizeof(char) == 1, "");
    assert(i < buffer.size());
    T x{};

    std::memcpy(static_cast<void *>(&x), buffer.data() + i, sizeof(T));
    i += sizeof(T);
    return x;
}

inline HiCFileStream::HiCFileStream(std::string url)
    : _fs(std::make_shared<filestream::FileStream>(HiCFileStream::openStream(std::move(url)))),
      _header(std::make_shared<const HiCHeader>(HiCFileStream::readHeader(*_fs))) {}

inline filestream::FileStream HiCFileStream::openStream(std::string url) {
    const auto isRemoteFile = StartsWith(url, "http");
    try {
        try {
            if (isRemoteFile) {
                constexpr std::size_t defaulChunkSize = 64 * 1024;
                return filestream::FileStream::Remote(url, defaulChunkSize, "straw");
            }
            return filestream::FileStream::Local(url);

        } catch (const std::system_error &e) {
            throw std::runtime_error(std::string(e.what()) + ": " + e.code().message());
        }
    } catch (const std::exception &e) {
        const auto file_type = isRemoteFile ? "remote file" : "file";
        throw std::runtime_error(std::string("Failed to open ") + file_type + " " + url + ": " +
                                 e.what());
    }
}

inline const std::string &HiCFileStream::url() const noexcept { return _fs->url(); }
inline const HiCHeader &HiCFileStream::header() const noexcept { return *_header; }

inline bool HiCFileStream::isLocal() const noexcept { return _fs->is_local(); }

inline bool HiCFileStream::isRemote() const noexcept { return _fs->is_remote(); }

inline std::int32_t HiCFileStream::version() const noexcept {
    assert(_header->version != -1);
    return _header->version;
}

inline void HiCFileStream::discardExpectedVector(std::int64_t nValues) {
    const auto elementSize = version() > 8 ? sizeof(float) : sizeof(double);
    _fs->seekg(nValues * elementSize, std::ios::cur);
}

inline std::vector<double> HiCFileStream::readExpectedVector(std::int64_t nValues) {
    std::vector<double> initialExpectedValues(nValues);
    if (version() > 8) {
        std::vector<float> tmpbuff(nValues);
        _fs->read(tmpbuff);
        std::transform(tmpbuff.begin(), tmpbuff.end(), initialExpectedValues.begin(),
                       [](float n) { return static_cast<double>(n); });
    } else {
        _fs->read(initialExpectedValues);
    }

    // This seems to be copying initialValues into finalResult at the moment
    // std::int32_t window = 5000000 / resolution;
    // rollingMedian(initialExpectedValues, _expectedValues, window);
    return initialExpectedValues;
}

inline std::vector<double> HiCFileStream::readNormalizationFactors(std::int32_t wantedChrom) {
    const auto nFactors = _fs->read<std::int32_t>();
    std::vector<double> normFactors{};
    auto readFactor = [this]() -> double {
        if (version() > 8) {
            return _fs->read<float>();
        }
        return _fs->read<double>();
    };

    for (auto i = 0; i < nFactors; ++i) {
        const auto foundChrom = _fs->read<std::int32_t>();
        const auto v = readFactor();
        if (foundChrom == wantedChrom) {
            normFactors.push_back(v);
        }
    }
    return normFactors;
}

inline void HiCFileStream::applyNormalizationFactors(std::vector<double> &expectedValues,
                                                     const std::vector<double> &normFactors) {
    if (normFactors.empty() || expectedValues.empty()) {
        return;
    }
    for (const auto factor : normFactors) {
        std::transform(expectedValues.begin(), expectedValues.end(), expectedValues.begin(),
                       [&](auto n) { return n / factor; });
    }
}
inline std::vector<double> HiCFileStream::readNormalizationVector(indexEntry cNormEntry) {
    _fs->seekg(cNormEntry.position);
    const auto numValues = static_cast<std::size_t>(readNValues());

    std::vector<double> buffer(numValues);
    if (version() > 8) {
        std::vector<float> tmpbuffer(numValues);
        _fs->read(tmpbuffer);
        std::transform(tmpbuffer.begin(), tmpbuffer.end(), buffer.begin(),
                       [](float n) { return static_cast<double>(n); });

    } else {
        _fs->read(buffer);
    }
    return buffer;
}

inline void HiCFileStream::discardNormalizationFactors(std::int32_t wantedChrom) {
    (void)readNormalizationFactors(wantedChrom);
}

inline MatrixType HiCFileStream::readMatrixType(filestream::FileStream &fs, std::string &buff) {
    fs.getline(buff, '\0');
    return ParseMatrixTypeStr(buff);
}

inline NormalizationMethod HiCFileStream::readNormalizationMethod(filestream::FileStream &fs,
                                                                  std::string &buff) {
    fs.getline(buff, '\0');
    return ParseNormStr(buff);
}

inline MatrixUnit HiCFileStream::readMatrixUnit(filestream::FileStream &fs, std::string &buff) {
    fs.getline(buff, '\0');
    return ParseUnitStr(buff);
}

inline MatrixType HiCFileStream::readMatrixType() {
    return HiCFileStream::readMatrixType(*_fs, _strbuff);
}

inline NormalizationMethod HiCFileStream::readNormalizationMethod() {
    return HiCFileStream::readNormalizationMethod(*_fs, _strbuff);
}

inline MatrixUnit HiCFileStream::readMatrixUnit() {
    return HiCFileStream::readMatrixUnit(*_fs, _strbuff);
}

inline std::int64_t HiCFileStream::readNValues() {
    if (version() > 8) {
        return _fs->read<std::int64_t>();
    }
    return _fs->read<std::int32_t>();
}

inline bool HiCFileStream::checkMagicString(filestream::FileStream &fs) {
    return fs.getline('\0') == "HIC";
}

inline std::int64_t HiCFileStream::masterOffset() const noexcept {
    return _header->masterIndexOffset;
}

inline const char *strawZError(int status) {
#ifdef STRAW_USE_ZLIBNG
    return zng_zError(status);
#else
    return zError(status);
#endif
}

inline auto HiCFileStream::initZStream() -> ZStream {
#ifdef STRAW_USE_ZLIBNG
    ZStream zs(new zng_stream, &zng_inflateEnd);
#else
    ZStream zs(new z_stream, &inflateEnd);
#endif
    zs->zalloc = Z_NULL;
    zs->zfree = Z_NULL;
    zs->opaque = Z_NULL;
    // Signal no input data is provided for initialization
    zs->avail_in = 0;
    zs->next_in = Z_NULL;

#ifdef STRAW_USE_ZLIBNG
    const auto status = zng_inflateInit(zs.get());
#else
    const auto status = inflateInit(zs.get());
#endif
    if (status != Z_OK) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("failed to initialize zlib decompression stream: {}"), strawZError(status)));
    }

    return zs;
}

inline void HiCFileStream::readBlockMap(std::int64_t fileOffset, const chromosome &chrom1,
                                        const chromosome &chrom2, MatrixUnit wantedUnit,
                                        std::int32_t wantedResolution, BlockMap &buffer) {
    _fs->seekg(fileOffset);
    auto &blockMap = buffer.blocks;
    blockMap.clear();

    const auto c1i = _fs->read<std::int32_t>();
    const auto c2i = _fs->read<std::int32_t>();
    const auto numResolutions = _fs->read<std::int32_t>();

    assert(c1i == chrom1.index);
    assert(c2i == chrom2.index);

    for (std::int32_t i = 0; i < numResolutions; ++i) {
        const auto foundUnit = readMatrixUnit();
        (void)_fs->read<std::int32_t>();  // oldIndex
        const auto sumCount = _fs->read<float>();
        (void)_fs->read<float>();  // occupiedCellCount
        (void)_fs->read<float>();  // stdDev
        (void)_fs->read<float>();  // percent95

        const auto foundResolution = _fs->read<std::int32_t>();
        const auto blockBinCount = _fs->read<std::int32_t>();
        const auto blockColumnCount = _fs->read<std::int32_t>();

        const auto nBlocks = static_cast<std::size_t>(_fs->read<std::int32_t>());

        if (wantedUnit == foundUnit && wantedResolution == foundResolution) {
            for (std::size_t j = 0; j < nBlocks; ++j) {
                const auto key = _fs->read<std::int32_t>();
                indexEntry index{_fs->read<std::int64_t>(), _fs->read<std::int32_t>()};
                assert(index.position + index.size < _fs->size());
                blockMap.emplace(key, std::move(index));
            }
            buffer.blockBinCount = blockBinCount;
            buffer.blockColumnCount = blockColumnCount;
            buffer.sumCount = static_cast<double>(sumCount);
            return;
        }

        constexpr auto blockSize = sizeof(int32_t) + sizeof(int64_t) + sizeof(int32_t);
        _fs->seekg(nBlocks * blockSize, std::ios::cur);
    }

    throw std::runtime_error(
        fmt::format(FMT_STRING("Unable to find block map for unit {} and resolution {}"),
                    wantedUnit, wantedResolution));
}

inline bool HiCFileStream::checkMagicString() { return checkMagicString(*_fs); }

// reads the header, storing the positions of the normalization vectors and returning the
// masterIndexPosition pointer
inline HiCHeader HiCFileStream::readHeader(filestream::FileStream &fs) {
    if (!checkMagicString(fs)) {
        throw std::runtime_error("Hi-C magic string is missing. " + fs.url() +
                                 " does not appear to be a hic file");
    }

    HiCHeader header{fs.url()};

    fs.read(header.version);
    if (header.version < 6) {
        throw std::runtime_error("Version " + std::to_string(header.version) +
                                 " no longer supported");
    }
    fs.read(header.masterIndexOffset);
    if (header.masterIndexOffset < 0 || header.masterIndexOffset >= fs.size()) {
        throw std::runtime_error(
            "Invalid masterOffset index offset. Expected offset between 0 and " +
            std::to_string(fs.size()) + ", found " + std::to_string(header.masterIndexOffset));
    }

    fs.getline(header.genomeID, '\0');
    if (header.genomeID.empty()) {
        header.genomeID = "unknown";
    }

    if (header.version > 8) {
        fs.read(header.nviPosition);
        fs.read(header.nviLength);
    }

    const auto nAttributes = fs.read<std::int32_t>();

    // reading and ignoring attribute-value dictionary
    for (std::int32_t i = 0; i < nAttributes; i++) {
        (void)fs.getline('\0');  // key
        (void)fs.getline('\0');  // value
    }

    // Read chromosomes
    const auto numChromosomes = fs.read<std::int32_t>();
    std::generate_n(std::inserter(header.chromosomes, header.chromosomes.begin()), numChromosomes,
                    [&]() {
                        chromosome chrom{};
                        chrom.index = static_cast<std::int32_t>(header.chromosomes.size());
                        fs.getline(chrom.name, '\0');
                        if (header.version > 8) {
                            fs.read(chrom.length);
                        } else {
                            chrom.length = static_cast<std::int64_t>(fs.read<std::int32_t>());
                        }

                        return std::make_pair(chrom.name, chrom);
                    });

    if (header.chromosomes.empty()) {
        throw std::runtime_error("unable to read chromosomes");
    }

    // Read resolutions
    const auto numResolutions = fs.read<std::int32_t>();
    header.resolutions.resize(numResolutions);
    std::generate(header.resolutions.begin(), header.resolutions.end(), [&]() {
        const auto res = fs.read<std::int32_t>();
        assert(res > 0);
        return res;
    });

    if (header.resolutions.empty()) {
        throw std::runtime_error("unable to read the list of available resolutions");
    }

    return header;
}

inline void HiCFileStream::readAndInflate(indexEntry idx, std::string &plainTextBuffer) {
    try {
        // _strbuff is used to store compressed data
        // plainTextBuffer is used to store decompressed data
        assert(_zlibstream);
        assert(idx.size > 0);

        _fs->seekg(idx.position);
        _fs->read(_strbuff, idx.size);

        _zlibstream->avail_in = static_cast<uInt>(_strbuff.size());
        _zlibstream->next_in = reinterpret_cast<Bytef *>(&*(_strbuff.begin()));

#ifdef STRAW_USE_ZLIBNG
        auto status = zng_inflateReset(_zlibstream.get());
#else
        auto status = inflateReset(_zlibstream.get());
#endif
        if (status != Z_OK) {
            plainTextBuffer.clear();
            throw std::runtime_error(strawZError(status));
        }

        plainTextBuffer.reserve(idx.size * 3);
        plainTextBuffer.resize(plainTextBuffer.capacity());
        auto current_size = 0;

        while (true) {
            _zlibstream->avail_out = static_cast<uInt>(plainTextBuffer.size() - current_size);
            _zlibstream->next_out =
                reinterpret_cast<Bytef *>(&*(plainTextBuffer.begin() + current_size));

#ifdef STRAW_USE_ZLIBNG
            status = zng_inflate(_zlibstream.get(), Z_NO_FLUSH);
#else
            status = inflate(_zlibstream.get(), Z_NO_FLUSH);
#endif
            if (status == Z_STREAM_END) {
                // assert(_zlibstream->avail_in == 0);
                break;
            }
            if (status != Z_OK) {
                plainTextBuffer.clear();
                throw std::runtime_error(strawZError(status));
            }

            current_size = plainTextBuffer.size();
            plainTextBuffer.resize(current_size + idx.size);
        }

        // assert(plainTextBuffer.size() >= _zlibstream->total_out);
        plainTextBuffer.resize(static_cast<std::size_t>(_zlibstream->total_out));
    } catch (const std::exception &e) {
        throw std::runtime_error(fmt::format(FMT_STRING("failed to decompress block at pos {}: {}"),
                                             idx.position, e.what()));
    }
}

inline HiCFooter HiCFileStream::readFooter(const std::int32_t chromId1, const std::int32_t chromId2,
                                           const MatrixType wantedMatrixType,
                                           const NormalizationMethod wantedNorm,
                                           const MatrixUnit wantedUnit,
                                           const std::int32_t wantedResolution) {
    assert(chromId1 <= chromId2);
    assert(std::find(_header->resolutions.begin(), _header->resolutions.end(), wantedResolution) !=
           _header->resolutions.end());

    using MT = MatrixType;
    using NM = NormalizationMethod;

    // clang-format off
    HiCFooter footer{
        HiCFooterMetadata{_fs->url(),
                          wantedMatrixType,
                          wantedNorm,
                          wantedUnit,
                          wantedResolution,
                          _header->getChromosome(chromId1),
                         _header->getChromosome(chromId2)}
        };
    // clang-format on

    auto &metadata = footer.metadata();
    auto &expectedValues = footer.expectedValues();
    auto &c1Norm = footer.c1Norm();
    auto &c2Norm = footer.c2Norm();

    const auto key = std::to_string(chromId1) + "_" + std::to_string(chromId2);

    _fs->seekg(masterOffset());
    (void)readNValues();  // nBytes

    auto nEntries = _fs->read<std::int32_t>();
    for (int i = 0; i < nEntries; i++) {
        const auto strbuff = _fs->getline('\0');
        const auto fpos = _fs->read<std::int64_t>();
        (void)_fs->read<std::int32_t>();  // sizeInBytes
        if (strbuff == key) {
            metadata.fileOffset = fpos;
        }
    }
    if (metadata.fileOffset == -1) {
        throw std::runtime_error("File doesn't have the given chr_chr map " + key);
    }

    if ((wantedMatrixType == MT::observed && wantedNorm == NM::NONE) ||
        ((wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) &&
         wantedNorm == NM::NONE && chromId1 != chromId2)) {
        return footer;  // no need to read wantedNorm vector index
    }

    // read in and ignore expected value maps; don't store; reading these to
    // get to wantedNorm vector index
    auto nExpectedValues = _fs->read<std::int32_t>();
    for (std::int32_t i = 0; i < nExpectedValues; i++) {
        const auto foundUnit = readMatrixUnit();
        const auto foundResolution = _fs->read<std::int32_t>();
        const auto nValues = readNValues();

        bool store = chromId1 == chromId2 &&
                     (wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) &&
                     wantedNorm == NM::NONE && foundUnit == wantedUnit &&
                     foundResolution == wantedResolution;

        if (store) {
            expectedValues = readExpectedVector(nValues);
            const auto normFactors = readNormalizationFactors(chromId1);
            applyNormalizationFactors(expectedValues, normFactors);

        } else {
            discardExpectedVector(nValues);
            discardNormalizationFactors(chromId1);
        }
    }

    if (chromId1 == chromId2 && (wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) &&
        wantedNorm == NM::NONE) {
        if (expectedValues.empty()) {
            throw std::runtime_error(fmt::format(
                FMT_STRING(
                    "File did not contain expected values vectors for unit {} at resolution {}"),
                wantedUnit, wantedResolution));
        }
        return footer;
    }

    nExpectedValues = _fs->read<std::int32_t>();
    for (std::int32_t i = 0; i < nExpectedValues; i++) {
        const auto foundNorm = readNormalizationMethod();
        const auto foundUnit = readMatrixUnit();
        const auto foundResolution = _fs->read<std::int32_t>();

        const auto nValues = readNValues();
        bool store = chromId1 == chromId2 &&
                     (wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) &&
                     foundNorm == wantedNorm && foundUnit == wantedUnit &&
                     foundResolution == wantedResolution;

        if (store) {
            expectedValues = readExpectedVector(nValues);
            const auto normFactors = readNormalizationFactors(chromId1);
            applyNormalizationFactors(expectedValues, normFactors);
        } else {
            discardExpectedVector(nValues);
            discardNormalizationFactors(chromId1);
        }
    }

    if (chromId1 == chromId2 && (wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) &&
        wantedNorm != NM::NONE) {
        if (expectedValues.empty()) {
            throw std::runtime_error(fmt::format(
                FMT_STRING(
                    "File did not contain normalized expected values for unit {} at resolution {}"),
                wantedUnit, wantedResolution));
        }
    }

    // Index of normalization vectors
    nEntries = _fs->read<std::int32_t>();
    for (std::int32_t i = 0; i < nEntries; i++) {
        const auto foundNorm = readNormalizationMethod();
        const auto foundChrom = _fs->read<std::int32_t>();
        const auto foundUnit = readMatrixUnit();

        const auto foundResolution = _fs->read<std::int32_t>();
        const auto filePosition = _fs->read<std::int64_t>();
        const auto sizeInBytes = version() > 8
                                     ? _fs->read<std::int64_t>()
                                     : static_cast<std::int64_t>(_fs->read<std::int32_t>());
        if (foundChrom == chromId1 && foundNorm == wantedNorm && foundUnit == wantedUnit &&
            foundResolution == wantedResolution) {
            const auto currentPos = this->_fs->tellg();
            c1Norm = readNormalizationVector(indexEntry{filePosition, sizeInBytes});
            _fs->seekg(currentPos);
        }
        if (chromId1 != chromId2 && foundChrom == chromId2 && foundNorm == wantedNorm &&
            foundUnit == wantedUnit && foundResolution == wantedResolution) {
            const auto currentPos = this->_fs->tellg();
            c2Norm = readNormalizationVector(indexEntry{filePosition, sizeInBytes});
            _fs->seekg(currentPos);
        }
    }
    if (footer.c1Norm().empty() || footer.c2Norm().empty()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("File did not contain {} normalization vectors for one or both "
                                   "chromosomes at resolution {} and unit {}"),
                        wantedNorm, wantedResolution, wantedUnit));
    }

    return footer;
}
/*
inline void HiCFileStream::removeCachedFooter(const HiCFooter &footer) {
    auto it = this->_footers.find(footer);
    if (it == this->_footers.end()) {
        throw std::out_of_range("unable to find footer in cache");
    }
    this->_footers.erase(it);
}

inline void HiCFileStream::purgeFooterCache() { this->_footers.clear(); }
 */

}  // namespace internal
