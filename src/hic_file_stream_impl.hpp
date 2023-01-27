// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

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

constexpr HiCFooter::operator bool() const noexcept { return fileOffset >= 0; }

inline bool HiCFooter::operator==(const HiCFooter &other) const noexcept {
    return url == other.url && matrixType == other.matrixType &&
           normalization == other.normalization && unit == other.unit &&
           resolution == other.resolution;
}

inline bool HiCFooter::operator!=(const HiCFooter &other) const noexcept {
    return !(*this == other);
}

inline HiCFileStream::HiCFileStream(std::string url)
    : fs_(HiCFileStream::openStream(std::move(url))), header_(HiCFileStream::readHeader(fs_)) {}

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

inline const std::string &HiCFileStream::url() const noexcept { return fs_.url(); }
constexpr const HiCHeader &HiCFileStream::header() const noexcept { return header_; }

inline bool HiCFileStream::isLocal() const noexcept { return fs_.is_local(); }

inline bool HiCFileStream::isRemote() const noexcept { return fs_.is_remote(); }

inline std::int32_t HiCFileStream::version() const noexcept {
    assert(header_.version != -1);
    return header_.version;
}

inline void HiCFileStream::readCompressedBytes(indexEntry idx, std::string &buffer) {
    fs_.seekg(idx.position);
    fs_.read(buffer, idx.size);
}

inline void HiCFileStream::discardExpectedVector(std::int64_t nValues) {
    const auto elementSize = version() > 8 ? sizeof(float) : sizeof(double);
    fs_.seekg(nValues * elementSize, std::ios::cur);
}

inline std::vector<double> HiCFileStream::readExpectedVector(std::int64_t nValues) {
    std::vector<double> initialExpectedValues(nValues);
    if (version() > 8) {
        std::vector<float> tmpbuff(nValues);
        fs_.read(tmpbuff);
        std::transform(tmpbuff.begin(), tmpbuff.end(), initialExpectedValues.begin(),
                       [](float n) { return static_cast<double>(n); });
    } else {
        fs_.read(initialExpectedValues);
    }

    // This seems to be copying initialValues into finalResult at the moment
    // std::int32_t window = 5000000 / resolution;
    // rollingMedian(initialExpectedValues, expectedValues, window);
    return initialExpectedValues;
}

inline std::vector<double> HiCFileStream::readNormalizationFactors(std::int32_t wantedChrom) {
    const auto nFactors = fs_.read<std::int32_t>();
    std::vector<double> normFactors{};
    auto readFactor = [this]() -> double {
        if (version() > 8) {
            return fs_.read<float>();
        }
        return fs_.read<double>();
    };

    for (auto i = 0; i < nFactors; ++i) {
        const auto foundChrom = fs_.read<std::int32_t>();
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

inline void HiCFileStream::discardNormalizationFactors(std::int32_t wantedChrom) {
    (void)readNormalizationFactors(wantedChrom);
}

inline MatrixType HiCFileStream::readMatrixType() {
    std::string strbuff;
    fs_.getline(strbuff, '\0');
    return ParseMatrixTypeStr(strbuff);
}

inline Normalization HiCFileStream::readNormalization() {
    std::string strbuff;
    fs_.getline(strbuff, '\0');
    return ParseNormStr(strbuff);
}

inline Unit HiCFileStream::readUnit() {
    std::string strbuff;
    fs_.getline(strbuff, '\0');
    return ParseUnitStr(strbuff);
}

inline std::int64_t HiCFileStream::readNValues() {
    if (version() > 8) {
        return fs_.read<std::int64_t>();
    }
    return fs_.read<std::int32_t>();
}

inline bool HiCFileStream::checkMagicString(filestream::FileStream &fs) {
    return fs.getline('\0') == "HIC";
}

inline std::int64_t HiCFileStream::masterOffset() const noexcept {
    return header_.masterIndexOffset;
}

inline bool HiCFileStream::checkMagicString() { return checkMagicString(fs_); }

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

inline HiCFooter HiCFileStream::readFooter(const chromosome &chrom1, const chromosome &chrom2,
                                           MatrixType matrixType, Normalization norm, Unit unit,
                                           std::int32_t wantedResolution) {
    return readFooter(chrom1.index, chrom2.index, matrixType, norm, unit, wantedResolution);
}
inline HiCFooter HiCFileStream::readFooter(const std::string &chrom1, const std::string &chrom2,
                                           MatrixType matrixType, Normalization norm, Unit unit,
                                           std::int32_t wantedResolution) {
    auto it1 = this->header_.chromosomes.find(chrom1);
    if (it1 == this->header_.chromosomes.end()) {
        throw std::out_of_range("no such chromosome \"" + chrom1 + "\"");
    }
    auto it2 = this->header_.chromosomes.find(chrom2);
    if (it2 == this->header_.chromosomes.end()) {
        throw std::out_of_range("no such chromosome \"" + chrom2 + "\"");
    }

    return readFooter(it1->second, it2->second, matrixType, norm, unit, wantedResolution);
}

inline HiCFooter HiCFileStream::readFooter(const chromosome &chrom, MatrixType matrixType,
                                           Normalization norm, Unit unit,
                                           std::int32_t wantedResolution) {
    return readFooter(chrom, chrom, matrixType, norm, unit, wantedResolution);
}
inline HiCFooter HiCFileStream::readFooter(const std::string &chrom, MatrixType matrixType,
                                           Normalization norm, Unit unit,
                                           std::int32_t wantedResolution) {
    return readFooter(chrom, chrom, matrixType, norm, unit, wantedResolution);
}
inline HiCFooter HiCFileStream::readFooter(std::int32_t chromId, MatrixType matrixType,
                                           Normalization norm, Unit unit,
                                           std::int32_t wantedResolution) {
    return readFooter(chromId, chromId, matrixType, norm, unit, wantedResolution);
}

inline HiCFooter HiCFileStream::readFooter(std::int32_t chromId1, std::int32_t chromId2,
                                           const MatrixType wantedMatrixType,
                                           const Normalization wantedNorm, const Unit wantedUnit,
                                           std::int32_t wantedResolution) {
    using MT = MatrixType;
    using NM = Normalization;

    HiCFooter footer{fs_.url(), wantedMatrixType, wantedNorm, wantedUnit, wantedResolution};

    auto it = footers_.find(footer);
    if (it != footers_.end()) {
        return *it;
    }

    fs_.seekg(masterOffset());

    auto insert_footer_and_return = [&]() {
        const auto node = footers_.emplace(std::move(footer));
        assert(node.second);
        return *node.first;
    };

    if (chromId1 > chromId2) {
        std::swap(chromId1, chromId2);
    }

    const auto key = std::to_string(chromId1) + "_" + std::to_string(chromId2);

    (void)readNValues();  // nBytes

    auto nEntries = fs_.read<std::int32_t>();
    for (int i = 0; i < nEntries; i++) {
        const auto strbuff = fs_.getline('\0');
        const auto fpos = fs_.read<std::int64_t>();
        (void)fs_.read<std::int32_t>();  // sizeInBytes
        if (strbuff == key) {
            footer.fileOffset = fpos;
        }
    }
    if (footer.fileOffset == -1) {
        throw std::runtime_error("File doesn't have the given chr_chr map " + key);
    }

    if ((wantedMatrixType == MT::observed && wantedNorm == NM::NONE) ||
        ((wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) &&
         wantedNorm == NM::NONE && chromId1 != chromId2)) {
        return insert_footer_and_return();  // no need to read wantedNorm vector index
    }

    // read in and ignore expected value maps; don't store; reading these to
    // get to wantedNorm vector index
    auto nExpectedValues = fs_.read<std::int32_t>();
    for (std::int32_t i = 0; i < nExpectedValues; i++) {
        const auto foundUnit = readUnit();
        const auto foundResolution = fs_.read<std::int32_t>();
        const auto nValues = readNValues();

        bool store = chromId1 == chromId2 &&
                     (wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) &&
                     wantedNorm == NM::NONE && foundUnit == wantedUnit &&
                     foundResolution == wantedResolution;

        if (store) {
            footer.expectedValues = readExpectedVector(nValues);
            const auto normFactors = readNormalizationFactors(chromId1);
            applyNormalizationFactors(footer.expectedValues, normFactors);

        } else {
            discardExpectedVector(nValues);
            discardNormalizationFactors(chromId1);
        }
    }

    if (chromId1 == chromId2 && (wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) &&
        wantedNorm == NM::NONE) {
        if (footer.expectedValues.empty()) {
            throw std::runtime_error("File did not contain expected values vectors at " +
                                     std::to_string(wantedResolution) + " " +
                                     to_string(wantedUnit));
        }
        return insert_footer_and_return();
    }

    nExpectedValues = fs_.read<std::int32_t>();
    for (std::int32_t i = 0; i < nExpectedValues; i++) {
        const auto foundNorm = readNormalization();
        const auto foundUnit = readUnit();
        const auto foundResolution = fs_.read<std::int32_t>();

        const auto nValues = readNValues();
        bool store = chromId1 == chromId2 &&
                     (wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) &&
                     foundNorm == wantedNorm && foundUnit == wantedUnit &&
                     foundResolution == wantedResolution;

        if (store) {
            footer.expectedValues = readExpectedVector(nValues);
            const auto normFactors = readNormalizationFactors(chromId1);
            applyNormalizationFactors(footer.expectedValues, normFactors);
        } else {
            discardExpectedVector(nValues);
            discardNormalizationFactors(chromId1);
        }
    }

    if (chromId1 == chromId2 && (wantedMatrixType == MT::oe || wantedMatrixType == MT::expected) &&
        wantedNorm != NM::NONE) {
        if (footer.expectedValues.empty()) {
            throw std::runtime_error("File did not contain normalized expected values vectors at " +
                                     std::to_string(wantedResolution) + " " +
                                     to_string(wantedUnit));
        }
    }

    // Index of normalization vectors
    nEntries = fs_.read<std::int32_t>();
    for (std::int32_t i = 0; i < nEntries; i++) {
        const auto foundNorm = readNormalization();
        const auto foundChrom = fs_.read<std::int32_t>();
        const auto foundUnit = readUnit();

        const auto foundResolution = fs_.read<std::int32_t>();
        const auto filePosition = fs_.read<std::int64_t>();
        const auto sizeInBytes = version() > 8
                                     ? fs_.read<std::int64_t>()
                                     : static_cast<std::int64_t>(fs_.read<std::int32_t>());
        if (foundChrom == chromId1 && foundNorm == wantedNorm && foundUnit == wantedUnit &&
            foundResolution == wantedResolution) {
            footer.c1NormEntry = indexEntry{filePosition, sizeInBytes};
        }
        if (foundChrom == chromId2 && foundNorm == wantedNorm && foundUnit == wantedUnit &&
            foundResolution == wantedResolution) {
            footer.c2NormEntry = indexEntry{filePosition, sizeInBytes};
        }
    }
    if (!footer.c1NormEntry || !footer.c2NormEntry) {
        throw std::runtime_error("File did not contain " + to_string(wantedNorm) +
                                 " normalization vectors for one or both chromosomes at " +
                                 std::to_string(wantedResolution) + " " + to_string(wantedUnit));
    }
    return insert_footer_and_return();
}

inline void HiCFileStream::removeCachedFooter(const HiCFooter &footer) {
    auto it = this->footers_.find(footer);
    if (it == this->footers_.end()) {
        throw std::out_of_range("unable to find footer in cache");
    }
    this->footers_.erase(it);
}

inline void HiCFileStream::purgeFooterCache() { this->footers_.clear(); }

}  // namespace internal
