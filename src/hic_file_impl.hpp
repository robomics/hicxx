// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "straw/internal/common.h"

inline HiCFile::HiCFile(std::string url_)
    : _fs(std::make_shared<internal::HiCFileStream>(std::move(url_))) {}

inline const std::string& HiCFile::url() const noexcept { return _fs->url(); }

inline const std::string& HiCFile::name() const noexcept { return url(); }

inline std::int32_t HiCFile::version() const noexcept { return _fs->version(); }

inline const ChromosomeMap& HiCFile::chromosomes() const noexcept {
    return _fs->header().chromosomes;
}

inline const std::string& HiCFile::genomeID() const noexcept { return _fs->header().genomeID; }

inline const std::vector<std::int32_t>& HiCFile::resolutions() const noexcept {
    return _fs->header().resolutions;
}

inline std::shared_ptr<const internal::HiCFooter> HiCFile::getFooter(
    std::int32_t chromId1, std::int32_t chromId2, MatrixType matrixType, NormalizationMethod norm,
    MatrixUnit unit, std::int32_t resolution) {
    const internal::HiCFooterMetadata metadata{url(), matrixType, norm, unit, resolution};
    auto it = _footers.find(metadata);
    if (it != _footers.end()) {
        return it->second;
    }

    auto footer = std::make_shared<const internal::HiCFooter>(
        _fs->readFooter(chromId1, chromId2, matrixType, norm, unit, resolution));
    auto node = _footers.emplace(std::move(metadata), std::move(footer));

    assert(node.second);
    return node.first->second;
}

inline internal::MatrixZoomData HiCFile::getMatrixZoomData(const chromosome& chrom,
                                                           MatrixType matrixType,
                                                           NormalizationMethod norm,
                                                           MatrixUnit unit,
                                                           std::int32_t resolution) {
    return getMatrixZoomData(chrom, chrom, matrixType, norm, unit, resolution);
}
inline internal::MatrixZoomData HiCFile::getMatrixZoomData(const std::string& chromName,
                                                           MatrixType matrixType,
                                                           NormalizationMethod norm,
                                                           MatrixUnit unit,
                                                           std::int32_t resolution) {
    return getMatrixZoomData(chromName, chromName, matrixType, norm, unit, resolution);
}
inline internal::MatrixZoomData HiCFile::getMatrixZoomData(std::int32_t chromId,
                                                           MatrixType matrixType,
                                                           NormalizationMethod norm,
                                                           MatrixUnit unit,
                                                           std::int32_t resolution) {
    return getMatrixZoomData(chromId, chromId, matrixType, norm, unit, resolution);
}

inline internal::MatrixZoomData HiCFile::getMatrixZoomData(
    const chromosome& chrom1, const chromosome& chrom2, MatrixType matrixType,
    NormalizationMethod norm, MatrixUnit unit, std::int32_t resolution) {
    return getMatrixZoomData(chrom1.index, chrom2.index, matrixType, norm, unit, resolution);
}

inline internal::MatrixZoomData HiCFile::getMatrixZoomData(
    const std::string& chromName1, const std::string& chromName2, MatrixType matrixType,
    NormalizationMethod norm, MatrixUnit unit, std::int32_t resolution) {
    const auto it1 = chromosomes().find(chromName1);
    if (it1 == chromosomes().end()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("unable to find chromosome named {}"), chromName1));
    }
    if (chromName1 == chromName2) {
        return getMatrixZoomData(it1->second, it1->second, matrixType, norm, unit, resolution);
    }

    const auto it2 = chromosomes().find(chromName2);
    if (it2 == chromosomes().end()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("unable to find chromosome named {}"), chromName2));
    }

    return getMatrixZoomData(it1->second, it2->second, matrixType, norm, unit, resolution);
}

inline internal::MatrixZoomData HiCFile::getMatrixZoomData(
    std::int32_t chromId1, std::int32_t chromId2, MatrixType matrixType, NormalizationMethod norm,
    MatrixUnit unit, std::int32_t resolution) {
    if (chromId1 >= chromosomes().size()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("unable to find chromosome corresponding to ID {}"), chromId1));
    }
    if (chromId2 >= chromosomes().size()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("unable to find chromosome corresponding to ID {}"), chromId2));
    }

    if (chromId1 > chromId2) {
        std::swap(chromId1, chromId2);
    }

    if (matrixType == MatrixType::expected && norm != NormalizationMethod::NONE) {
        throw std::logic_error(
            fmt::format(FMT_STRING("matrix type {} is incompatible with normalization method {}"),
                        matrixType, norm));
    }

    const auto it = std::find(resolutions().begin(), resolutions().end(), resolution);
    if (it == resolutions().end()) {
        throw std::runtime_error(fmt::format(
            FMT_STRING(
                "matrix does not have interactions for resolution {}. Available resolutions: {}"),
            resolution, fmt::join(_fs->header().resolutions, ", ")));
    }

    return internal::MatrixZoomData(
        _fs, getFooter(chromId1, chromId2, matrixType, norm, unit, resolution));
}

inline GenomicCoordinates GenomicCoordinates::fromString(std::string coord, bool noChromName) {
    GenomicCoordinates gc{};

    const auto original_coords = coord;

    if (!noChromName) {
        auto pos = coord.find(':');
        if (pos == std::string::npos) {
            gc.chrom = coord;
            return gc;
        }

        gc.chrom = coord.substr(0, pos);
        coord = coord.substr(pos + 1);
    }

    auto pos = coord.find('-');
    if (pos == std::string::npos) {
        pos = coord.find(':');
    }
    if (pos == std::string::npos) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("unable to parse coordinate \"{}\""), coord));
    }

    try {
        gc.start = std::stoi(coord.substr(0, pos));
        gc.end = std::stoi(coord.substr(pos + 1));
        if (gc.start >= gc.end) {
            throw std::runtime_error(fmt::format(
                FMT_STRING("invalid coordinate {}: start position >= end position"), coord));
        }
        if (gc.start < 0) {
            throw std::runtime_error(fmt::format(
                FMT_STRING("invalid coordinate {}: start position is negative"), coord));
        }
    } catch (const std::exception& e) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("unable to parse coordinate \"{}\": {}"), coord, e.what()));
    }
    return gc;
}

inline std::size_t HiCFile::numCachedFooters() const noexcept { return _footers.size(); }

inline void HiCFile::purgeFooterCache() { _footers.clear(); }
