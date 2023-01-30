// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "straw/internal/common.hpp"

namespace internal {
struct HiCFooterMetadata {
    std::string url{};
    MatrixType matrixType{MatrixType::observed};
    NormalizationMethod normalization{NormalizationMethod::NONE};
    MatrixUnit unit{MatrixUnit::BP};
    std::int32_t resolution{-1};
    chromosome chrom1{};
    chromosome chrom2{};
    std::int64_t fileOffset{-1};

    constexpr explicit operator bool() const noexcept;
    bool operator==(const HiCFooterMetadata &other) const noexcept;
    bool operator!=(const HiCFooterMetadata &other) const noexcept;
};

class HiCFooter {
    HiCFooterMetadata _metadata{};
    std::vector<double> _expectedValues{};
    std::vector<double> _c1Norm{};
    std::vector<double> _c2Norm{};

   public:
    HiCFooter() = default;
    explicit HiCFooter(HiCFooterMetadata metadata_) noexcept;

    constexpr explicit operator bool() const noexcept;
    bool operator==(const HiCFooter &other) const noexcept;
    bool operator!=(const HiCFooter &other) const noexcept;

    constexpr const HiCFooterMetadata &metadata() const noexcept;
    constexpr HiCFooterMetadata &metadata() noexcept;

    constexpr const std::string &url() const noexcept;
    constexpr MatrixType matrixType() const noexcept;
    constexpr NormalizationMethod normalization() const noexcept;
    constexpr MatrixUnit unit() const noexcept;
    constexpr std::int32_t resolution() const noexcept;
    constexpr const chromosome &chrom1() const noexcept;
    constexpr const chromosome &chrom2() const noexcept;
    constexpr std::int64_t fileOffset() const noexcept;

    constexpr const std::vector<double> &expectedValues() const noexcept;
    constexpr const std::vector<double> &c1Norm() const noexcept;
    constexpr const std::vector<double> &c2Norm() const noexcept;

    constexpr std::vector<double> &expectedValues() noexcept;
    constexpr std::vector<double> &c1Norm() noexcept;
    constexpr std::vector<double> &c2Norm() noexcept;
};
}  // namespace internal

#include "../../../hic_footer_impl.hpp"
