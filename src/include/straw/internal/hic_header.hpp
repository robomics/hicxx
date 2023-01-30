// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "straw/internal/common.hpp"

namespace internal {

struct HiCHeader {
    std::string url{};
    std::int32_t version{-1};
    std::int64_t masterIndexOffset{-1};
    std::string genomeID{};
    std::int64_t nviPosition{-1};
    std::int64_t nviLength{-1};
    ChromosomeMap chromosomes{};
    std::vector<std::int32_t> resolutions{};

    constexpr explicit operator bool() const noexcept;
    bool operator==(const HiCHeader &other) const noexcept;
    bool operator!=(const HiCHeader &other) const noexcept;

    std::int32_t nChromosomes() const noexcept;
    std::int32_t nResolutions() const noexcept;

    const chromosome &getChromosome(std::int32_t id) const noexcept;
};

}  // namespace internal

template <>
struct std::hash<internal::HiCHeader> {
    inline std::size_t operator()(internal::HiCHeader const &h) const noexcept {
        return internal::hash_combine(0, h.url, h.masterIndexOffset);
    }
};

#include "../../../hic_header_impl.hpp"
