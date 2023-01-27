// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cassert>
#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

// pointer structure for reading blocks or matrices, holds the size and position
struct indexEntry {
    std::int64_t position{-1};
    std::int64_t size{-1};

    constexpr explicit operator bool() const noexcept { return size >= 0 && position >= 0; }
};

// sparse matrixType entry
struct contactRecord {
    std::int32_t binX{};
    std::int32_t binY{};
    float counts{};
};

// chromosome
struct chromosome {
    std::string name{};
    std::int32_t index{};
    std::int64_t length{};

    constexpr bool operator==(const chromosome &other) const noexcept {
        return index == other.index;
    }
    constexpr bool operator!=(const chromosome &other) const noexcept { return !(*this == other); }
    constexpr bool operator<(const chromosome &other) const noexcept { return index < other.index; }
    constexpr bool operator<=(const chromosome &other) const noexcept {
        return index <= other.index;
    }
    constexpr bool operator>(const chromosome &other) const noexcept { return index > other.index; }
    constexpr bool operator>=(const chromosome &other) const noexcept {
        return index >= other.index;
    }
};

using ChromosomeMap = std::map<std::string, chromosome>;

enum class Normalization { NONE, VC, VC_SQRT, KR, SCALE };
enum class MatrixType { observed, oe, expected };
enum class Unit { BP, FRAG };

inline Normalization ParseNormStr(const std::string &s) {
    if (s == "NONE") {
        return Normalization::NONE;
    }
    if (s == "VC") {
        return Normalization::VC;
    }

    if (s == "VC_SQRT") {
        return Normalization::VC_SQRT;
    }

    if (s == "KR") {
        return Normalization::KR;
    }

    if (s == "SCALE") {
        return Normalization::SCALE;
    }

    throw std::runtime_error("Invalid normalization \"" + s + "\"");
}

inline MatrixType ParseMatrixTypeStr(const std::string &s) {
    if (s == "observed") {
        return MatrixType::observed;
    }
    if (s == "oe") {
        return MatrixType::oe;
    }
    if (s == "expected") {
        return MatrixType::expected;
    }

    throw std::runtime_error("Invalid matrix type \"" + s + "\"");
}

inline Unit ParseUnitStr(const std::string &s) {
    if (s == "BP") {
        return Unit::BP;
    }
    if (s == "FRAG") {
        return Unit::FRAG;
    }

    throw std::runtime_error("Invalid unit \"" + s + "\"");
}

inline std::string to_string(Normalization n) {
    switch (n) {
        case Normalization::NONE:
            return "NONE";
        case Normalization::VC:
            return "VC";
        case Normalization::VC_SQRT:
            return "VC_SQRT";
        case Normalization::KR:
            return "KR";
        case Normalization::SCALE:
            return "SCALE";
    }
    assert(false);
    std::abort();
}

inline std::string to_string(MatrixType t) {
    switch (t) {
        case MatrixType::observed:
            return "observed";
        case MatrixType::oe:
            return "oe";
        case MatrixType::expected:
            return "expected";
    }
    assert(false);
    std::abort();
}

inline std::string to_string(Unit u) {
    switch (u) {
        case Unit::BP:
            return "BP";
        case Unit::FRAG:
            return "FRAG";
    }
    assert(false);
    std::abort();
}

namespace internal {
// Adapted from:
// https://www.boost.org/doc/libs/1_37_0/doc/html/hash/reference.html#boost.hash_combine

template <typename T>
inline std::size_t hash_combine(std::size_t seed, const T &v) {
    seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
    return seed;
}
template <typename T, typename... Args>
inline std::size_t hash_combine(std::size_t seed, const T &v, Args... args) {
    seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
    return hash_combine(seed, args...);
}

inline bool StartsWith(const std::string &s, const std::string &prefix) {
    return s.find(prefix) == 0;
}

}  // namespace internal

inline void convertGenomeToBinPos(const std::int64_t origRegionIndices[4],
                                  std::int64_t regionIndices[4], std::int32_t resolution) {
    for (std::uint16_t q = 0; q < 4; q++) {
        // used to find the blocks we need to access
        regionIndices[q] = origRegionIndices[q] / resolution;
    }
}
