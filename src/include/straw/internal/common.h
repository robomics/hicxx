// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
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

template <>
struct std::hash<chromosome> {
    inline std::size_t operator()(chromosome const &c) const noexcept {
        return std::hash<std::int32_t>{}(c.index);
    }
};

using ChromosomeMap = std::map<std::string, chromosome>;

enum class NormalizationMethod { NONE, VC, VC_SQRT, KR, SCALE };
enum class MatrixType { observed, oe, expected };
enum class MatrixUnit { BP, FRAG };

inline NormalizationMethod ParseNormStr(const std::string &s) {
    if (s == "NONE") {
        return NormalizationMethod::NONE;
    }
    if (s == "VC") {
        return NormalizationMethod::VC;
    }

    if (s == "VC_SQRT") {
        return NormalizationMethod::VC_SQRT;
    }

    if (s == "KR") {
        return NormalizationMethod::KR;
    }

    if (s == "SCALE") {
        return NormalizationMethod::SCALE;
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

inline MatrixUnit ParseUnitStr(const std::string &s) {
    if (s == "BP") {
        return MatrixUnit::BP;
    }
    if (s == "FRAG") {
        return MatrixUnit::FRAG;
    }

    throw std::runtime_error("Invalid unit \"" + s + "\"");
}
/*
inline std::string to_string(NormalizationMethod n) {
    switch (n) {
        case NormalizationMethod::NONE:
            return "NONE";
        case NormalizationMethod::VC:
            return "VC";
        case NormalizationMethod::VC_SQRT:
            return "VC_SQRT";
        case NormalizationMethod::KR:
            return "KR";
        case NormalizationMethod::SCALE:
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

inline std::string to_string(MatrixUnit u) {
    switch (u) {
        case MatrixUnit::BP:
            return "BP";
        case MatrixUnit::FRAG:
            return "FRAG";
    }
    assert(false);
    std::abort();
}
 */

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

template <>
struct fmt::formatter<NormalizationMethod> {
    static constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
        if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
            throw fmt::format_error("invalid format");
        }
        return ctx.end();
    }

    template <class FormatContext>
    static auto format(const NormalizationMethod n, FormatContext &ctx) -> decltype(ctx.out()) {
        switch (n) {
            case NormalizationMethod::NONE:
                return fmt::format_to(ctx.out(), FMT_STRING("NONE"));
            case NormalizationMethod::VC:
                return fmt::format_to(ctx.out(), FMT_STRING("VC"));
            case NormalizationMethod::VC_SQRT:
                return fmt::format_to(ctx.out(), FMT_STRING("VC_SQRT"));
            case NormalizationMethod::KR:
                return fmt::format_to(ctx.out(), FMT_STRING("KR"));
            case NormalizationMethod::SCALE:
                return fmt::format_to(ctx.out(), FMT_STRING("SCALE"));
        }
        assert(false);
        std::abort();
    }
};

template <>
struct fmt::formatter<MatrixType> {
    static constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
        if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
            throw fmt::format_error("invalid format");
        }
        return ctx.end();
    }

    template <class FormatContext>
    static auto format(const MatrixType t, FormatContext &ctx) -> decltype(ctx.out()) {
        switch (t) {
            case MatrixType::observed:
                return fmt::format_to(ctx.out(), FMT_STRING("observed"));
            case MatrixType::oe:
                return fmt::format_to(ctx.out(), FMT_STRING("oe"));
            case MatrixType::expected:
                return fmt::format_to(ctx.out(), FMT_STRING("expected"));
        }
        assert(false);
        std::abort();
    }
};

template <>
struct fmt::formatter<MatrixUnit> {
    static constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
        if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
            throw fmt::format_error("invalid format");
        }
        return ctx.end();
    }

    template <class FormatContext>
    static auto format(const MatrixUnit u, FormatContext &ctx) -> decltype(ctx.out()) {
        switch (u) {
            case MatrixUnit::BP:
                return fmt::format_to(ctx.out(), FMT_STRING("BP"));
            case MatrixUnit::FRAG:
                return fmt::format_to(ctx.out(), FMT_STRING("FRAG"));
        }
        assert(false);
        std::abort();
    }
};

template <typename T>
using UniquePtrWithDeleter = std::unique_ptr<T, std::function<void(T *)>>;
