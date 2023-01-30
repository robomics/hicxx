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

namespace hicxx {

// pointer structure for reading blocks or matrices, holds the size and position
struct indexEntry {
    std::int64_t position{-1};
    std::int64_t size{-1};

    constexpr explicit operator bool() const noexcept { return size >= 0 && position >= 0; }
    constexpr bool operator<(const indexEntry &other) const noexcept {
        return position < other.position;
    }
    constexpr bool operator==(const indexEntry &other) const noexcept {
        return position == other.position && size == other.size;
    }
    constexpr bool operator!=(const indexEntry &other) const noexcept { return !(*this == other); }
};

// sparse matrixType entry
struct contactRecord {
    std::int32_t bin1_start{};
    std::int32_t bin2_start{};
    float count{};

    constexpr bool operator<(const contactRecord &other) const noexcept {
        if (bin2_start == other.bin2_start) {
            return bin1_start < other.bin1_start;
        }
        return bin2_start < other.bin2_start;
    }
    constexpr bool operator==(const contactRecord &other) const noexcept {
        return bin1_start == other.bin1_start && bin2_start == other.bin2_start &&
               count == other.count;
    }
    constexpr bool operator!=(const contactRecord &other) const noexcept {
        return !(*this == other);
    }
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
}  // namespace hicxx

template <>
struct std::hash<hicxx::chromosome> {
    inline std::size_t operator()(hicxx::chromosome const &c) const noexcept {
        return std::hash<std::int32_t>{}(c.index);
    }
};

namespace hicxx {

using ChromosomeMap = std::map<std::string, chromosome>;

enum class NormalizationMethod { NONE, VC, VC_SQRT, KR, SCALE };
enum class MatrixType { observed, oe, expected };
enum class MatrixUnit { BP, FRAG };

[[nodiscard]] inline NormalizationMethod ParseNormStr(const std::string &s) {
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

[[nodiscard]] inline MatrixType ParseMatrixTypeStr(const std::string &s) {
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

[[nodiscard]] inline MatrixUnit ParseUnitStr(const std::string &s) {
    if (s == "BP") {
        return MatrixUnit::BP;
    }
    if (s == "FRAG") {
        return MatrixUnit::FRAG;
    }

    throw std::runtime_error("Invalid unit \"" + s + "\"");
}

namespace internal {
// Adapted from:
// https://www.boost.org/doc/libs/1_37_0/doc/html/hash/reference.html#boost.hash_combine

template <typename T>
[[nodiscard]] inline std::size_t hash_combine(std::size_t seed, const T &v) {
    seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
    return seed;
}
template <typename T, typename... Args>
[[nodiscard]] inline std::size_t hash_combine(std::size_t seed, const T &v, Args... args) {
    seed ^= std::hash<T>{}(v) + 0x9e3779b9 + (seed << 6U) + (seed >> 2U);
    return hash_combine(seed, args...);
}

inline bool StartsWith(const std::string &s, const std::string &prefix) {
    return s.find(prefix) == 0;
}

}  // namespace internal

// to avoid useless casts (see https://github.com/nlohmann/json/issues/2893#issuecomment-889152324)
template <class T, class U>
[[maybe_unused]] [[nodiscard]] constexpr T conditional_static_cast(U value) {
    if constexpr (std::is_same_v<T, U>) {
        return value;
    } else {
        return static_cast<T>(value);
    }
}
}  // namespace hicxx

template <>
struct fmt::formatter<hicxx::NormalizationMethod> {
    static constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
        if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
            throw fmt::format_error("invalid format");
        }
        return ctx.end();
    }

    template <class FormatContext>
    static auto format(const hicxx::NormalizationMethod n, FormatContext &ctx)
        -> decltype(ctx.out()) {
        switch (n) {
            case hicxx::NormalizationMethod::NONE:
                return fmt::format_to(ctx.out(), FMT_STRING("NONE"));
            case hicxx::NormalizationMethod::VC:
                return fmt::format_to(ctx.out(), FMT_STRING("VC"));
            case hicxx::NormalizationMethod::VC_SQRT:
                return fmt::format_to(ctx.out(), FMT_STRING("VC_SQRT"));
            case hicxx::NormalizationMethod::KR:
                return fmt::format_to(ctx.out(), FMT_STRING("KR"));
            case hicxx::NormalizationMethod::SCALE:
                return fmt::format_to(ctx.out(), FMT_STRING("SCALE"));
        }
        assert(false);
        std::abort();
    }
};

template <>
struct fmt::formatter<hicxx::MatrixType> {
    static constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
        if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
            throw fmt::format_error("invalid format");
        }
        return ctx.end();
    }

    template <class FormatContext>
    static auto format(const hicxx::MatrixType t, FormatContext &ctx) -> decltype(ctx.out()) {
        switch (t) {
            case hicxx::MatrixType::observed:
                return fmt::format_to(ctx.out(), FMT_STRING("observed"));
            case hicxx::MatrixType::oe:
                return fmt::format_to(ctx.out(), FMT_STRING("oe"));
            case hicxx::MatrixType::expected:
                return fmt::format_to(ctx.out(), FMT_STRING("expected"));
        }
        assert(false);
        std::abort();
    }
};

template <>
struct fmt::formatter<hicxx::MatrixUnit> {
    static constexpr auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
        if (ctx.begin() != ctx.end() && *ctx.begin() != '}') {
            throw fmt::format_error("invalid format");
        }
        return ctx.end();
    }

    template <class FormatContext>
    static auto format(const hicxx::MatrixUnit u, FormatContext &ctx) -> decltype(ctx.out()) {
        switch (u) {
            case hicxx::MatrixUnit::BP:
                return fmt::format_to(ctx.out(), FMT_STRING("BP"));
            case hicxx::MatrixUnit::FRAG:
                return fmt::format_to(ctx.out(), FMT_STRING("FRAG"));
        }
        assert(false);
        std::abort();
    }
};

template <typename T>
using UniquePtrWithDeleter = std::unique_ptr<T, std::function<void(T *)>>;

struct GenomicCoordinates {
    std::string chrom;
    std::int32_t start;
    std::int32_t end;

    [[nodiscard]] inline static GenomicCoordinates fromString(std::string coord,
                                                              bool noChromName = false) {
        GenomicCoordinates gc{};

        const auto original_coord = coord;

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
            std::size_t tail{0};
            gc.start = std::stoi(coord.substr(0, pos));
            gc.end = std::stoi(coord.substr(pos + 1), &tail);
            if (gc.start >= gc.end) {
                throw std::runtime_error(fmt::format(
                    FMT_STRING("invalid coordinate {}: start position >= end position"), coord));
            }
            if (gc.start < 0) {
                throw std::runtime_error(fmt::format(
                    FMT_STRING("invalid coordinate {}: start position is negative"), coord));
            }
            coord = coord.substr(pos + 1);
            if (tail != coord.size()) {
                throw std::runtime_error(fmt::format(FMT_STRING("unable to parse \"{}\""), coord));
            }
        } catch (const std::exception &e) {
            throw std::runtime_error(fmt::format(
                FMT_STRING("unable to parse coordinate \"{}\": {}"), original_coord, e.what()));
        }
        return gc;
    }
};
