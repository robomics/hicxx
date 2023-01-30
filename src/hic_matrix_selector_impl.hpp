// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ios>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "straw/internal/common.hpp"

namespace hicxx::internal {

template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type *>
inline T MatrixZoomData::BinaryBuffer::read() {
    static_assert(sizeof(char) == 1, "");
    assert(i < buffer.size());
    T x{};

    std::memcpy(static_cast<void *>(&x), buffer.data() + i, sizeof(T));
    i += sizeof(T);
    return x;
}

inline MatrixZoomData::MatrixZoomData(std::shared_ptr<HiCFileStream> fs,
                                      std::shared_ptr<const HiCFooter> footer)
    : _fs(std::move(fs)), _footer(std::move(footer)), _blockMap(readBlockMap(*_fs, *_footer)) {}

inline const chromosome &MatrixZoomData::chrom1() const noexcept { return _footer->chrom1(); }

inline const chromosome &MatrixZoomData::chrom2() const noexcept { return _footer->chrom2(); }

inline std::int64_t MatrixZoomData::resolution() const noexcept { return _footer->resolution(); }

inline MatrixType MatrixZoomData::matrixType() const noexcept { return _footer->matrixType(); }

inline NormalizationMethod MatrixZoomData::normalizationMethod() const noexcept {
    return _footer->normalization();
}

inline MatrixUnit MatrixZoomData::matrixUnit() const noexcept { return _footer->unit(); }

inline std::int64_t MatrixZoomData::numBins1() const noexcept {
    return (chrom1().length + resolution() - 1) / resolution();
}

inline std::int64_t MatrixZoomData::numBins2() const noexcept {
    return (chrom2().length + resolution() - 1) / resolution();
}

inline bool MatrixZoomData::isIntra() const noexcept { return chrom1() == chrom2(); }

inline bool MatrixZoomData::isInter() const noexcept { return !isIntra(); }

inline const std::vector<double> &MatrixZoomData::chrom1Norm() const noexcept {
    return _footer->c1Norm();
}

inline const std::vector<double> &MatrixZoomData::chrom2Norm() const noexcept {
    return _footer->c2Norm();
}

inline double MatrixZoomData::avgCount() const {
    if (isInter()) {
        return _blockMap.sumCount / static_cast<double>(numBins1() * numBins2());
    }
    throw std::domain_error(
        "MatrixZoomData::avgCount is not implemented for intra-chromosomal matrices");
}

inline void MatrixZoomData::fetch(std::vector<contactRecord> &buffer) {
    return fetch(0, chrom1().length, 0, chrom2().length, buffer);
}

inline void MatrixZoomData::fetch(const std::string &coord, std::vector<contactRecord> &buffer) {
    return fetch(coord, coord, buffer);
}

inline void MatrixZoomData::fetch(const std::string &coord1, const std::string &coord2,
                                  std::vector<contactRecord> &buffer) {
    auto coord1_ = GenomicCoordinates::fromString(coord1, true);
    auto coord2_ = GenomicCoordinates::fromString(coord2, true);

    return fetch(coord1_.start, coord1_.end, coord2_.start, coord2_.end, buffer);
}

inline void MatrixZoomData::fetch(std::int64_t start, std::int64_t end,
                                  std::vector<contactRecord> &buffer) {
    return fetch(start, end, start, end, buffer);
}

inline void MatrixZoomData::fetch(std::int64_t start1, std::int64_t end1, std::int64_t start2,
                                  std::int64_t end2, std::vector<contactRecord> &buffer) {
    if (start1 > end1) {
        throw std::invalid_argument(
            fmt::format(FMT_STRING("start1 > end1: {} > {}"), start1, end1));
    }
    if (start2 > end2) {
        throw std::invalid_argument(
            fmt::format(FMT_STRING("start2 > end2: {} > {}"), start2, end2));
    }

    if (start1 < 0 || end1 > chrom1().length) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("query extends past chromosome {}: interval {}-{} lies outside of 0-{}"),
            chrom1().name, start1, end1, chrom1().length));
    }

    if (start2 < 0 || end2 > chrom1().length) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("query extends past chromosome {}: interval {}-{} lies outside of 0-{}"),
            chrom2().name, start2, end2, chrom2().length));
    }

    const auto bin1 = start1 / resolution();
    const auto bin2 = end1 / resolution();
    const auto bin3 = start2 / resolution();
    const auto bin4 = end2 / resolution();

    if (_fs->version() > 8 && isIntra()) {
        readBlockNumbersV9Intra(bin1, bin2, bin3, bin4, _blockNumberBuff);
    } else {
        readBlockNumbers(bin1, bin2, bin3, bin4, _blockNumberBuff);
    }

    buffer.clear();
    for (auto blockNumber : _blockNumberBuff) {
        readBlockOfInteractions(_blockMap.blocks[blockNumber], _contactRecordBuff);

        for (auto &&record : _contactRecordBuff) {
            record.bin1_start *= std::int32_t(resolution());
            record.bin2_start *= std::int32_t(resolution());

            const auto pos1 = record.bin1_start;
            const auto pos2 = record.bin2_start;

            // Obs we use open-closed interval instead of open-open like is done in straw
            const auto overlapsInter =
                pos1 >= start1 && pos1 < end1 && pos2 >= start2 && pos2 < end2;
            const auto overlapsIntra =
                isIntra() && (pos2 >= start1 && pos2 < end1 && pos1 >= start2 && pos1 < end2);

            if (overlapsInter || overlapsIntra) {
                processInteraction(record);
                if (!std::isnan(record.count) && !std::isinf(record.count)) {
                    buffer.emplace_back(std::move(record));
                }
            }
        }
    }
}

inline void MatrixZoomData::fetch(std::vector<std::vector<float>> &buffer) {
    return fetch(0, chrom1().length, 0, chrom2().length, buffer);
}

inline void MatrixZoomData::fetch(const std::string &coord,
                                  std::vector<std::vector<float>> &buffer) {
    return fetch(coord, coord, buffer);
}

inline void MatrixZoomData::fetch(std::int64_t start, std::int64_t end,
                                  std::vector<std::vector<float>> &buffer) {
    return fetch(start, end, start, end, buffer);
}

inline void MatrixZoomData::fetch(const std::string &coord1, const std::string &coord2,
                                  std::vector<std::vector<float>> &buffer) {
    const auto coord1_ = GenomicCoordinates::fromString(coord1, true);
    const auto coord2_ = GenomicCoordinates::fromString(coord2, true);

    return fetch(coord1_.start, coord1_.end, coord2_.start, coord2_.end, buffer);
}

inline void MatrixZoomData::fetch(std::int64_t start1, std::int64_t end1, std::int64_t start2,
                                  std::int64_t end2, std::vector<std::vector<float>> &buffer) {
    const auto records = [&]() {
        std::vector<contactRecord> records_;
        fetch(start1, end1, start2, end2, records_);
        return records_;
    }();

    // We resize the buffer here so that we let fetch() deal with invalid queries
    const auto nRows = static_cast<std::size_t>((end1 - start1 + resolution() - 1) / resolution());
    const auto nCols = static_cast<std::size_t>((end2 - start2 + resolution() - 1) / resolution());

    buffer.resize(nRows);
    for (auto &row : buffer) {
        row.clear();
        row.resize(nCols, 0.0F);
    }

    if (records.empty()) {
        return;
    }

    const auto rowOffset = static_cast<std::size_t>(start1 / resolution());
    const auto colOffset = static_cast<std::size_t>(start2 / resolution());
    for (const auto &record : records) {
        assert(std::size_t(record.bin2_start) >= rowOffset);
        assert(std::size_t(record.bin1_start) >= colOffset);
        const auto i = std::size_t(record.bin2_start) - rowOffset;
        const auto j = std::size_t(record.bin1_start) - colOffset;

        buffer[i][j] = record.count;
    }
}

inline void MatrixZoomData::processInteraction(contactRecord &record) {
    const auto &c1Norm = _footer->c1Norm();
    const auto &c2Norm = _footer->c2Norm();
    const auto &expected = _footer->expectedValues();

    const auto skipNormalization =
        normalizationMethod() == NormalizationMethod::NONE || matrixType() == MatrixType::expected;

    if (!skipNormalization) {
        const auto bin1 = static_cast<std::size_t>(record.bin1_start / resolution());
        const auto bin2 = static_cast<std::size_t>(record.bin2_start / resolution());
        record.count /= static_cast<float>(c1Norm[bin1] * c2Norm[bin2]);
    }

    if (matrixType() == MatrixType::observed) {
        return;
    }

    const auto expectedCount = [&]() {
        if (isInter()) {
            return float(avgCount());
        }

        const auto i = static_cast<std::size_t>(std::abs(record.bin1_start - record.bin2_start) /
                                                resolution());
        assert(i < expected.size());
        return float(expected[i]);
    }();

    if (matrixType() == MatrixType::expected) {
        record.count = expectedCount;
        return;
    }

    assert(matrixType() == MatrixType::oe);
    record.count /= expectedCount;
}

inline void MatrixZoomData::readBlockOfInteractionsV6(BinaryBuffer &src,
                                                      std::vector<contactRecord> &dest) {
    assert(src.i == sizeof(std::int32_t));

    constexpr auto recordSize = sizeof(std::int32_t) + sizeof(std::int32_t) + sizeof(float);

    const auto srcSize = sizeof(std::int32_t) + (src.buffer.size() * sizeof(char));
    const auto destSize = recordSize * dest.size();

    if (srcSize != destSize) {
        throw std::runtime_error(fmt::format(
            FMT_STRING("binary buffer appears to be corrupted: expected {}B, found {}B"), destSize,
            srcSize));
    }

    std::generate(dest.begin(), dest.end(), [&]() {
        // clang-format off
        return contactRecord{src.read<std::int32_t>(),
                             src.read<std::int32_t>(),
                             src.read<float>()};
        // clang-format on
    });
    return;
}

inline void MatrixZoomData::readBlockOfInteractions(indexEntry idx,
                                                    std::vector<contactRecord> &buffer) {
    buffer.clear();
    if (idx.size <= 0) {
        return;
    }

    _fs->readAndInflate(idx, _buffer.buffer);
    _buffer.i = 0;

    const auto nRecords = static_cast<std::size_t>(_buffer.read<std::int32_t>());
    buffer.resize(nRecords);

    if (_fs->version() == 6) {
        return readBlockOfInteractionsV6(_buffer, buffer);
    }

    const auto bin1Offset = _buffer.read<std::int32_t>();
    const auto bin2Offset = _buffer.read<std::int32_t>();

    const auto i16Counts = _buffer.read<char>() == 0;

    auto readUseShortBinFlag = [&]() {
        if (_fs->version() > 8) {
            return _buffer.read<char>() == 0;
        }
        return true;
    };

    const auto i16Bin1 = readUseShortBinFlag();
    const auto i16Bin2 = readUseShortBinFlag();

    const auto type = static_cast<std::int8_t>(_buffer.read<char>());
    if (type != 1 && type != 2) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("uknown interaction type \"{}\". Supported types: 1, 2"), type));
    }

    switch (type) {
        case 1:
            readBlockOfInteractionsType1Dispatcher(i16Bin1, i16Bin2, i16Counts, bin1Offset,
                                                   bin2Offset, _buffer, buffer);
            return;
        case 2:
            if (i16Counts) {
                readBlockOfInteractionsType2<std::int16_t>(bin1Offset, bin2Offset, _buffer, buffer);
                return;
            }
            readBlockOfInteractionsType2<float>(bin1Offset, bin2Offset, _buffer, buffer);
            return;
        default:
            assert(false);
            std::abort();
    }
}

inline void MatrixZoomData::readBlockOfInteractionsType1Dispatcher(
    bool i16Bin1, bool i16Bin2, bool i16Counts, std::int32_t bin1Offset, std::int32_t bin2Offset,
    BinaryBuffer &src, std::vector<contactRecord> &dest) noexcept {
    using BS = std::int16_t;  // Short type for bins
    using CS = std::int16_t;  // Short type for count

    using BL = std::int32_t;  // Long type for bins
    using CL = float;         // Long type for count

    if (i16Bin1 && i16Bin2 && i16Counts) {
        readBlockOfInteractionsType1<BS, BS, CS>(bin1Offset, bin2Offset, src, dest);
        return;
    }
    if (!i16Bin1 && i16Bin2 && i16Counts) {
        readBlockOfInteractionsType1<BL, BS, CS>(bin1Offset, bin2Offset, src, dest);
        return;
    }
    if (i16Bin1 && !i16Bin2 && i16Counts) {
        readBlockOfInteractionsType1<BS, BL, CS>(bin1Offset, bin2Offset, src, dest);
        return;
    }
    if (i16Bin1 && i16Bin2 && !i16Counts) {
        readBlockOfInteractionsType1<BS, BS, CL>(bin1Offset, bin2Offset, src, dest);
        return;
    }
    if (!i16Bin1 && !i16Bin2 && i16Counts) {
        readBlockOfInteractionsType1<BL, BL, CS>(bin1Offset, bin2Offset, src, dest);
        return;
    }
    if (!i16Bin1 && i16Bin2 && !i16Counts) {
        readBlockOfInteractionsType1<BL, BS, CL>(bin1Offset, bin2Offset, src, dest);
        return;
    }
    if (i16Bin1 && !i16Bin2 && !i16Counts) {
        readBlockOfInteractionsType1<BS, BL, CL>(bin1Offset, bin2Offset, src, dest);
        return;
    }
    assert(!i16Bin1 && !i16Bin2 && !i16Counts);
    readBlockOfInteractionsType1<BL, BL, CL>(bin1Offset, bin2Offset, src, dest);
}

template <typename Bin1Type, typename Bin2Type, typename CountType>
inline void MatrixZoomData::readBlockOfInteractionsType1(
    std::int32_t bin1Offset, std::int32_t bin2Offset, BinaryBuffer &src,
    std::vector<contactRecord> &dest) noexcept {
    using i16 = std::int16_t;
    using i32 = std::int32_t;
    using f32 = float;
    static_assert(std::is_same<i16, Bin1Type>::value || std::is_same<i32, Bin1Type>::value, "");
    static_assert(std::is_same<i16, Bin2Type>::value || std::is_same<i32, Bin2Type>::value, "");
    static_assert(std::is_same<i16, CountType>::value || std::is_same<f32, CountType>::value, "");

    constexpr auto expectedOffsetV7 = (3 * sizeof(i32)) + (2 * sizeof(char));
    constexpr auto expectedOffsetV8plus = expectedOffsetV7 + (2 * sizeof(char));
    std::ignore = expectedOffsetV7;
    std::ignore = expectedOffsetV8plus;
    assert(src.i == expectedOffsetV7 || src.i == expectedOffsetV8plus);

    const auto expectedNumRecords = dest.size();
    dest.clear();
    const auto numRows = static_cast<i32>(src.read<Bin2Type>());
    for (i32 i = 0; i < numRows; ++i) {
        const auto bin2 = bin2Offset + static_cast<i32>(src.read<Bin2Type>());

        const auto numCols = static_cast<i32>(src.read<Bin1Type>());
        for (i32 j = 0; j < numCols; ++j) {
            const auto bin1 = bin1Offset + static_cast<i32>(src.read<Bin1Type>());

            const auto counts = static_cast<f32>(src.read<CountType>());
            dest.push_back(contactRecord{bin1, bin2, counts});
        }
    }

    std::ignore = expectedNumRecords;
    assert(expectedNumRecords == dest.size());
}

template <typename CountType>
inline void MatrixZoomData::readBlockOfInteractionsType2(
    std::int32_t bin1Offset, std::int32_t bin2Offset, BinaryBuffer &src,
    std::vector<contactRecord> &dest) noexcept {
    using i16 = std::int16_t;
    using i32 = std::int32_t;
    using f32 = float;
    static_assert(std::is_same<i16, CountType>::value || std::is_same<f32, CountType>::value, "");

    const auto nPts = src.read<i32>();
    const auto w = static_cast<i32>(src.read<i16>());

    constexpr auto i16Sentinel = (std::numeric_limits<i16>::lowest)();
    constexpr auto i16Counts = std::is_same<i16, CountType>::value;

    auto isValid = [&](CountType n) {
        return (i16Counts && static_cast<i16>(n) != i16Sentinel) &&
               (!i16Counts && !std::isnan(static_cast<f32>(n)));
    };

    dest.reserve(static_cast<std::size_t>(nPts));
    dest.clear();
    for (i32 i = 0; i < nPts; ++i) {
        const auto count = src.read<CountType>();
        if (!isValid(count)) {
            continue;
        }
        const auto row = i / w;
        const auto col = i - row * w;
        const auto bin1 = bin1Offset + col;
        const auto bin2 = bin2Offset + row;

        dest.emplace_back(contactRecord{bin1, bin2, static_cast<f32>(count)});
    }
}

inline BlockMap MatrixZoomData::readBlockMap(HiCFileStream &fs, const HiCFooter &footer) {
    BlockMap buffer{};
    fs.readBlockMap(footer.fileOffset(), footer.chrom1(), footer.chrom2(), footer.unit(),
                    footer.resolution(), buffer);
    return buffer;
}

inline void MatrixZoomData::readBlockNumbers(std::int64_t bin1, std::int64_t bin2,
                                             std::int64_t bin3, std::int64_t bin4,
                                             std::set<std::size_t> &buffer) const {
    const auto blockBinCount = _blockMap.blockBinCount;
    const auto blockColumnCount = _blockMap.blockColumnCount;

    const auto col1 = bin1 / blockBinCount;
    const auto col2 = (bin2 + 1) / blockBinCount;
    const auto row1 = bin3 / blockBinCount;
    const auto row2 = (bin4 + 1) / blockBinCount;

    // check region part that overlaps with lower left triangle but only if intrachromosomal
    const auto checkLowerLeftTri = isIntra();
    buffer.clear();
    // first check the upper triangular matrixType
    for (auto row = row1; row <= row2; ++row) {
        for (auto col = col1; col <= col2; ++col) {
            buffer.insert(static_cast<std::size_t>(row * blockColumnCount + col));
            if (checkLowerLeftTri) {
                buffer.insert(static_cast<std::size_t>(col * blockColumnCount + row));
            }
        }
    }
}

inline void MatrixZoomData::readBlockNumbersV9Intra(std::int64_t bin1, std::int64_t bin2,
                                                    std::int64_t bin3, std::int64_t bin4,
                                                    std::set<std::size_t> &buffer) const {
    const auto blockBinCount = _blockMap.blockBinCount;
    const auto blockColumnCount = _blockMap.blockColumnCount;

    const auto translatedLowerPAD = (bin1 + bin3) / 2 / blockBinCount;
    const auto translatedHigherPAD = (bin2 + bin4) / 2 / blockBinCount + 1;
    const auto translatedNearerDepth = static_cast<std::int64_t>(
        std::log2(1.0 + double(std::abs(bin1 - bin4)) / std::sqrt(2.0) / blockBinCount));
    const auto translatedFurtherDepth = static_cast<std::int64_t>(
        std::log2(1.0 + double(std::abs(bin2 - bin3)) / std::sqrt(2.0) / blockBinCount));

    // because code above assumes above diagonal; but we could be below diagonal
    const auto nearerDepth = [&]() -> std::int64_t {
        if ((bin1 > bin4 && bin2 < bin3) || (bin2 > bin3 && bin1 < bin4)) {
            return 0;
        }
        return std::min(translatedNearerDepth, translatedFurtherDepth);
    }();

    // +1; integer divide rounds down
    const auto furtherDepth = std::max(translatedNearerDepth, translatedFurtherDepth) + 1;

    buffer.clear();
    for (auto depth = nearerDepth; depth <= furtherDepth; ++depth) {
        for (auto pad = translatedLowerPAD; pad <= translatedHigherPAD; ++pad) {
            buffer.insert(static_cast<std::size_t>(depth * blockColumnCount + pad));
        }
    }
}

}  // namespace hicxx::internal
