// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <fmt/format.h>
#include <zlib.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <ios>
#include <istream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "straw/internal/common.h"

namespace internal {

inline MatrixZoomData::MatrixZoomData(std::shared_ptr<HiCFileStream> fs,
                                      std::shared_ptr<const HiCFooter> footer)
    : _fs(std::move(fs)), _footer(std::move(footer)), _blockMap(readBlockMap(*_fs, *_footer)) {}

inline const chromosome &MatrixZoomData::chrom1() const noexcept { return _footer->chrom1(); }

inline const chromosome &MatrixZoomData::chrom2() const noexcept { return _footer->chrom2(); }

inline std::int32_t MatrixZoomData::resolution() const noexcept { return _footer->resolution(); }

inline MatrixType MatrixZoomData::matrixType() const noexcept { return _footer->matrixType(); }

inline NormalizationMethod MatrixZoomData::normalizationMethod() const noexcept {
    return _footer->normalization();
}

inline MatrixUnit MatrixZoomData::matrixUnit() const noexcept { return _footer->unit(); }

inline std::int32_t MatrixZoomData::numBins1() const noexcept {
    return (chrom1().length + resolution() - 1) / resolution();
}

inline std::int32_t MatrixZoomData::numBins2() const noexcept {
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

inline double MatrixZoomData::avgCount() const noexcept {
    if (isIntra()) {
        return 0;
    }
    return (_sumCount / numBins1()) / numBins2();  // <= trying to avoid overflows
}

inline void MatrixZoomData::fetch(const std::string &coord, std::vector<contactRecord> &buffer) {
    return fetch(coord, coord, buffer);
}

inline void MatrixZoomData::fetch(const std::string &coord1, const std::string &coord2,
                                  std::vector<contactRecord> &buffer) {
    auto parse_coords =
        [](const std::string &coord) {  // TODO this should be a free-standing function
            try {
                auto sep_pos = std::min(coord.find('-'), coord.find(':'));
                if (sep_pos == std::string::npos) {
                    throw std::invalid_argument("missing position delimiter");
                }

                const auto pos1 = std::stoll(coord.substr(0, sep_pos));
                const auto pos2 = std::stoll(coord.substr(sep_pos + 1));
                if (pos1 > pos2) {
                    throw std::invalid_argument("pos1 > pos2");
                }
                return std::make_pair(pos1, pos2);
            } catch (const std::exception &e) {
                throw std::runtime_error(fmt::format(
                    FMT_STRING("Unable to parse genomic coordinates \"{}\": {}"), coord, e.what()));
            }
        };

    const auto coord1_ = parse_coords(coord1);
    const auto coord2_ = parse_coords(coord2);

    return fetch(coord1_.first, coord1_.second, coord2_.first, coord2_.second, buffer);
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
            const auto pos1 = record.binX * resolution();
            const auto pos2 = record.binY * resolution();

            if (pos1 < start1 || pos1 > end1 || pos2 < start2 || pos2 > end2) {
                continue;
            }

            if (isIntra() && (pos2 < start1 || pos2 > end1 || pos1 < start2 || pos1 > end2)) {
                continue;
            }

            processInteraction(record, pos1, pos2);
            if (!std::isnan(record.counts) && !std::isinf(record.counts)) {
                buffer.emplace_back(std::move(record));
            }
        }
    }
}

inline void MatrixZoomData::processInteraction(contactRecord &record, std::int32_t pos1,
                                               std::int32_t pos2) {
    const auto &c1Norm = _footer->c1Norm();
    const auto &c2Norm = _footer->c2Norm();
    const auto &expected = _footer->expectedValues();

    const auto skipNormalization =
        normalizationMethod() == NormalizationMethod::NONE || matrixType() == MatrixType::expected;

    if (!skipNormalization) {
        record.counts /= c1Norm[record.binX] * c2Norm[record.binY];
    }

    if (matrixType() == MatrixType::observed) {
        return;
    }

    const auto expectedCount = [&]() {
        if (isInter()) {
            return avgCount();
        }

        const auto i = std::abs(pos1 - pos2) / resolution();
        assert(i < expected.size());
        return expected[i];
    }();

    if (matrixType() == MatrixType::expected) {
        record.counts = expectedCount;
        return;
    }

    assert(matrixType() == MatrixType::oe);
    record.counts /= expectedCount;
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
    using CS = std::int16_t;  // Short type for counts

    using BL = std::int32_t;  // Long type for bins
    using CL = float;         // Long type for counts

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
    assert(src.i == expectedOffsetV7 || src.i == expectedOffsetV8plus);

    std::size_t i = 0;
    const auto rowCount = static_cast<i32>(src.read<Bin2Type>());
    for (Bin2Type r = 0; r < rowCount; ++r) {
        const auto bin2 = bin2Offset + static_cast<i32>(src.read<Bin2Type>());

        const auto colCount = static_cast<i32>(src.read<Bin1Type>());
        for (Bin1Type c = 0; c < colCount; ++c) {
            const auto bin1 = bin1Offset + static_cast<i32>(src.read<Bin1Type>());

            const auto counts = static_cast<f32>(src.read<CountType>());
            assert(i < dest.size());
            dest[i++] = contactRecord{bin1, bin2, counts};
        }
    }

    assert(dest.size() == i);
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

    dest.reserve(nPts);
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
                                             std::set<std::int32_t> &buffer) const {
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
            buffer.insert(row * blockColumnCount + col);
            if (checkLowerLeftTri) {
                buffer.insert(col * blockColumnCount + row);
            }
        }
    }
}

inline void MatrixZoomData::readBlockNumbersV9Intra(std::int64_t bin1, std::int64_t bin2,
                                                    std::int64_t bin3, std::int64_t bin4,
                                                    std::set<std::int32_t> &buffer) const {
    const auto blockBinCount = _blockMap.blockBinCount;
    const auto blockColumnCount = _blockMap.blockColumnCount;

    const auto translatedLowerPAD = (bin1 + bin3) / 2 / blockBinCount;
    const auto translatedHigherPAD = (bin2 + bin4) / 2 / blockBinCount + 1;
    const auto translatedNearerDepth = static_cast<std::int64_t>(
        std::log2(1.0 + std::abs(bin1 - bin4) / std::sqrt(2.0) / blockBinCount));
    const auto translatedFurtherDepth = static_cast<std::int64_t>(
        std::log2(1.0 + std::abs(bin2 - bin3) / std::sqrt(2.0) / blockBinCount));

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
            buffer.insert(depth * blockColumnCount + pad);
        }
    }
}

}  // namespace internal
