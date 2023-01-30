// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef STRAW_USE_ZLIBNG
#include <zlib-ng.h>
#else
#include <zlib.h>
#endif

#include <cstdint>
#include <filestream/filestream.hpp>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "straw/internal/common.h"

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

template <>
struct std::hash<internal::HiCHeader> {
    inline std::size_t operator()(internal::HiCHeader const &h) const noexcept {
        return internal::hash_combine(0, h.url, h.masterIndexOffset);
    }
};

template <>
struct std::hash<internal::HiCFooterMetadata> {
    inline std::size_t operator()(internal::HiCFooterMetadata const &m) const noexcept {
        return internal::hash_combine(0, m.url, m.matrixType, m.normalization, m.unit, m.resolution,
                                      m.chrom1, m.chrom2);
    }
};

namespace internal {

struct BlockMap {
    std::map<std::int32_t, indexEntry> blocks{};
    std::int32_t blockBinCount{};
    std::int32_t blockColumnCount{};
    double sumCount{};
};

struct BinaryBuffer {
    std::string buffer{};
    std::size_t i{};

    template <typename T, typename std::enable_if<std::is_fundamental<T>::value>::type * = nullptr>
    T read();
};

class HiCFileStream {
#ifdef STRAW_USE_ZLIBNG
    using ZStream = UniquePtrWithDeleter<zng_stream>;
#else
    using ZStream = UniquePtrWithDeleter<z_stream>;
#endif
    std::shared_ptr<filestream::FileStream> _fs{};
    std::shared_ptr<const HiCHeader> _header{};
    std::string _strbuff{};
    ZStream _zlibstream{initZStream()};

   public:
    HiCFileStream() = default;
    explicit HiCFileStream(std::string url);
    inline const std::string &url() const noexcept;
    const HiCHeader &header() const noexcept;

    bool isLocal() const noexcept;
    bool isRemote() const noexcept;

    std::int32_t version() const noexcept;

    // reads the footer given a pair of chromosomes, norm, unit (BP or FRAG) and resolution.
    HiCFooter readFooter(std::int32_t chromId1, std::int32_t chromId2, MatrixType matrixType,
                         NormalizationMethod norm, MatrixUnit unit, std::int32_t wantedResolution);

    static MatrixType readMatrixType(filestream::FileStream &fs, std::string &buff);
    static NormalizationMethod readNormalizationMethod(filestream::FileStream &fs,
                                                       std::string &buff);
    static MatrixUnit readMatrixUnit(filestream::FileStream &fs, std::string &buff);

    void readBlockMap(std::int64_t fileOffset, const chromosome &chrom1, const chromosome &chrom2,
                      MatrixUnit wantedUnit, std::int32_t wantedResolution, BlockMap &buffer);
    void readAndInflate(indexEntry idx, std::string &plainTextBuffer);

   private:
    static filestream::FileStream openStream(std::string url);
    // reads the header, storing the positions of the normalization vectors and returning the
    // masterIndexPosition pointer
    static HiCHeader readHeader(filestream::FileStream &fs);

    std::vector<double> readExpectedVector(std::int64_t nValues);
    std::vector<double> readNormalizationFactors(std::int32_t wantedChrom);
    void applyNormalizationFactors(std::vector<double> &expectedValues,
                                   const std::vector<double> &normFactors);
    std::vector<double> readNormalizationVector(indexEntry cNormEntry);

    void discardExpectedVector(std::int64_t nValues);
    void discardNormalizationFactors(std::int32_t wantedChrom);

    MatrixType readMatrixType();
    NormalizationMethod readNormalizationMethod();
    MatrixUnit readMatrixUnit();

    std::int64_t readNValues();
    bool checkMagicString();
    static bool checkMagicString(filestream::FileStream &fs);
    std::int64_t masterOffset() const noexcept;

    static auto initZStream() -> ZStream;
};

class MatrixZoomData {
    std::shared_ptr<HiCFileStream> _fs;
    std::shared_ptr<const HiCFooter> _footer;
    BlockMap _blockMap{};
    std::set<std::int32_t> _blockNumberBuff{};
    std::vector<contactRecord> _contactRecordBuff{};
    BinaryBuffer _buffer{};

   public:
    MatrixZoomData() = delete;
    MatrixZoomData(std::shared_ptr<HiCFileStream> fs, std::shared_ptr<const HiCFooter> footer);

    const chromosome &chrom1() const noexcept;
    const chromosome &chrom2() const noexcept;

    std::int32_t resolution() const noexcept;
    MatrixType matrixType() const noexcept;
    NormalizationMethod normalizationMethod() const noexcept;
    MatrixUnit matrixUnit() const noexcept;

    std::int32_t numBins1() const noexcept;
    std::int32_t numBins2() const noexcept;

    bool isIntra() const noexcept;
    bool isInter() const noexcept;

    const std::vector<double> &chrom1Norm() const noexcept;
    const std::vector<double> &chrom2Norm() const noexcept;

    inline double avgCount() const;

    void fetch(std::vector<contactRecord> &buffer);
    void fetch(const std::string &coord, std::vector<contactRecord> &buffer);
    void fetch(std::int64_t start, std::int64_t end, std::vector<contactRecord> &buffer);

    void fetch(const std::string &coord1, const std::string &coord2,
               std::vector<contactRecord> &buffer);
    void fetch(std::int64_t start1, std::int64_t end1, std::int64_t start2, std::int64_t end2,
               std::vector<contactRecord> &buffer);

    void fetch(std::vector<std::vector<float>> &buffer);
    void fetch(const std::string &coord, std::vector<std::vector<float>> &buffer);
    void fetch(std::int64_t start, std::int64_t end, std::vector<std::vector<float>> &buffer);

    void fetch(const std::string &coord1, const std::string &coord2,
               std::vector<std::vector<float>> &buffer);
    void fetch(std::int64_t start1, std::int64_t end1, std::int64_t start2, std::int64_t end2,
               std::vector<std::vector<float>> &buffer);

   private:
    static BlockMap readBlockMap(HiCFileStream &fs, const HiCFooter &footer);

    void readBlockNumbers(std::int64_t bin1, std::int64_t bin2, std::int64_t bin3,
                          std::int64_t bin4, std::set<std::int32_t> &buffer) const;
    void readBlockNumbersV9Intra(std::int64_t bin1, std::int64_t bin2, std::int64_t bin3,
                                 std::int64_t bin4, std::set<std::int32_t> &buffer) const;
    void readBlockOfInteractions(indexEntry idx, std::vector<contactRecord> &buffer);
    void processInteraction(contactRecord &record, std::int32_t pos1, std::int32_t pos2);
    static void readBlockOfInteractionsV6(BinaryBuffer &src, std::vector<contactRecord> &dest);

    static void readBlockOfInteractionsType1Dispatcher(bool i16Bin1, bool i16Bin2, bool i16Counts,
                                                       std::int32_t bin1Offset,
                                                       std::int32_t bin2Offset, BinaryBuffer &src,
                                                       std::vector<contactRecord> &dest) noexcept;
    template <typename Bin1Type, typename Bin2Type, typename CountType>
    static void readBlockOfInteractionsType1(std::int32_t bin1Offset, std::int32_t bin2Offset,
                                             BinaryBuffer &src,
                                             std::vector<contactRecord> &dest) noexcept;

    template <typename CountType>
    static void readBlockOfInteractionsType2(std::int32_t bin1Offset, std::int32_t bin2Offset,
                                             BinaryBuffer &src,
                                             std::vector<contactRecord> &dest) noexcept;
};

}  // namespace internal

class HiCFile {
    // clang-format off
    using FooterCacheT =
        std::unordered_map<internal::HiCFooterMetadata,
                           std::shared_ptr<const internal::HiCFooter>>;
    // clang-format on
    std::shared_ptr<internal::HiCFileStream> _fs{};
    FooterCacheT _footers{};

   public:
    explicit HiCFile(std::string url_);

    const std::string &url() const noexcept;
    const std::string &name() const noexcept;
    std::int32_t version() const noexcept;
    const ChromosomeMap &chromosomes() const noexcept;
    const std::string &genomeID() const noexcept;
    const std::vector<std::int32_t> &resolutions() const noexcept;

    internal::MatrixZoomData getMatrixZoomData(const chromosome &chrom, MatrixType matrixType,
                                               NormalizationMethod norm, MatrixUnit unit,
                                               std::int32_t resolution);
    internal::MatrixZoomData getMatrixZoomData(const std::string &chromName, MatrixType matrixType,
                                               NormalizationMethod norm, MatrixUnit unit,
                                               std::int32_t resolution);
    internal::MatrixZoomData getMatrixZoomData(std::int32_t chromId, MatrixType matrixType,
                                               NormalizationMethod norm, MatrixUnit unit,
                                               std::int32_t resolution);

    internal::MatrixZoomData getMatrixZoomData(const chromosome &chrom1, const chromosome &chrom2,
                                               MatrixType matrixType, NormalizationMethod norm,
                                               MatrixUnit unit, std::int32_t resolution);
    internal::MatrixZoomData getMatrixZoomData(const std::string &chromName1,
                                               const std::string &chromName2, MatrixType matrixType,
                                               NormalizationMethod norm, MatrixUnit unit,
                                               std::int32_t resolution);
    internal::MatrixZoomData getMatrixZoomData(std::int32_t chromId1, std::int32_t chromId2,
                                               MatrixType matrixType, NormalizationMethod norm,
                                               MatrixUnit unit, std::int32_t resolution);

    std::size_t numCachedFooters() const noexcept;
    void purgeFooterCache();

   private:
    std::shared_ptr<const internal::HiCFooter> getFooter(std::int32_t chromId1,
                                                         std::int32_t chromId2,
                                                         MatrixType matrixType,
                                                         NormalizationMethod norm, MatrixUnit unit,
                                                         std::int32_t resolution);
};

struct GenomicCoordinates {
    std::string chrom;
    std::int32_t start;
    std::int32_t end;

    static GenomicCoordinates fromString(std::string coord, bool noChromName = false);
};

std::map<std::int32_t, indexEntry> readMatrixZoomData(std::istream &fin, const std::string &myunit,
                                                      std::int32_t mybinsize, float &mySumCounts,
                                                      std::int32_t &myBlockBinCount,
                                                      std::int32_t &myBlockColumnCount,
                                                      bool &found);

std::vector<double> readNormalizationVector(std::istream &fin, indexEntry entry);

std::vector<contactRecord> straw(const std::string &matrixType, const std::string &norm,
                                 const std::string &fname, const std::string &chr1loc,
                                 const std::string &chr2loc, const std::string &unit,
                                 std::int32_t binsize);

std::vector<contactRecord> straw(MatrixType matrixType, NormalizationMethod norm,
                                 const std::string &fname, const std::string &chr1loc,
                                 const std::string &chr2loc, MatrixUnit unit, std::int32_t binsize);

std::int64_t getNumRecordsForFile(const std::string &filename, std::int32_t binsize,
                                  bool interOnly);

#include "../../hic_file_impl.hpp"
#include "../../hic_file_stream_impl.hpp"
#include "../../matrix_zoom_data_impl.hpp"
