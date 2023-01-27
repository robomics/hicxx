// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <filestream/filestream.hpp>
#include <memory>
#include <set>
#include <string>
#include <unordered_set>
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
};

struct HiCFooter {
    std::string url{};
    MatrixType matrixType{MatrixType::observed};
    Normalization normalization{Normalization::NONE};
    Unit unit{Unit::BP};
    std::int32_t resolution{-1};
    std::int64_t fileOffset{-1};
    indexEntry c1NormEntry{-1, -1};
    indexEntry c2NormEntry{-1, -1};
    std::vector<double> expectedValues{};

    constexpr explicit operator bool() const noexcept;
    bool operator==(const HiCFooter &other) const noexcept;
    bool operator!=(const HiCFooter &other) const noexcept;
};
}  // namespace internal

template <>
struct std::hash<internal::HiCHeader> {
    inline std::size_t operator()(internal::HiCHeader const &h) const noexcept {
        return internal::hash_combine(0, h.url, h.masterIndexOffset);
    }
};

template <>
struct std::hash<internal::HiCFooter> {
    inline std::size_t operator()(internal::HiCFooter const &f) const noexcept {
        return internal::hash_combine(0, f.url, f.matrixType, f.normalization, f.unit,
                                      f.resolution);
    }
};

namespace internal {
class HiCFileStream {
    filestream::FileStream fs_{};
    HiCHeader header_{};
    std::unordered_set<HiCFooter> footers_{};

   public:
    HiCFileStream() = default;
    explicit HiCFileStream(std::string url);
    inline const std::string &url() const noexcept;
    constexpr const HiCHeader &header() const noexcept;

    bool isLocal() const noexcept;
    bool isRemote() const noexcept;

    // reads the footer given a pair of chromosomes, norm, unit (BP or FRAG) and resolution.
    HiCFooter readFooter(const chromosome &chrom1, const chromosome &chrom2, MatrixType matrixType,
                         Normalization norm, Unit unit, std::int32_t wantedResolution);
    HiCFooter readFooter(const std::string &chrom1, const std::string &chrom2,
                         MatrixType matrixType, Normalization norm, Unit unit,
                         std::int32_t wantedResolution);
    HiCFooter readFooter(std::int32_t chromId1, std::int32_t chromId2, MatrixType matrixType,
                         Normalization norm, Unit unit, std::int32_t wantedResolution);
    HiCFooter readFooter(const chromosome &chrom, MatrixType matrixType, Normalization norm,
                         Unit unit, std::int32_t wantedResolution);
    HiCFooter readFooter(const std::string &chrom, MatrixType matrixType, Normalization norm,
                         Unit unit, std::int32_t wantedResolution);
    HiCFooter readFooter(std::int32_t chromId, MatrixType matrixType, Normalization norm, Unit unit,
                         std::int32_t wantedResolution);

    void removeCachedFooter(const HiCFooter &footer);
    void purgeFooterCache();

   private:
    static filestream::FileStream openStream(std::string url);
    // reads the header, storing the positions of the normalization vectors and returning the
    // masterIndexPosition pointer
    static HiCHeader readHeader(filestream::FileStream &fs);

    std::int32_t version() const noexcept;
    void readCompressedBytes(indexEntry idx, std::string &buffer);

    std::vector<double> readExpectedVector(std::int64_t nValues);
    std::vector<double> readNormalizationFactors(std::int32_t wantedChrom);
    void applyNormalizationFactors(std::vector<double> &expectedValues,
                                   const std::vector<double> &normFactors);

    void discardExpectedVector(std::int64_t nValues);
    void discardNormalizationFactors(std::int32_t wantedChrom);

    MatrixType readMatrixType();
    Normalization readNormalization();
    Unit readUnit();

    std::int64_t readNValues();
    bool checkMagicString();
    static bool checkMagicString(filestream::FileStream &fs);
    std::int64_t masterOffset() const noexcept;
};

class MatrixZoomData {
    std::shared_ptr<HiCFileStream> stream;
    std::int64_t myFilePos = 0LL;
    std::vector<double> expectedValues;
    std::vector<double> c1Norm;
    std::vector<double> c2Norm;
    std::int32_t c1 = 0;
    std::int32_t c2 = 0;
    MatrixType matrixType;
    Normalization norm;
    std::int32_t version = 0;
    std::int32_t resolution = 0;
    std::int32_t numBins1 = 0;
    std::int32_t numBins2 = 0;
    float sumCounts;
    std::int32_t blockBinCount;
    std::int32_t blockColumnCount;
    std::map<std::int32_t, indexEntry> blockMap;
    double avgCount;

   public:
    MatrixZoomData() = delete;
    MatrixZoomData(std::shared_ptr<HiCFileStream> stream_, const chromosome &chrom1,
                   const chromosome &chrom2, const std::string &matrixType, const std::string &norm,
                   const std::string &unit, std::int32_t resolution, std::int32_t &version,
                   std::int64_t &master, std::int64_t &totalFileSize);

    MatrixZoomData(std::shared_ptr<HiCFileStream> stream_, const chromosome &chrom1,
                   const chromosome &chrom2, MatrixType matrixType, Normalization norm, Unit unit,
                   std::int32_t resolution, std::int32_t &version, std::int64_t &master,
                   std::int64_t &totalFileSize);

    std::vector<contactRecord> getRecords(std::int64_t gx0, std::int64_t gx1, std::int64_t gy0,
                                          std::int64_t gy1);

    std::vector<std::vector<float>> getRecordsAsMatrix(std::int64_t gx0, std::int64_t gx1,
                                                       std::int64_t gy0, std::int64_t gy1);

    std::int64_t getNumberOfTotalRecords();
    bool isIntra() const noexcept;

   private:
    static std::vector<double> readNormalizationVectorFromFooter(indexEntry cNormEntry,
                                                                 std::int32_t &version,
                                                                 const std::string &fileName);

    static bool isInRange(std::int32_t r, std::int32_t c, std::int32_t numRows,
                          std::int32_t numCols);

    std::set<std::int32_t> getBlockNumbers(std::int64_t *regionIndices) const;

    std::vector<double> getNormVector(std::int32_t index);

    std::vector<double> getExpectedValues();
};

}  // namespace internal

class HiCFile {
    std::shared_ptr<internal::HiCFileStream> stream;
    std::int64_t master = 0LL;
    ChromosomeMap chromosomes;
    std::string genomeID;
    std::int32_t numChromosomes = 0;
    std::int32_t version = 0;
    std::int64_t nviPosition = 0LL;
    std::int64_t nviLength = 0LL;
    std::vector<std::int32_t> resolutions;
    std::int64_t totalFileSize;

   public:
    explicit HiCFile(std::string fileName_);

    const std::string &url() const noexcept;
    const std::string &getGenomeID() const noexcept;

    const std::vector<std::int32_t> &getResolutions() const noexcept;

    std::vector<chromosome> getChromosomes() const;
    auto getChromosomeMap() const noexcept -> const ChromosomeMap &;

    internal::MatrixZoomData getMatrixZoomData(const std::string &chr1, const std::string &chr2,
                                               const std::string &matrixType,
                                               const std::string &norm, const std::string &unit,
                                               std::int32_t resolution);
    internal::MatrixZoomData getMatrixZoomData(const std::string &chr1, const std::string &chr2,
                                               MatrixType matrixType, Normalization norm, Unit unit,
                                               std::int32_t resolution);

   private:
    static std::int64_t readTotalFileSize(const std::string &url);

    auto readHeader(std::istream &fin, std::int64_t &masterIndexPosition, std::string &genomeID,
                    std::int32_t &numChromosomes, std::int32_t &version, std::int64_t &nviPosition,
                    std::int64_t &nviLength) -> ChromosomeMap;
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

std::vector<contactRecord> straw(MatrixType matrixType, Normalization norm,
                                 const std::string &fname, const std::string &chr1loc,
                                 const std::string &chr2loc, Unit unit, std::int32_t binsize);

std::int64_t getNumRecordsForFile(const std::string &filename, std::int32_t binsize,
                                  bool interOnly);

// #include "../../hic_file_impl.hpp"
#include "../../hic_file_stream_impl.hpp"
// #include "../../matrix_zoom_data_impl.hpp"
// #include "../../straw_impl.hpp"
