// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#include <cstdint>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "straw/internal/common.hpp"
#include "straw/internal/filestream.hpp"
#include "straw/internal/hic_file_stream.hpp"
#include "straw/internal/hic_footer.hpp"
#include "straw/internal/hic_header.hpp"
#include "straw/internal/hic_matrix_selector.hpp"

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
