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

#include "hicxx/internal/common.hpp"
#include "hicxx/internal/filestream.hpp"
#include "hicxx/internal/hic_file_stream.hpp"
#include "hicxx/internal/hic_footer.hpp"
#include "hicxx/internal/hic_header.hpp"
#include "hicxx/internal/hic_matrix_selector.hpp"

namespace hicxx {
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

    [[nodiscard]] const std::string &url() const noexcept;
    [[nodiscard]] const std::string &name() const noexcept;
    [[nodiscard]] std::int32_t version() const noexcept;
    [[nodiscard]] const ChromosomeMap &chromosomes() const noexcept;
    [[nodiscard]] const std::string &genomeID() const noexcept;
    [[nodiscard]] const std::vector<std::int32_t> &resolutions() const noexcept;

    [[nodiscard]] internal::MatrixZoomData getMatrixZoomData(const chromosome &chrom,
                                                             MatrixType matrixType,
                                                             NormalizationMethod norm,
                                                             MatrixUnit unit,
                                                             std::int32_t resolution);
    [[nodiscard]] internal::MatrixZoomData getMatrixZoomData(const std::string &chromName,
                                                             MatrixType matrixType,
                                                             NormalizationMethod norm,
                                                             MatrixUnit unit,
                                                             std::int32_t resolution);
    [[nodiscard]] internal::MatrixZoomData getMatrixZoomData(std::int32_t chromId,
                                                             MatrixType matrixType,
                                                             NormalizationMethod norm,
                                                             MatrixUnit unit,
                                                             std::int32_t resolution);

    [[nodiscard]] internal::MatrixZoomData getMatrixZoomData(
        const chromosome &chrom1, const chromosome &chrom2, MatrixType matrixType,
        NormalizationMethod norm, MatrixUnit unit, std::int32_t resolution);
    [[nodiscard]] internal::MatrixZoomData getMatrixZoomData(
        const std::string &chromName1, const std::string &chromName2, MatrixType matrixType,
        NormalizationMethod norm, MatrixUnit unit, std::int32_t resolution);
    [[nodiscard]] internal::MatrixZoomData getMatrixZoomData(
        std::int32_t chromId1, std::int32_t chromId2, MatrixType matrixType,
        NormalizationMethod norm, MatrixUnit unit, std::int32_t resolution);

    [[nodiscard]] std::size_t numCachedFooters() const noexcept;
    void purgeFooterCache();

   private:
    [[nodiscard]] std::shared_ptr<const internal::HiCFooter> getFooter(
        std::int32_t chromId1, std::int32_t chromId2, MatrixType matrixType,
        NormalizationMethod norm, MatrixUnit unit, std::int32_t resolution);
};
}  // namespace hicxx

#include "../../hic_file_impl.hpp"
