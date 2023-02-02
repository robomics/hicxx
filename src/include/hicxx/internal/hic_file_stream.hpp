// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#pragma once

#ifdef HICXX_USE_ZLIBNG
#include <zlib-ng.h>
#else
#include <zlib.h>
#endif

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "hicxx/internal/common.hpp"
#include "hicxx/internal/filestream.hpp"
#include "hicxx/internal/hic_footer.hpp"
#include "hicxx/internal/hic_header.hpp"

namespace hicxx::internal {

struct BlockMap {
    std::map<std::size_t, indexEntry> blocks{};
    std::int32_t blockBinCount{};
    std::int32_t blockColumnCount{};
    double sumCount{};
};

class HiCFileStream {
#ifdef HICXX_USE_ZLIBNG
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
    [[nodiscard]] inline const std::string &url() const noexcept;
    [[nodiscard]] const HiCHeader &header() const noexcept;

    [[nodiscard]] bool isLocal() const noexcept;
    [[nodiscard]] bool isRemote() const noexcept;

    [[nodiscard]] std::int32_t version() const noexcept;

    // reads the footer given a pair of chromosomes, norm, unit (BP or FRAG) and resolution.
    [[nodiscard]] HiCFooter readFooter(std::int32_t chromId1, std::int32_t chromId2,
                                       MatrixType matrixType, NormalizationMethod norm,
                                       MatrixUnit unit, std::int32_t wantedResolution);

    [[nodiscard]] static MatrixType readMatrixType(filestream::FileStream &fs, std::string &buff);
    [[nodiscard]] static NormalizationMethod readNormalizationMethod(filestream::FileStream &fs,
                                                                     std::string &buff);
    [[nodiscard]] static MatrixUnit readMatrixUnit(filestream::FileStream &fs, std::string &buff);

    void readBlockMap(std::int64_t fileOffset, const chromosome &chrom1, const chromosome &chrom2,
                      MatrixUnit wantedUnit, std::int64_t wantedResolution, BlockMap &buffer);
    void readAndInflate(indexEntry idx, std::string &plainTextBuffer);

   private:
    [[nodiscard]] static filestream::FileStream openStream(std::string url);
    // reads the header, storing the positions of the normalization vectors and returning the
    // masterIndexPosition pointer
    [[nodiscard]] static HiCHeader readHeader(filestream::FileStream &fs);

    [[nodiscard]] std::vector<double> readExpectedVector(std::int64_t nValues);
    [[nodiscard]] std::vector<double> readNormalizationFactors(std::int32_t wantedChrom);
    void applyNormalizationFactors(std::vector<double> &expectedValues,
                                   const std::vector<double> &normFactors);
    [[nodiscard]] std::vector<double> readNormalizationVector(indexEntry cNormEntry,
                                                              std::size_t numValuesExpected);

    void discardExpectedVector(std::int64_t nValues);
    void discardNormalizationFactors(std::int32_t wantedChrom);

    [[nodiscard]] MatrixType readMatrixType();
    [[nodiscard]] NormalizationMethod readNormalizationMethod();
    [[nodiscard]] MatrixUnit readMatrixUnit();

    [[nodiscard]] std::int64_t readNValues();
    [[nodiscard]] bool checkMagicString();
    [[nodiscard]] static bool checkMagicString(filestream::FileStream &fs);
    [[nodiscard]] std::int64_t masterOffset() const noexcept;

    [[nodiscard]] static auto initZStream() -> ZStream;
};
}  // namespace hicxx::internal

#include "../../../hic_file_stream_impl.hpp"
