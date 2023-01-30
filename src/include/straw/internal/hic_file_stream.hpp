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
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "straw/internal/common.hpp"
#include "straw/internal/hic_footer.hpp"
#include "straw/internal/hic_header.hpp"

namespace internal {

struct BlockMap {
    std::map<std::int32_t, indexEntry> blocks{};
    std::int32_t blockBinCount{};
    std::int32_t blockColumnCount{};
    double sumCount{};
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
}  // namespace internal

#include "../../../hic_file_stream_impl.hpp"
