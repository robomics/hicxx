// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <fmt/compile.h>
#include <fmt/format.h>

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "straw/straw.hpp"

using namespace hicxx;

static std::int32_t getChromSize(const HiCFile& hic, const std::string& chrom_name) {
    auto it = hic.chromosomes().find(chrom_name);
    if (it == hic.chromosomes().end()) {
        throw std::runtime_error(
            fmt::format(FMT_STRING("Unable to find a chromosome named \"{}\""), chrom_name));
    }
    return static_cast<std::int32_t>(it->second.length);
}

int main(int argc, char** argv) noexcept {
    try {
        if (argc != 7 && argc != 8) {
            fmt::print(
                stderr,
                FMT_STRING("Incorrect arguments\n"
                           "Usage: straw [observed/oe/expected] <NONE/VC/VC_SQRT/KR> <hicFile(s)> "
                           "<chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>\n"));
            return 1;
        }

        auto matrixType = MatrixType::observed;
        auto* args = argv + 1;
        if (argc == 8) {
            matrixType = ParseMatrixTypeStr(*args++);
        }
        const auto norm = ParseNormStr(*args++);
        const std::string url(*args++);
        auto coord1 = GenomicCoordinates::fromString(*args++);
        auto coord2 = GenomicCoordinates::fromString(*args++);
        const auto unit = ParseUnitStr(*args++);
        const auto resolution = std::stoi(*args++);

        HiCFile hic(url);

        if (coord1.end == 0) {
            coord1.end = getChromSize(hic, coord1.chrom);
        }
        if (coord2.end == 0) {
            coord2.end = getChromSize(hic, coord2.chrom);
        }

        auto selector =
            hic.getMatrixZoomData(coord1.chrom, coord2.chrom, matrixType, norm, unit, resolution);
        std::vector<contactRecord> buffer{};
        selector.fetch(coord1.start, coord1.end, coord2.start, coord2.end, buffer);

        for (const auto& record : buffer) {
            const auto start1 = record.bin1_start;
            const auto start2 = record.bin2_start;
            fmt::print(FMT_COMPILE("{}\t{}\t{}\n"), start1, start2, record.count);
        }
    } catch (const std::exception& e) {
        const auto* url = argc == 8 ? *(argv + 3) : *(argv + 2);
        fmt::print(
            stderr,
            FMT_STRING("straw encountered the following error while processing file \"{}\": {}\n"),
            url, e.what());
        return 1;
    } catch (...) {
        fmt::print(stderr, FMT_STRING("straw encountered an unknown error!\n"));
        return 1;
    }
}
