/*
  The MIT License (MIT)

  Copyright (c) 2011-2016 Broad Institute, Aiden Lab

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
*/
#include <fmt/compile.h>
#include <fmt/format.h>

#include <cassert>
#include <cstdint>
#include <string>
#include <utility>

#include "straw/straw.h"

struct GenomicCoord {
    std::string chrom{};
    std::int32_t start{};
    std::int32_t end{};

    GenomicCoord(std::string chrom_, std::int32_t start_, std::int32_t end_)
        : chrom(std::move(chrom_)), start(start_), end(end_) {
        if (start >= end) {
            throw std::runtime_error(fmt::format(
                FMT_STRING("invalid coordinate {}:{}-{}: start position >= end position"), chrom,
                start, end));
        }
        assert(start <= end);
    }

    explicit GenomicCoord(const std::string& coord) {
        auto pos = coord.find(':');
        if (pos == std::string::npos) {
            chrom = coord;
            return;
        }

        chrom = coord.substr(0, pos);
        const auto coord_ = coord.substr(pos + 1);

        pos = coord_.find('-');
        if (pos == std::string::npos) {
            pos = coord_.find(':');
        }
        if (pos == std::string::npos) {
            throw std::runtime_error(
                fmt::format(FMT_STRING("unable to parse coordinate \"{}\""), coord));
        }

        try {
            start = std::stoi(coord_.substr(0, pos));
            end = std::stoi(coord_.substr(pos + 1));
            if (start >= end) {
                throw std::runtime_error(fmt::format(
                    FMT_STRING("invalid coordinate {}: start position >= end position"), coord));
            }
        } catch (const std::exception& e) {
            throw std::runtime_error(
                fmt::format(FMT_STRING("unable to parse coordinate \"{}\": {}"), coord, e.what()));
        }
    }
};

int main(int argc, char** argv) {
    if (argc != 7 && argc != 8) {
        fmt::print(
            stderr,
            FMT_STRING("Incorrect arguments\n"
                       "Usage: straw [observed/oe/expected] <NONE/VC/VC_SQRT/KR> <hicFile(s)> "
                       "<chr1>[:x1:x2] <chr2>[:y1:y2] <BP/FRAG> <binsize>\n"));
        return 1;
    }
    try {
        auto matrixType = MatrixType::observed;
        auto* args = argv + 1;
        if (argc == 8) {
            matrixType = ParseMatrixTypeStr(*args++);
        }
        const auto norm = ParseNormStr(*args++);
        const std::string url(*args++);
        GenomicCoord coord1(*args++);
        GenomicCoord coord2(*args++);
        const auto unit = ParseUnitStr(*args++);
        const auto resolution = std::stoi(*args++);

        HiCFile f(url);

        if (coord1.end == 0) {
            coord1.end = f.chromosomes().at(coord1.chrom).length;
        }
        if (coord2.end == 0) {
            coord2.end = f.chromosomes().at(coord2.chrom).length;
        }

        auto selector =
            f.getMatrixZoomData(coord1.chrom, coord2.chrom, matrixType, norm, unit, resolution);
        std::vector<contactRecord> buffer{};
        selector.fetch(coord1.start, coord1.end, coord2.start, coord2.end, buffer);

        const auto chrom1 = f.chromosomes().at(coord1.chrom);
        const auto chrom2 = f.chromosomes().at(coord2.chrom);

        for (const auto& record : buffer) {
            const auto start1 = record.binX * resolution;
            const auto end1 =
                std::min(start1 + resolution, static_cast<std::int32_t>(chrom1.length));
            const auto start2 = record.binY * resolution;
            const auto end2 =
                std::min(start2 + resolution, static_cast<std::int32_t>(chrom2.length));
            fmt::print(FMT_COMPILE("{}\t{}\t{}\n"), start1, start2, record.counts);
            // fmt::print(FMT_COMPILE("{:s}\t{:d}\t{:d}\t{:s}\t{:d}\t{:d}\t{:g}\n"), chrom1.name,
            //            start1, end1, chrom2.name, start2, end2, record.counts);
        }
    }

    catch (const std::exception& e) {
        const auto* url = argc == 8 ? *(argv + 3) : *(argv + 2);
        fmt::print(
            stderr,
            FMT_STRING("straw encountered the following error while processing file \"{}\": {}\n"),
            url, e.what());
        return 1;
    }
}
