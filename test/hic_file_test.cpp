// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdint>
#include <string>

#include "straw/straw.h"

constexpr auto* urlv8 = "test/data/4DNFIZ1ZVXC8.hic8";
constexpr auto* urlv9 = "test/data/4DNFIZ1ZVXC8.hic9";

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiCFile accessors") {
    HiCFile f(urlv8);

    CHECK(f.url() == urlv8);
    CHECK(f.name() == urlv8);
    CHECK(f.version() == 8);
    CHECK(f.chromosomes().size() == 9);
    CHECK(f.genomeID() == "dm6");

    CHECK(f.resolutions().size() == 10);
    CHECK(f.resolutions().front() == 2500000);
    CHECK(f.resolutions().back() == 1000);
}

TEST_CASE("HiCFile footer cache") {
    HiCFile f(urlv8);

    REQUIRE(f.resolutions().size() == 10);

    CHECK(f.numCachedFooters() == 0);
    for (const auto res : f.resolutions()) {
        (void)f.getMatrixZoomData("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                  MatrixUnit::BP, res);
    }

    CHECK(f.numCachedFooters() == f.resolutions().size());

    const auto sel1 = f.getMatrixZoomData("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                          MatrixUnit::BP, 2500000);
    const auto sel2 = f.getMatrixZoomData("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                          MatrixUnit::BP, 2500000);

    // this check relies on the fact that chrom1Norm are stored in the footer, and that footers are
    // looked up in the cache when creating matrix selectors
    CHECK(&sel1.chrom1Norm() == &sel2.chrom1Norm());

    f.purgeFooterCache();
    CHECK(f.numCachedFooters() == 0);

    const auto sel3 = f.getMatrixZoomData("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                          MatrixUnit::BP, 2500000);

    CHECK(f.numCachedFooters() == 1);
    CHECK(&sel1.chrom1Norm() != &sel3.chrom1Norm());
}

TEST_CASE("HiCFile getMatrixZoomData") {
    HiCFile f(urlv8);

    REQUIRE(f.chromosomes().size() == 9);

    constexpr auto mt = MatrixType::observed;
    constexpr auto norm = NormalizationMethod::NONE;
    constexpr auto unit = MatrixUnit::BP;
    constexpr std::int32_t res = 2500000;

    const auto chrom1 = f.chromosomes().at("chr2L");
    const auto chrom2 = f.chromosomes().at("chr2R");

    SECTION("intra-chromosomal") {
        auto sel = f.getMatrixZoomData(chrom1.name, mt, norm, unit, res);
        CHECK(sel.chrom1() == chrom1);
        CHECK(sel.isIntra());

        sel = f.getMatrixZoomData(chrom1.index, mt, norm, unit, res);
        CHECK(sel.chrom1() == chrom1);
        CHECK(sel.isIntra());

        sel = f.getMatrixZoomData(chrom1, mt, norm, unit, res);
        CHECK(sel.chrom1() == chrom1);
        CHECK(sel.isIntra());
    }

    SECTION("inter-chromosomal") {
        auto sel = f.getMatrixZoomData(chrom1.name, chrom2.name, mt, norm, unit, res);
        CHECK(sel.chrom1() == chrom1);
        CHECK(sel.chrom2() == chrom2);

        sel = f.getMatrixZoomData(chrom1.index, chrom2.index, mt, norm, unit, res);
        CHECK(sel.chrom1() == chrom1);
        CHECK(sel.chrom2() == chrom2);

        sel = f.getMatrixZoomData(chrom1, chrom2, mt, norm, unit, res);
        CHECK(sel.chrom1() == chrom1);
        CHECK(sel.chrom2() == chrom2);

        sel = f.getMatrixZoomData(chrom2, chrom1, mt, norm, unit, res);
        CHECK(sel.chrom1() == chrom1);
        CHECK(sel.chrom2() == chrom2);
    }

    SECTION("invalid chromosome") {
        CHECK_THROWS(f.getMatrixZoomData("not-a-chromosome", mt, norm, unit, res));
        CHECK_THROWS(f.getMatrixZoomData(chrom1.name, "not-a-chromosome", mt, norm, unit, res));
        CHECK_THROWS(f.getMatrixZoomData(999, mt, norm, unit, res));
        CHECK_THROWS(f.getMatrixZoomData(chrom1.index, 999, mt, norm, unit, res));
    }

    SECTION("malformed") {
        CHECK_THROWS(f.getMatrixZoomData(chrom1, mt, norm, unit, 123));
        CHECK_THROWS(
            f.getMatrixZoomData(chrom1, MatrixType::expected, NormalizationMethod::VC, unit, res));

        // Matrix does not have contacts for fragments
        CHECK_THROWS(f.getMatrixZoomData(chrom1, mt, norm, MatrixUnit::FRAG, res));
    }
}

TEST_CASE("GenomicCoordinates") {
    SECTION("chrom-only") {
        const auto coord = GenomicCoordinates::fromString("chr1");
        CHECK(coord.chrom == "chr1");
        CHECK(coord.start == 0);
        CHECK(coord.end == 0);
    }

    SECTION("valid") {
        auto coord = GenomicCoordinates::fromString("chr1:0-1000");
        CHECK(coord.chrom == "chr1");
        CHECK(coord.start == 0);
        CHECK(coord.end == 1000);

        coord = GenomicCoordinates::fromString("chr1:0:1000");
        CHECK(coord.chrom == "chr1");
        CHECK(coord.start == 0);
        CHECK(coord.end == 1000);
    }

    SECTION("invalid") {
        CHECK_THROWS(GenomicCoordinates::fromString("chr1:0"));
        CHECK_THROWS(GenomicCoordinates::fromString("chr1:0:"));

        CHECK_THROWS(GenomicCoordinates::fromString("chr1:0-"));
        CHECK_THROWS(GenomicCoordinates::fromString("chr1:-"));

        CHECK_THROWS(GenomicCoordinates::fromString("chr1::"));
        CHECK_THROWS(GenomicCoordinates::fromString("chr1:a:b"));
        CHECK_THROWS(GenomicCoordinates::fromString("chr1:a-b"));

        CHECK_THROWS(GenomicCoordinates::fromString("chr1:100-0"));
        CHECK_THROWS(GenomicCoordinates::fromString("chr1:0-100suffix"));
    }
}
