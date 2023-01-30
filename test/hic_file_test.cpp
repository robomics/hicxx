// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <catch2/catch_test_macros.hpp>
#include <cstdint>
#include <string>

#include "straw/straw.hpp"

constexpr auto* pathV8 = "test/data/4DNFIZ1ZVXC8.hic8";
#ifdef STRAW_USE_CURL
constexpr auto* urlV8 = "https://www.dropbox.com/s/zt62d0d3fhbkha0/4DNFIZ1ZVXC8.hic8?dl=1";
#endif

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiCFile accessors") {
    SECTION("local") {
        HiCFile f(pathV8);

        CHECK(f.url() == pathV8);
        CHECK(f.name() == pathV8);
        CHECK(f.version() == 8);
        CHECK(f.chromosomes().size() == 9);
        CHECK(f.genomeID() == "dm6");

        CHECK(f.resolutions().size() == 10);
        CHECK(f.resolutions().front() == 2500000);
        CHECK(f.resolutions().back() == 1000);
    }
#ifdef STRAW_USE_CURL
    SECTION("remote") {
        HiCFile f(urlV8);

        CHECK(f.url() == urlV8);
        CHECK(f.name() == urlV8);
        CHECK(f.version() == 8);
        CHECK(f.chromosomes().size() == 9);
        CHECK(f.genomeID() == "dm6");

        CHECK(f.resolutions().size() == 10);
        CHECK(f.resolutions().front() == 2500000);
        CHECK(f.resolutions().back() == 1000);
    }
#endif

    SECTION("invalid") {
        CHECK_THROWS(HiCFile("non-existing-file"));
        CHECK_THROWS(HiCFile("https://localhost:non-existing-url"));
        CHECK_THROWS(HiCFile("test/CMakeLists.txt"));
    }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiCFile footer cache") {
    HiCFile f(pathV8);

    REQUIRE(f.resolutions().size() == 10);

    CHECK(f.numCachedFooters() == 0);
    for (const auto res : f.resolutions()) {
        std::ignore = f.getMatrixZoomData("chr2L", MatrixType::observed, NormalizationMethod::NONE,
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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("HiCFile getMatrixZoomData") {
    HiCFile f(pathV8);

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

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
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
