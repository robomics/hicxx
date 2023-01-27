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
TEST_CASE("readHeader (v8)", "[v8]") {
    constexpr std::array<std::int32_t, 10> resolutions{2500000, 1000000, 500000, 250000, 100000,
                                                       50000,   25000,   10000,  5000,   1000};
    constexpr auto* genomeID = "dm6";
    constexpr auto nChromosomes = 9;

    const auto header = internal::HiCFileStream(urlv8).header();
    CHECK(header.url == urlv8);
    CHECK(header.masterIndexOffset == 131515430);
    CHECK(header.genomeID == genomeID);
    CHECK(header.nChromosomes() == nChromosomes);
    CHECK(header.version == 8);
    CHECK(header.nviPosition == -1);
    CHECK(header.nviLength == -1);

    REQUIRE(header.nResolutions() == resolutions.size());
    for (std::size_t i = 0; i < resolutions.size(); ++i) {
        CHECK(resolutions[i] == header.resolutions[i]);
    }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("readHeader (v9)", "[v9]") {
    constexpr std::array<std::int32_t, 10> resolutions{2500000, 1000000, 500000, 250000, 100000,
                                                       50000,   25000,   10000,  5000,   1000};
    constexpr auto* genomeID = "dm6";
    constexpr auto nChromosomes = 9;

    const auto header = internal::HiCFileStream(urlv9).header();

    CHECK(header.url == urlv9);
    CHECK(header.masterIndexOffset == 130706734);
    CHECK(header.genomeID == genomeID);
    CHECK(header.nChromosomes() == nChromosomes);
    CHECK(header.version == 9);
    CHECK(header.nviPosition == 131417220);
    CHECK(header.nviLength == 6600);

    REQUIRE(header.nResolutions() == resolutions.size());
    for (std::size_t i = 0; i < resolutions.size(); ++i) {
        CHECK(resolutions[i] == header.resolutions[i]);
    }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("readFooter (v8)", "[v8]") {
    internal::HiCFileStream s(urlv8);
    // first 5 expected values
    constexpr std::array<double, 5> expected1{864.6735714977542, 620.9907283534235,
                                              311.1254999778368, 203.9822974509631,
                                              147.9273228359822};
    // last 5 expected values
    constexpr std::array<double, 5> expected2{0.008417076032024847, 0.008417076032024847,
                                              0.008417076032024847, 0.008417076032024847,
                                              0.008417076032024847};

    SECTION("observed NONE BP 5000") {
        const auto f =
            s.readFooter("chr2L", MatrixType::observed, Normalization::NONE, Unit::BP, 5000);

        CHECK(f.matrixType == MatrixType::observed);
        CHECK(f.normalization == Normalization::NONE);
        CHECK(f.unit == Unit::BP);
        CHECK(f.resolution == 5000);
        CHECK(f.fileOffset == 340697);
        CHECK(f.c1NormEntry.position == -1);
        CHECK(f.c1NormEntry.size == -1);
        CHECK(f.c2NormEntry.position == -1);
        CHECK(f.c2NormEntry.size == -1);
        CHECK(f.expectedValues.empty());
    }

    SECTION("observed VC BP 5000") {
        const auto f =
            s.readFooter("chr2L", "chr2R", MatrixType::observed, Normalization::VC, Unit::BP, 5000);

        CHECK(f.matrixType == MatrixType::observed);
        CHECK(f.normalization == Normalization::VC);
        CHECK(f.unit == Unit::BP);
        CHECK(f.resolution == 5000);
        CHECK(f.fileOffset == 11389664);
        CHECK(f.c1NormEntry.position == 134081028);
        CHECK(f.c1NormEntry.size == 37644);
        CHECK(f.c2NormEntry.position == 134231604);
        CHECK(f.c2NormEntry.size == 40516);
        CHECK(f.expectedValues.empty());
    }

    SECTION("observed VC_SQRT BP 5000") {
        const auto f = s.readFooter("chr2L", "chr2R", MatrixType::observed, Normalization::VC_SQRT,
                                    Unit::BP, 5000);

        CHECK(f.matrixType == MatrixType::observed);
        CHECK(f.normalization == Normalization::VC_SQRT);
        CHECK(f.unit == Unit::BP);
        CHECK(f.resolution == 5000);
        CHECK(f.fileOffset == 11389664);
        CHECK(f.c1NormEntry.position == 134118672);
        CHECK(f.c1NormEntry.size == 37644);
        CHECK(f.c2NormEntry.position == 134272120);
        CHECK(f.c2NormEntry.size == 40516);
        CHECK(f.expectedValues.empty());
    }

    SECTION("observed KR BP 5000") {
        const auto f =
            s.readFooter("chr2L", "chr2R", MatrixType::observed, Normalization::KR, Unit::BP, 5000);

        CHECK(f.matrixType == MatrixType::observed);
        CHECK(f.normalization == Normalization::KR);
        CHECK(f.unit == Unit::BP);
        CHECK(f.resolution == 5000);
        CHECK(f.fileOffset == 11389664);
        CHECK(f.c1NormEntry.position == 134156316);
        CHECK(f.c1NormEntry.size == 37644);
        CHECK(f.c2NormEntry.position == 134312636);
        CHECK(f.c2NormEntry.size == 40516);
        CHECK(f.expectedValues.empty());
    }

    SECTION("observed SCALE BP 5000") {
        const auto f = s.readFooter("chr2L", "chr2R", MatrixType::observed, Normalization::SCALE,
                                    Unit::BP, 5000);

        CHECK(f.matrixType == MatrixType::observed);
        CHECK(f.normalization == Normalization::SCALE);
        CHECK(f.unit == Unit::BP);
        CHECK(f.resolution == 5000);
        CHECK(f.fileOffset == 11389664);
        CHECK(f.c1NormEntry.position == 134193960);
        CHECK(f.c1NormEntry.size == 37644);
        CHECK(f.c2NormEntry.position == 134353152);
        CHECK(f.c2NormEntry.size == 40516);
        CHECK(f.expectedValues.empty());
    }

    SECTION("oe NONE BP 5000") {
        const auto f = s.readFooter("chr2L", MatrixType::oe, Normalization::NONE, Unit::BP, 5000);

        CHECK(f.matrixType == MatrixType::oe);
        CHECK(f.normalization == Normalization::NONE);
        CHECK(f.unit == Unit::BP);
        CHECK(f.resolution == 5000);
        CHECK(f.fileOffset == 340697);
        CHECK(f.c1NormEntry.position == -1);
        CHECK(f.c1NormEntry.size == -1);
        CHECK(f.c2NormEntry.position == -1);
        CHECK(f.c2NormEntry.size == -1);
        REQUIRE(f.expectedValues.size() == 6415);

        for (std::size_t i = 0; i < expected1.size(); ++i) {
            const auto j = f.expectedValues.size() - (expected2.size() - i);
            CHECK_THAT(expected1[i], Catch::Matchers::WithinRel(f.expectedValues[i]));
            CHECK_THAT(expected2[i], Catch::Matchers::WithinRel(f.expectedValues[j]));
        }
    }

    SECTION("expected NONE BP 5000") {
        const auto f =
            s.readFooter("chr2L", MatrixType::expected, Normalization::NONE, Unit::BP, 5000);

        CHECK(f.matrixType == MatrixType::expected);
        CHECK(f.normalization == Normalization::NONE);
        CHECK(f.unit == Unit::BP);
        CHECK(f.resolution == 5000);
        CHECK(f.fileOffset == 340697);
        CHECK(f.c1NormEntry.position == -1);
        CHECK(f.c1NormEntry.size == -1);
        CHECK(f.c2NormEntry.position == -1);
        CHECK(f.c2NormEntry.size == -1);
        REQUIRE(f.expectedValues.size() == 6415);

        for (std::size_t i = 0; i < expected1.size(); ++i) {
            const auto j = f.expectedValues.size() - (expected2.size() - i);
            CHECK_THAT(expected1[i], Catch::Matchers::WithinRel(f.expectedValues[i]));
            CHECK_THAT(expected2[i], Catch::Matchers::WithinRel(f.expectedValues[j]));
        }
    }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("readFooter (v9)", "[v9]") {
    internal::HiCFileStream s(urlv9);
    // first 5 expected values
    constexpr std::array<double, 5> expected1{864.6735708339686, 620.990715491172,
                                              311.1255023627755, 203.9822882714327,
                                              147.9273192507429};
    // last 5 expected values
    constexpr std::array<double, 5> expected2{0.008417075820557469, 0.008417075820557469,
                                              0.008417075820557469, 0.008417075820557469,
                                              0.008417075820557469};

    SECTION("observed NONE BP 5000") {
        const auto f =
            s.readFooter("chr2L", MatrixType::observed, Normalization::NONE, Unit::BP, 5000);

        CHECK(f.matrixType == MatrixType::observed);
        CHECK(f.normalization == Normalization::NONE);
        CHECK(f.unit == Unit::BP);
        CHECK(f.resolution == 5000);
        CHECK(f.fileOffset == 340696);
        CHECK(f.c1NormEntry.position == -1);
        CHECK(f.c1NormEntry.size == -1);
        CHECK(f.c2NormEntry.position == -1);
        CHECK(f.c2NormEntry.size == -1);
        CHECK(f.expectedValues.empty());
    }

    SECTION("observed VC BP 5000") {
        const auto f =
            s.readFooter("chr2L", "chr2R", MatrixType::observed, Normalization::VC, Unit::BP, 5000);

        CHECK(f.matrixType == MatrixType::observed);
        CHECK(f.normalization == Normalization::VC);
        CHECK(f.unit == Unit::BP);
        CHECK(f.resolution == 5000);
        CHECK(f.fileOffset == 11625116);
        CHECK(f.c1NormEntry.position == 131715708);
        CHECK(f.c1NormEntry.size == 18820);
        CHECK(f.c2NormEntry.position == 131772168);
        CHECK(f.c2NormEntry.size == 20240);
        CHECK(f.expectedValues.empty());
    }

    SECTION("observed VC_SQRT BP 5000") {
        const auto f = s.readFooter("chr2L", "chr2R", MatrixType::observed, Normalization::VC_SQRT,
                                    Unit::BP, 5000);

        CHECK(f.matrixType == MatrixType::observed);
        CHECK(f.normalization == Normalization::VC_SQRT);
        CHECK(f.unit == Unit::BP);
        CHECK(f.resolution == 5000);
        CHECK(f.fileOffset == 11625116);
        CHECK(f.c1NormEntry.position == 131734528);
        CHECK(f.c1NormEntry.size == 18820);
        CHECK(f.c2NormEntry.position == 131792408);
        CHECK(f.c2NormEntry.size == 20240);
        CHECK(f.expectedValues.empty());
    }

    /*  TODO: for some reason KR normalization is missing
    SECTION("observed KR BP 5000") {
        const auto f = s.readFooter("chr2L", "chr2R", MatrixType::observed, Normalization::KR,
                                    Unit::BP, 5000);

        CHECK(f.matrixType == MatrixType::observed);
        CHECK(f.normalization == Normalization::KR);
        CHECK(f.unit == Unit::BP);
        CHECK(f.resolution == 5000);
        CHECK(f.fileOffset == 11625116);
        CHECK(f.c1NormEntry.position == -1);  // TODO
        CHECK(f.c1NormEntry.size == -1);      // TODO
        CHECK(f.c2NormEntry.position == -1);  // TODO
        CHECK(f.c2NormEntry.size == -1);      // TODO
        CHECK(f.expectedValues.empty());
    } */

    SECTION("observed SCALE BP 5000") {
        const auto f = s.readFooter("chr2L", "chr2R", MatrixType::observed, Normalization::SCALE,
                                    Unit::BP, 5000);

        CHECK(f.matrixType == MatrixType::observed);
        CHECK(f.normalization == Normalization::SCALE);
        CHECK(f.unit == Unit::BP);
        CHECK(f.resolution == 5000);
        CHECK(f.fileOffset == 11625116);
        CHECK(f.c1NormEntry.position == 131753348);
        CHECK(f.c1NormEntry.size == 18820);
        CHECK(f.c2NormEntry.position == 131812648);
        CHECK(f.c2NormEntry.size == 20240);
        CHECK(f.expectedValues.empty());
    }

    SECTION("oe NONE BP 5000") {
        const auto f = s.readFooter("chr2L", MatrixType::oe, Normalization::NONE, Unit::BP, 5000);

        CHECK(f.matrixType == MatrixType::oe);
        CHECK(f.normalization == Normalization::NONE);
        CHECK(f.unit == Unit::BP);
        CHECK(f.resolution == 5000);
        CHECK(f.fileOffset == 340696);
        CHECK(f.c1NormEntry.position == -1);
        CHECK(f.c1NormEntry.size == -1);
        CHECK(f.c2NormEntry.position == -1);
        CHECK(f.c2NormEntry.size == -1);
        REQUIRE(f.expectedValues.size() == 6415);

        for (std::size_t i = 0; i < expected1.size(); ++i) {
            const auto j = f.expectedValues.size() - (expected2.size() - i);
            CHECK_THAT(expected1[i], Catch::Matchers::WithinRel(f.expectedValues[i]));
            CHECK_THAT(expected2[i], Catch::Matchers::WithinRel(f.expectedValues[j]));
        }
    }

    SECTION("expected NONE BP 5000") {
        const auto f =
            s.readFooter("chr2L", MatrixType::expected, Normalization::NONE, Unit::BP, 5000);

        CHECK(f.matrixType == MatrixType::expected);
        CHECK(f.normalization == Normalization::NONE);
        CHECK(f.unit == Unit::BP);
        CHECK(f.resolution == 5000);
        CHECK(f.fileOffset == 340696);
        CHECK(f.c1NormEntry.position == -1);
        CHECK(f.c1NormEntry.size == -1);
        CHECK(f.c2NormEntry.position == -1);
        CHECK(f.c2NormEntry.size == -1);
        REQUIRE(f.expectedValues.size() == 6415);

        for (std::size_t i = 0; i < expected1.size(); ++i) {
            const auto j = f.expectedValues.size() - (expected2.size() - i);
            CHECK_THAT(expected1[i], Catch::Matchers::WithinRel(f.expectedValues[i]));
            CHECK_THAT(expected2[i], Catch::Matchers::WithinRel(f.expectedValues[j]));
        }
    }
}
