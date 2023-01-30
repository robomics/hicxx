// Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
//
// SPDX-License-Identifier: MIT

#include <algorithm>
#include <array>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cstdint>
#include <numeric>
#include <string>

#include "straw/straw.hpp"

constexpr auto* pathV8 = "test/data/4DNFIZ1ZVXC8.hic8";
constexpr auto* pathV9 = "test/data/4DNFIZ1ZVXC8.hic9";

static std::vector<contactRecord> head(const std::vector<contactRecord>& buffer,
                                       std::size_t n = 5) {
    REQUIRE(buffer.size() >= n);

    std::vector<contactRecord> slice(n);
    std::copy_n(buffer.begin(), n, slice.begin());
    return slice;
}

static std::vector<contactRecord> tail(const std::vector<contactRecord>& buffer,
                                       std::size_t n = 5) {
    REQUIRE(buffer.size() >= n);

    std::vector<contactRecord> slice(n);
    std::copy_n(buffer.end() - std::int32_t(n), n, slice.begin());
    return slice;
}

template <typename N>
static N sumCounts(const std::vector<contactRecord>& buffer) {
    return std::accumulate(buffer.begin(), buffer.end(), N(0),
                           [](N accumulator, const contactRecord& r) {
                               return accumulator + static_cast<N>(r.count);
                           });
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("matrixZoomData accessors") {
    const auto mzd = HiCFile(pathV8).getMatrixZoomData(
        "chr2L", MatrixType::observed, NormalizationMethod::NONE, MatrixUnit::BP, 2500000);

    CHECK(mzd.chrom1().name == "chr2L");
    CHECK(mzd.chrom2().name == "chr2L");
    CHECK(mzd.matrixType() == MatrixType::observed);
    CHECK(mzd.normalizationMethod() == NormalizationMethod::NONE);
    CHECK(mzd.matrixUnit() == MatrixUnit::BP);
    CHECK(mzd.resolution() == 2500000);

    REQUIRE(mzd.chrom1().length == 23513712);
    CHECK(mzd.numBins1() == 10);
    CHECK(mzd.numBins2() == 10);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("matrixData fetch (observed NONE BP 10000)") {
    std::vector<contactRecord> buffer;
    SECTION("v8") {
        SECTION("intra-chromosomal") {
            HiCFile(pathV8)
                .getMatrixZoomData("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 1433133);
            CHECK(sumCounts<std::int32_t>(buffer) == 19968156);

            constexpr std::size_t N = 5;
            constexpr std::array<float, N> head_expected{1745, 2844, 5719, 409, 1815};
            constexpr std::array<float, N> tail_expected{8, 21, 34, 53, 193};

            const auto h = head(buffer, N);
            const auto t = tail(buffer, N);

            for (std::size_t i = 0; i < N; ++i) {
                CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
                CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
            }

            CHECK_THAT(buffer[1227146].count, Catch::Matchers::WithinRel(1234.0F));
        }
        SECTION("inter-chromosomal") {
            HiCFile(pathV8)
                .getMatrixZoomData("chr2L", "chr4", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 56743);
            CHECK(sumCounts<std::int32_t>(buffer) == 70567);

            constexpr std::size_t N = 5;
            constexpr std::array<float, N> head_expected{1, 1, 1, 1, 1};
            constexpr std::array<float, N> tail_expected{1, 1, 1, 1, 1};

            const auto h = head(buffer, N);
            const auto t = tail(buffer, N);

            for (std::size_t i = 0; i < N; ++i) {
                CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
                CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
            }

            CHECK_THAT(buffer[20324].count, Catch::Matchers::WithinRel(12.0F));
        }
    }
    SECTION("v9") {
        SECTION("intra-chromosomal") {
            HiCFile(pathV9)
                .getMatrixZoomData("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 1433133);
            CHECK(sumCounts<std::int32_t>(buffer) == 19968156);

            constexpr std::size_t N = 5;
            constexpr std::array<float, N> head_expected{1745, 2844, 5719, 409, 1815};
            constexpr std::array<float, N> tail_expected{1, 1, 1, 1, 1};

            const auto h = head(buffer, N);
            const auto t = tail(buffer, N);

            for (std::size_t i = 0; i < N; ++i) {
                CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
                CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
            }

            CHECK_THAT(buffer[870834].count, Catch::Matchers::WithinRel(1234.0F));
        }
        SECTION("inter-chromosomal") {
            HiCFile(pathV9)
                .getMatrixZoomData("chr2L", "chr4", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 56743);
            CHECK(sumCounts<std::int32_t>(buffer) == 70567);

            constexpr std::size_t N = 5;
            constexpr std::array<float, N> head_expected{1, 1, 1, 1, 1};
            constexpr std::array<float, N> tail_expected{1, 1, 1, 1, 1};

            const auto h = head(buffer, N);
            const auto t = tail(buffer, N);

            for (std::size_t i = 0; i < N; ++i) {
                CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
                CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
            }

            CHECK_THAT(buffer[20295].count, Catch::Matchers::WithinRel(12.0F));
        }

        SECTION("inter-chromosomal") {
            HiCFile(pathV9)
                .getMatrixZoomData("chr2L", "chr4", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 56743);
            CHECK(sumCounts<std::int32_t>(buffer) == 70567);

            constexpr std::size_t N = 5;
            constexpr std::array<float, N> head_expected{1, 1, 1, 1, 1};
            constexpr std::array<float, N> tail_expected{1, 1, 1, 1, 1};

            const auto h = head(buffer, N);
            const auto t = tail(buffer, N);

            for (std::size_t i = 0; i < N; ++i) {
                CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
                CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
            }

            CHECK_THAT(buffer[20295].count, Catch::Matchers::WithinRel(12.0F));
        }
    }

    SECTION("sub-queries") {
        SECTION("single pixel") {
            HiCFile(pathV9)
                .getMatrixZoomData("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch("100000-100001", "100000-100001", buffer);
            REQUIRE(buffer.size() == 1);
            CHECK(buffer.front().count == 13895.0F);
        }

        SECTION("upper-triangle") {
            HiCFile(pathV9)
                .getMatrixZoomData("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(123456, 200001, 0, 200001, buffer);
            CHECK(buffer.size() == 140);
            CHECK(sumCounts<std::int32_t>(buffer) == 120811);
        }

        SECTION("lower-triangle") {
            HiCFile(pathV9)
                .getMatrixZoomData("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(0, 200001, 123456, 200001, buffer);
            CHECK(buffer.size() == 140);
            CHECK(sumCounts<std::int32_t>(buffer) == 120811);
        }

        SECTION("inter-chromosomal") {
            HiCFile(pathV9)
                .getMatrixZoomData("chr2L", "chr4", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(123456, 200001, 0, 200001, buffer);
            CHECK(buffer.size() == 55);
            CHECK(sumCounts<std::int32_t>(buffer) == 74);
        }
    }

    SECTION("invalid") {
        HiCFile hic(pathV9);

        SECTION("invalid chromosome") {
            CHECK_THROWS(hic.getMatrixZoomData("chr123", MatrixType::observed,
                                               NormalizationMethod::NONE, MatrixUnit::BP, 10000));
            CHECK_THROWS(hic.getMatrixZoomData(999, MatrixType::observed, NormalizationMethod::NONE,
                                               MatrixUnit::BP, 10000));
        }
        SECTION("invalid resolution") {
            CHECK_THROWS(hic.getMatrixZoomData("chr2L", MatrixType::observed,
                                               NormalizationMethod::NONE, MatrixUnit::BP, -1));
        }
        SECTION("invalid unit") {
            CHECK_THROWS(hic.getMatrixZoomData("chr2L", MatrixType::observed,
                                               NormalizationMethod::NONE, MatrixUnit::FRAG, 10000));
        }
        SECTION("expected + norm") {
            CHECK_THROWS(hic.getMatrixZoomData("chr2L", MatrixType::expected,
                                               NormalizationMethod::VC, MatrixUnit::BP, 10000));
        }
        SECTION("invalid range") {
            CHECK_THROWS(hic.getMatrixZoomData("chr2L", MatrixType::observed,
                                               NormalizationMethod::NONE, MatrixUnit::BP, 10000)
                             .fetch(1000, 0, buffer));
            CHECK_THROWS(hic.getMatrixZoomData("chr2L", MatrixType::observed,
                                               NormalizationMethod::NONE, MatrixUnit::BP, 10000)
                             .fetch(0, 1'000'000'000, buffer));
        }
    }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("matrixData fetch (observed VC BP 10000)") {
    std::vector<contactRecord> buffer;
    SECTION("v8") {
        SECTION("intra-chromosomal") {
            HiCFile(pathV8)
                .getMatrixZoomData("chr2L", MatrixType::observed, NormalizationMethod::VC,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 1433133);
            CHECK_THAT(sumCounts<double>(buffer),
                       Catch::Matchers::WithinRel(20391277.41514, 1.0e-6));
        }
        SECTION("inter-chromosomal") {
            HiCFile(pathV8)
                .getMatrixZoomData("chr2L", "chr4", MatrixType::observed, NormalizationMethod::VC,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 56743);
            CHECK_THAT(sumCounts<double>(buffer),
                       Catch::Matchers::WithinRel(96690.056244753, 1.0e-6));
        }
    }
    SECTION("v9") {
        SECTION("intra-chromosomal") {
            HiCFile(pathV9)
                .getMatrixZoomData("chr2L", MatrixType::observed, NormalizationMethod::VC,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 1433133);
            CHECK_THAT(sumCounts<double>(buffer),
                       Catch::Matchers::WithinRel(20391277.41514, 1.0e-6));
        }
        SECTION("inter-chromosomal") {
            HiCFile(pathV9)
                .getMatrixZoomData("chr2L", "chr4", MatrixType::observed, NormalizationMethod::VC,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 56743);
            CHECK_THAT(sumCounts<double>(buffer),
                       Catch::Matchers::WithinRel(96690.056244753, 1.0e-6));
        }
    }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("matrixData fetch (expected NONE BP 10000)") {
    std::vector<contactRecord> buffer;
    SECTION("v8") {
        SECTION("intra-chromosomal") {
            HiCFile(pathV8)
                .getMatrixZoomData("chr2L", MatrixType::expected, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 1433133);
            CHECK_THAT(sumCounts<double>(buffer),
                       Catch::Matchers::WithinRel(18314748.068024, 1.0e-6));
        }
        SECTION("inter-chromosomal") {
            HiCFile(pathV8)
                .getMatrixZoomData("chr2L", "chr4", MatrixType::expected, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 56743);
            CHECK_THAT(sumCounts<double>(buffer),
                       Catch::Matchers::WithinRel(12610.80619812, 1.0e-6));
        }
    }
    SECTION("v9") {
        SECTION("intra-chromosomal") {
            HiCFile(pathV9)
                .getMatrixZoomData("chr2L", MatrixType::expected, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 1433133);
            CHECK_THAT(sumCounts<double>(buffer),
                       Catch::Matchers::WithinRel(18314748.068024, 1.0e-6));
        }
        SECTION("inter-chromosomal") {
            HiCFile(pathV9)
                .getMatrixZoomData("chr2L", "chr4", MatrixType::expected, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 56743);
            CHECK_THAT(sumCounts<double>(buffer),
                       Catch::Matchers::WithinRel(12610.80619812, 1.0e-6));
        }
    }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("matrixData fetch (oe NONE BP 10000)") {
    std::vector<contactRecord> buffer;
    SECTION("v8") {
        SECTION("intra-chromosomal") {
            HiCFile(pathV8)
                .getMatrixZoomData("chr2L", MatrixType::oe, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 1433133);
            CHECK_THAT(sumCounts<double>(buffer),
                       Catch::Matchers::WithinRel(2785506.2274201, 1.0e-6));
        }
        SECTION("inter-chromosomal") {
            HiCFile(pathV8)
                .getMatrixZoomData("chr2L", "chr4", MatrixType::oe, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 56743);
            CHECK_THAT(sumCounts<double>(buffer),
                       Catch::Matchers::WithinRel(317520.00459671, 1.0e-6));
        }
    }
    SECTION("v9") {
        SECTION("intra-chromosomal") {
            HiCFile(pathV9)
                .getMatrixZoomData("chr2L", MatrixType::oe, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 1433133);
            CHECK_THAT(sumCounts<double>(buffer),
                       Catch::Matchers::WithinRel(2785506.2274201, 1.0e-6));
        }
        SECTION("inter-chromosomal") {
            HiCFile(pathV9)
                .getMatrixZoomData("chr2L", "chr4", MatrixType::oe, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer);
            CHECK(buffer.size() == 56743);
            CHECK_THAT(sumCounts<double>(buffer),
                       Catch::Matchers::WithinRel(317520.00459671, 1.0e-6));
        }
    }
}
