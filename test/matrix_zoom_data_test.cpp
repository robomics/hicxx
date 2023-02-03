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

#include "hicxx/hicxx.hpp"

using namespace hicxx;

constexpr auto* pathV8 = "test/data/4DNFIZ1ZVXC8.hic8";
constexpr auto* pathV9 = "test/data/4DNFIZ1ZVXC8.hic9";

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
static std::vector<contactRecord> head(const std::vector<contactRecord>& buffer,
                                       std::size_t n = 5) {
    REQUIRE(buffer.size() >= n);

    std::vector<contactRecord> slice(n);
    std::copy_n(buffer.begin(), n, slice.begin());
    return slice;
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
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
static void checkContactRecordsAreWithinBound(std::int32_t start1, std::int32_t end1,
                                              std::int32_t start2, std::int32_t end2,
                                              const std::vector<contactRecord>& buffer) {
    assert(start1 < end1);
    assert(start2 < end2);

    for (const auto& r : buffer) {
        CHECK(r.bin1_start >= std::min(start1, start2));
        CHECK(r.bin1_start < std::max(end1, end2));
        CHECK(r.bin2_start >= std::min(start1, start2));
        CHECK(r.bin2_start < std::max(end1, end2));
    }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
static void compareContactRecord(const contactRecord& r1, const contactRecord& r2) {
    CHECK(r1.bin1_start == r2.bin1_start);
    CHECK(r1.bin2_start == r2.bin2_start);
    CHECK_THAT(r1.count, Catch::Matchers::WithinRel(r2.count));
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MatrixSelector accessors") {
    const auto sel = HiCFile(pathV8).getMatrixSelector(
        "chr2L", MatrixType::observed, NormalizationMethod::NONE, MatrixUnit::BP, 2500000);

    CHECK(sel.chrom1().name == "chr2L");
    CHECK(sel.chrom2().name == "chr2L");
    CHECK(sel.matrixType() == MatrixType::observed);
    CHECK(sel.normalizationMethod() == NormalizationMethod::NONE);
    CHECK(sel.matrixUnit() == MatrixUnit::BP);
    CHECK(sel.resolution() == 2500000);

    REQUIRE(sel.chrom1().length == 23513712);
    CHECK(sel.numBins1() == 10);
    CHECK(sel.numBins2() == 10);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MatrixSelector LRU cache") {
    std::vector<contactRecord> buffer;
    HiCFile f(pathV8);

    auto sel = f.getMatrixSelector("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000);

    CHECK(sel.blockCacheHitRate() == 0.0);
    CHECK(sel.blockCacheSize() == 0);

    // Fill cache
    sel.fetch(buffer);
    CHECK(sel.blockCacheHitRate() == 0.0);

    sel.fetch(buffer);
    CHECK(sel.blockCacheHitRate() == 0.5);
    CHECK(sel.blockCacheSize() == 6);

    for (auto i = 0; i < 5; ++i) {
        sel.fetch(buffer);
    }
    CHECK(sel.blockCacheHitRate() == 6.0 / 7.0);
    CHECK(sel.blockCacheSize() == 6);

    sel.clearBlockCache();
    CHECK(sel.blockCacheHitRate() == 0);
    CHECK(sel.blockCacheSize() == 0);
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MatrixSelector fetch (observed NONE BP 10000)") {
    std::vector<contactRecord> buffer;
    SECTION("intra-chromosomal") {
        constexpr std::size_t expected_size = 1433133;
        constexpr std::int32_t expected_sum = 19968156;

        constexpr std::size_t N = 5;
        constexpr std::array<float, N> head_expected{1745, 2844, 5719, 409, 1815};
        constexpr std::array<float, N> tail_expected{8, 21, 34, 53, 193};

        constexpr auto expected_value =
            std::make_pair(std::size_t(770433), contactRecord{15770000, 15770000, 1234.0F});

        SECTION("v8") {
            HiCFile(pathV8)
                .getMatrixSelector("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK(sumCounts<std::int32_t>(buffer) == expected_sum);

            const auto h = head(buffer, N);
            const auto t = tail(buffer, N);

            for (std::size_t i = 0; i < N; ++i) {
                CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
                CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
            }

            compareContactRecord(buffer[expected_value.first], expected_value.second);
        }
        SECTION("v9") {
            HiCFile(pathV9)
                .getMatrixSelector("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK(sumCounts<std::int32_t>(buffer) == expected_sum);

            const auto h = head(buffer, N);
            const auto t = tail(buffer, N);

            for (std::size_t i = 0; i < N; ++i) {
                CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
                CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
            }

            compareContactRecord(buffer[expected_value.first], expected_value.second);
        }
    }

    SECTION("inter-chromosomal") {
        constexpr std::size_t expected_size = 56743;
        constexpr std::int32_t expected_sum = 70567;

        constexpr std::size_t N = 5;
        constexpr std::array<float, N> head_expected{1, 1, 1, 1, 1};
        constexpr std::array<float, N> tail_expected{1, 1, 1, 1, 1};

        constexpr auto expected_value =
            std::make_pair(std::size_t(54023), contactRecord{21770000, 1250000, 12.0F});

        SECTION("v8") {
            HiCFile(pathV8)
                .getMatrixSelector("chr2L", "chr4", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK(sumCounts<std::int32_t>(buffer) == expected_sum);

            const auto h = head(buffer, N);
            const auto t = tail(buffer, N);

            for (std::size_t i = 0; i < N; ++i) {
                CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
                CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
            }

            compareContactRecord(buffer[expected_value.first], expected_value.second);
        }

        SECTION("v9") {
            HiCFile(pathV9)
                .getMatrixSelector("chr2L", "chr4", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK(sumCounts<std::int32_t>(buffer) == expected_sum);

            const auto h = head(buffer, N);
            const auto t = tail(buffer, N);

            for (std::size_t i = 0; i < N; ++i) {
                CHECK_THAT(head_expected[i], Catch::Matchers::WithinRel(h[i].count));
                CHECK_THAT(tail_expected[i], Catch::Matchers::WithinRel(t[i].count));
            }

            compareContactRecord(buffer[expected_value.first], expected_value.second);
        }

        SECTION("cover type 2 interactions") {
            HiCFile(pathV8)
                .getMatrixSelector("chr2L", "chr2R", MatrixType::observed, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 2500000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == 110);
            CHECK(sumCounts<std::int32_t>(buffer) == 1483112);

            compareContactRecord(buffer[53], contactRecord{7500000, 12500000, 16512});
        }

        SECTION("sub-queries") {
            SECTION("single pixel") {
                HiCFile(pathV9)
                    .getMatrixSelector("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                       MatrixUnit::BP, 10000)
                    .fetch("100000-100001", "100000-100001", buffer);
                REQUIRE(buffer.size() == 1);
                compareContactRecord(buffer.front(), contactRecord{100000, 100000, 13895.0F});
            }

            SECTION("upper-triangle") {
                HiCFile(pathV9)
                    .getMatrixSelector("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                       MatrixUnit::BP, 10000)
                    .fetch(123456, 200000, 0, 200000, buffer);
                REQUIRE(buffer.size() == 132);
                CHECK(sumCounts<std::int32_t>(buffer) == 124561);
                compareContactRecord(buffer[17], contactRecord{40000, 130000, 148});
                checkContactRecordsAreWithinBound(123456, 200000, 0, 200000, buffer);
            }

            SECTION("lower-triangle") {
                HiCFile(pathV9)
                    .getMatrixSelector("chr2L", MatrixType::observed, NormalizationMethod::NONE,
                                       MatrixUnit::BP, 10000)
                    .fetch(0, 200000, 123456, 200000, buffer);
                REQUIRE(buffer.size() == 132);
                CHECK(sumCounts<std::int32_t>(buffer) == 124561);
                compareContactRecord(buffer[17], contactRecord{40000, 130000, 148});
                checkContactRecordsAreWithinBound(0, 200000, 123456, 200000, buffer);
            }

            SECTION("inter-chromosomal") {
                HiCFile(pathV9)
                    .getMatrixSelector("chr2L", "chr4", MatrixType::observed,
                                       NormalizationMethod::NONE, MatrixUnit::BP, 10000)
                    .fetch(123456, 200000, 0, 200000, buffer);
                REQUIRE(buffer.size() == 57);
                CHECK(sumCounts<std::int32_t>(buffer) == 74);
                checkContactRecordsAreWithinBound(123456, 200000, 0, 200000, buffer);
            }
        }

        SECTION("invalid") {
            HiCFile hic(pathV9);

            SECTION("invalid chromosome") {
                CHECK_THROWS(hic.getMatrixSelector("chr123", MatrixType::observed,
                                                   NormalizationMethod::NONE, MatrixUnit::BP,
                                                   10000));
                CHECK_THROWS(hic.getMatrixSelector(
                    999, MatrixType::observed, NormalizationMethod::NONE, MatrixUnit::BP, 10000));
            }
            SECTION("invalid resolution") {
                CHECK_THROWS(hic.getMatrixSelector("chr2L", MatrixType::observed,
                                                   NormalizationMethod::NONE, MatrixUnit::BP, -1));
            }
            SECTION("invalid unit") {
                CHECK_THROWS(hic.getMatrixSelector("chr2L", MatrixType::observed,
                                                   NormalizationMethod::NONE, MatrixUnit::FRAG,
                                                   10000));
            }
            SECTION("expected + norm") {
                CHECK_THROWS(hic.getMatrixSelector("chr2L", MatrixType::expected,
                                                   NormalizationMethod::VC, MatrixUnit::BP, 10000));
            }
            SECTION("invalid range") {
                CHECK_THROWS(hic.getMatrixSelector("chr2L", MatrixType::observed,
                                                   NormalizationMethod::NONE, MatrixUnit::BP, 10000)
                                 .fetch(1000, 0, buffer));
                CHECK_THROWS(hic.getMatrixSelector("chr2L", MatrixType::observed,
                                                   NormalizationMethod::NONE, MatrixUnit::BP, 10000)
                                 .fetch(0, 1'000'000'000, buffer));
            }
        }
    }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MatrixSelector fetch (observed VC BP 10000)") {
    std::vector<contactRecord> buffer;
    SECTION("intra-chromosomal") {
        constexpr std::size_t expected_size = 1433133;
        constexpr double expected_sum = 20391277.41514;
        SECTION("v8") {
            HiCFile(pathV8)
                .getMatrixSelector("chr2L", MatrixType::observed, NormalizationMethod::VC,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
        }
        SECTION("v9") {
            HiCFile(pathV9)
                .getMatrixSelector("chr2L", MatrixType::observed, NormalizationMethod::VC,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
        }
    }
    SECTION("inter-chromosomal") {
        constexpr std::size_t expected_size = 56743;
        constexpr double expected_sum = 96690.056244753;
        SECTION("v8") {
            HiCFile(pathV8)
                .getMatrixSelector("chr2L", "chr4", MatrixType::observed, NormalizationMethod::VC,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
        }

        SECTION("v9") {
            HiCFile(pathV9)
                .getMatrixSelector("chr2L", "chr4", MatrixType::observed, NormalizationMethod::VC,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
        }
    }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MatrixSelector fetch (expected NONE BP 10000)") {
    std::vector<contactRecord> buffer;
    SECTION("intra-chromosomal") {
        constexpr std::size_t expected_size = 1433133;
        constexpr double expected_sum = 18314748.068024;
        SECTION("v8") {
            HiCFile(pathV8)
                .getMatrixSelector("chr2L", MatrixType::expected, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
        }
        SECTION("v9") {
            HiCFile(pathV9)
                .getMatrixSelector("chr2L", MatrixType::expected, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
        }
    }
    SECTION("inter-chromosomal") {
        constexpr std::size_t expected_size = 56743;
        constexpr double expected_sum = 12610.80619812;
        SECTION("v8") {
            HiCFile(pathV8)
                .getMatrixSelector("chr2L", "chr4", MatrixType::expected, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
        }

        SECTION("v9") {
            HiCFile(pathV9)
                .getMatrixSelector("chr2L", "chr4", MatrixType::expected, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
        }
    }
}

// NOLINTNEXTLINE(readability-function-cognitive-complexity)
TEST_CASE("MatrixSelector fetch (oe NONE BP 10000)") {
    std::vector<contactRecord> buffer;
    SECTION("intra-chromosomal") {
        constexpr std::size_t expected_size = 1433133;
        constexpr double expected_sum = 2785506.2274201;
        SECTION("v8") {
            HiCFile(pathV8)
                .getMatrixSelector("chr2L", MatrixType::oe, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
        }
        SECTION("v9") {
            HiCFile(pathV9)
                .getMatrixSelector("chr2L", MatrixType::oe, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
        }
    }
    SECTION("inter-chromosomal") {
        constexpr std::size_t expected_size = 56743;
        constexpr double expected_sum = 317520.00459671;
        SECTION("v8") {
            HiCFile(pathV8)
                .getMatrixSelector("chr2L", "chr4", MatrixType::oe, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
        }

        SECTION("v9") {
            HiCFile(pathV9)
                .getMatrixSelector("chr2L", "chr4", MatrixType::oe, NormalizationMethod::NONE,
                                   MatrixUnit::BP, 10000)
                .fetch(buffer, true);
            REQUIRE(buffer.size() == expected_size);
            CHECK_THAT(sumCounts<double>(buffer), Catch::Matchers::WithinRel(expected_sum, 1.0e-6));
        }
    }
}
