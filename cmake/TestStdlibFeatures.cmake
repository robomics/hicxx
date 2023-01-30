# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

include(CheckCXXSourceCompiles)

set(CMAKE_REQUIRED_FLAGS, "-std=c++17")

#
# Detect <variant> availablilty
#
check_cxx_source_compiles(
    "#include <variant>

   int main(int argc, char* argv[]) {
      std::variant<int, long> x;
      (void)x;
      return 0;
   }"
    VARIANT_AVAILABLE)

unset(CMAKE_REQUIRED_FLAGS)
