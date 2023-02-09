# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT

# cmake-format: off
file(
  DOWNLOAD https://www.dropbox.com/s/zt62d0d3fhbkha0/4DNFIZ1ZVXC8.hic8?dl=1
  EXPECTED_HASH SHA256=4a42d793d0fe0dcbe9dff460f42dbd2e2d2ab16a22a2cc3271b6749dc3e6d8f6
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/4DNFIZ1ZVXC8.hic8)

file(
  DOWNLOAD https://www.dropbox.com/s/tgjy9mzc1u0opps/4DNFIZ1ZVXC8.hic9?dl=1
  EXPECTED_HASH SHA256=29ab3cf0c1bf14a500dd0ee6e2d1b1065c13ffb25e96e95ec4e880d8b4546643
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/4DNFIZ1ZVXC8.hic9)

file(
  DOWNLOAD https://www.dropbox.com/s/lza04jiynt68ouc/data.txt?dl=1
  EXPECTED_HASH SHA256=21c1cdce6ab0e5509b04d84a28000836c7a087cf786efe6f04877ebfff47232a
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/data.txt)

file(
  DOWNLOAD https://www.dropbox.com/s/1tzg3edxpzdmp3y/data.zip?dl=1
  EXPECTED_HASH SHA256=2c4dda145ee44a71e0b85ac2ac243851d639407a3c0f41619fee766a3006b5fc
  ${CMAKE_CURRENT_SOURCE_DIR}/test/data/data.zip)
# cmake-format: on
