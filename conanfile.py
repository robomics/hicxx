# Copyright (C) 2023 Roberto Rossini <roberros@uio.no>
#
# SPDX-License-Identifier: MIT


from conan import ConanFile
from conan.tools.build import check_min_cppstd

required_conan_version = ">=1.53.0"


class Hicxx(ConanFile):
    name = "hicxx"
    homepage = "https://github.com/robomics/hicxx"
    license = "MIT"
    author = "Roberto Rossini (roberros@uio.no)"
    settings = "os", "compiler", "build_type", "arch"
    requires = ["catch2/3.3.1",
                "fmt/9.1.0",
                "tsl-ordered-map/1.1.0",
                "zlib/1.2.13"]

    generators = "cmake", "cmake_find_package", "cmake_find_package_multi"

    def validate(self):
        if self.settings.compiler.get_safe("cppstd"):
            check_min_cppstd(self, 17)

    def configure(self):
        if self.settings.compiler in ["clang", "gcc"]:
            self.settings.compiler.libcxx = "libstdc++11"

        self.options["fmt"].header_only = True

    def imports(self):
        self.copy("license*", dst="licenses", folder=True, ignore_case=True)
