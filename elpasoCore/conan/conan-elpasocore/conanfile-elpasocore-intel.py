from conans import ConanFile, CMake, tools


class elpasoCoreConan(ConanFile):
    name = "elpaso-core"
    version = "23.01.1"
    license = "All rights reserved"
    author = "Harikrishnan Sreekumar"
    url = "https://www.tu-braunschweig.de/en/ina/institute/ina-tech/research-code-elpaso"
    description = "Lib assist for elPaSo build"
    topics = ("<Put some tag here>", "<here>", "<and here>")
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}
    generators = "cmake_find_package"

    def config_options(self):
        pass

    def source(self):
        pass

    def build(self):
        pass

    def package(self):
        self.copy("*.h", dst="include", keep_path=True)
        self.copy("*.so", dst="lib", keep_path=True)
        self.copy("*.so.*", dst="lib", keep_path=True)

    def package_info(self):
        self.cpp_info.libs = ["libelpasoCore-intel-cxx-complex-o.so", "libelpasoCore-intel-cxx-o"]
        pass

