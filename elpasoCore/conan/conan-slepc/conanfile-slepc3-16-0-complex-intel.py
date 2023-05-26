from conans import ConanFile, CMake, tools


class SlepcConan(ConanFile):
    name = "slepc-complex"
    version = "3.16.0"
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
        self.copy("*.h", src="intel-cxx-complex-o/include", dst="include", keep_path=True)
        self.copy("*.h", src="include", dst="include", keep_path=True)
        self.copy("*.so", src="intel-cxx-complex-o", dst="lib", keep_path=False)
        self.copy("*.so.*", src="intel-cxx-complex-o", dst="lib", keep_path=False)
        self.copy("*.a", src="intel-cxx-complex-o", dst="lib", keep_path=False)

    def package_info(self):
        self.cpp_info.libs = ["libslepc.so"]
        pass

