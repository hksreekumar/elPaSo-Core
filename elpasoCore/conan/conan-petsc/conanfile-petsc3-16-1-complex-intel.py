from conans import ConanFile, CMake, tools


class PetscConan(ConanFile):
    name = "petsc-complex"
    version = "3.16.1"
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
        self.copy("*.so", src="intel-cxx-complex-o/lib", dst="lib", keep_path=True)
        self.copy("*.so.*", src="intel-cxx-complex-o/lib", dst="lib", keep_path=True)
        self.copy("*.a", src="intel-cxx-complex-o/lib", dst="lib", keep_path=True)

    def package_info(self):
        self.cpp_info.libs = ["libpetsc.so","libscalapack.a","libparmetis.so","libpetsc.so.3.16","libdmumps.a","libzmumps.a","libcmumps.a","libsmumps.a","libmumps_common.a","libEl.so","libpmrrr.so","libmetis.so","libElSuiteSparse.so"]
        pass

