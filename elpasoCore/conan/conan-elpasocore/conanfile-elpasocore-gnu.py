from conans import ConanFile, CMake, tools


class elpasoCoreConan(ConanFile):
    name = "elpasocore"
    version = "23.05.1"
    license = "All rights reserved"
    author = "Harikrishnan Sreekumar"
    url = "https://www.tu-braunschweig.de/en/ina/institute/ina-tech/research-code-elpaso"
    description = "Prebuilt elpaso-core"
    topics = ("<Put some tag here>", "<here>", "<and here>")
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}

    def config_options(self):
        pass

    def source(self):
        pass

    def build(self):
        pass

    def package(self):
        self.copy("elpaso*", src=".", dst=".", keep_path=True)
        self.copy("*.h", src="include", dst="include", keep_path=True)
        self.copy("**", src="bin", dst="bin", keep_path=True)
        self.copy("*", src="lib", dst="lib", keep_path=True)
        self.copy("*", src="share", dst="share", keep_path=True)

    def package_info(self):
        pass