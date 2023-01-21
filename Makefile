elpaso-intel:
	cd ./project && rm -rf * && cmake .. -DCMAKE_C_COMPILER=icc -DCMAKE_CXX_COMPILER=icpc

elpaso-gnu:
	cd ./project && rm -rf * && cmake ..

elpasolib-build:
	cd ./project && rm -rf ../*cxx*-o && cmake --build . --target install -j40 

elpaso-build:
	cd ./project && make all -j40

cmake-gen:
	cd ./cmake && cp YOURCOMPUTER.config.cmake ${HOSTNAME}.config.cmake

cmake-gen-ci:
	cd ./cmake && cp YOURCOMPUTER.config-ci.cmake ${HOSTNAME}.config.cmake

cmake-gen-lib:
	cd ./cmake && cp YOURCOMPUTER.config-lib.cmake ${HOSTNAME}.config.cmake

cmake-gen-lint:
	cd ./cmake && cp YOURCOMPUTER.config-lint-clang.cmake ${HOSTNAME}.config.cmake