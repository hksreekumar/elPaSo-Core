#---------------------------------------------------------------------------#
# elPaSo - Configuration for additional dependencies
#
# 14.10.2020
# Harikrishnan Sreekumar
# Institute for Acoustics, Technische Universitaet Braunschweig
#---------------------------------------------------------------------------#

MACRO(DEPENDENCY_CONFIG)
	# ---- Master control for dynamic library or executable ----------------#
	OPTION(GEN_DLIB			"elPaSo DYNAMIC LIB"	OFF) # switch off for generating exectutable

	# ---- Master control for Lib installation from CONAN repos ------------#
	OPTION(LINK_CONAN		"CONAN LIBS"		ON)

	# ---- Base directory for additional dependencies ----------------------#
	SET(ELPASO_LIB_DIR		"/opt/elPaSo")
	SET(INTEL_DIR			"/software/intel")

	# ---- Set version (MPI is linked automatically) -----------------------#
	SET(PETSC_VERSION		3.19.1)
	SET(SLEPC_VERSION		3.19.1)
	SET(XMLIO_VERSION		0.92)
	SET(HDF5_VERSION		1.12.0)
	SET(OPENMPI_VERSION		3.1.6)
	SET(INTELMPI_VERSION	2019.6.166)

	# ---- Swicht for elPaSo Projects --------------------------------------#
	OPTION(LINK_ELPASO		"elpaso : Real"		ON)
	OPTION(LINK_ELPASOC		"elpasoC: Complex"	ON)
	OPTION(LINK_ELPASOT		"elpasoT: Test"		ON)
	
	# ---- Switch for different libraries ----------------------------------#
	OPTION(LINK_MKL			"MKL LIB"			ON)
	OPTION(LINK_PETSC		"PETSc LIB"			ON)
	OPTION(LINK_SLEPC		"SLEPc LIB"			ON)
	OPTION(LINK_ARPACK		"ARPACK LIB"		ON)
	OPTION(LINK_XMLIO		"XMLIO LIB"			OFF)
	OPTION(LINK_GZSTREAM	"GZSTREAM LIB"		ON)
	OPTION(LINK_HDF5		"HDF5 LIB"			ON)
	OPTION(LINK_MUMPS  		"MUMPS SOLVER"		ON)
	OPTION(LINK_PARDISO		"PARDISO SOLVER"	ON)
	OPTION(LINK_COV			"COV TOOLS"			OFF)
	OPTION(LINK_GTEST		"GTEST LIB"			ON)
	OPTION(LINK_CLANGTIDY	"CLANG-TIDY"		OFF)
	OPTION(LINK_BOOST		"BOOST LIB"			OFF)
	
	# ---- Switch for different FE module ----------------------------------#
	OPTION(INCLUDE_ANALYSIS_ESSENTIALS		"Analysis Essential"	ON)
	OPTION(INCLUDE_BC_ESSENTIALS			"BC Essential"			ON)
	OPTION(INCLUDE_ELEMENT_ESSENTIALS		"Elements Essential"	ON)
	OPTION(INCLUDE_MATERIAL_ESSENTIALS		"Materials Essentials"	ON)
	OPTION(INCLUDE_HDF5						"HDF5"					ON)

	# ---- Set LIBS to be linked by Cmake ----------------------------------#
	SET(PETSC_ARCH_REAL ${ELPASO_COMPILER_ID}-cxx-o)
	SET(PETSC_ARCH_COMPLEX ${ELPASO_COMPILER_ID}-cxx-complex-o)

ENDMACRO(DEPENDENCY_CONFIG)
