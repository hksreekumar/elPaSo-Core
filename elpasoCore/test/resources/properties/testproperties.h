/* Copyright (c) 2023. Authors listed in AUTHORS.md

 * This file is part of elPaSo-Core.

 * elPaSo-Core is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your option)
 * any later version.

 * elPaSo-Core is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.

 * You should have received a copy of the GNU Lesser General Public License along
 * with elPaSo-Core (COPYING.txt and COPYING.LESSER.txt). If not, see
 * <https://www.gnu.org/licenses/>. 
 */

//! @brief Globally accessible properties for unit testing
//! @author Harikrishnan Sreekumar
//! @date 23.04.2021

// XML file to test frequency analysis
#define XML_FREQUENCYANALYSIS_FILE "plate_SS_ssd.ak3"

// XML file to test mor-offline analysis
#define XML_MORANALYSIS_FILE	   "plate_SS_mor.ak3"

// hdf5 file to test mor-offline analysis
#define HDF5_MORANALYSIS_FILE	   "plate_small_mor.hdf5"

// HDF5 file to test basic operations
#define HDF5_PARSER_TEST_FILE		"TestContainer.hdf5"

// HDF5 file to test elements, nodes, analysis, element nodes, nodal nodes, materials
#define HDF5_PARSER_BASIC_FILE		"plate.hdf5"

// HDF5 file to test interface and fluid elements
#define HDF5_PARSER_INTERFACEFLUIDELEM_FILE		"plate2_interface_and_fluid.hdf5"

// HDF5 file to test nodal constraints
#define HDF5_PARSER_NODECONSTRAINTS_FILE		"plate1_structBC.hdf5"

// HDF5 file to test sample outputs
#define HDF5_PARSER_OUTPUT_EXAMPLE_FILE		"elpasoT_TEMPORARY.hdf5"

// Log file to test sample outputs with OpenMPI
#define LOG_OUTPUT_OMPI_EXAMPLE_FILE		"test_logfile_ompi.log.0"

// Log file to test sample outputs with IntelMPI
#define LOG_OUTPUT_IMPI_EXAMPLE_FILE		"test_logfile_intelmpi.log.0"

// Path to elpasoT executable from ELPASO_TEST_RESOURCE_DIR
#define PATH_TO_ELPASOT_EXECUTABLE		"../../../bin/"

// HDF5 file to test Structural isotropic linear visco-elastic material
#define HDF5_STR_LIN_VIS_ISO_DIR		"plate_dsg4_2elem_STR_LIN_VIS_ISO_DIR.hdf5"

// HDF5 file to test frequency dependent structural isotropic linear visco-elastic material
#define HDF5_FREQ_STR_LIN_VIS_ISO_DIR		"plate_dsg4_2elemf_STR_LIN_VIS_ISO_DIR.hdf5"

// HDF5 file to test Structural isotropic linear elastic material
#define HDF5_STR_LIN_ELA_ISO_DIR		"plate_dsg4_2elem_STR_LIN_ELA_ISO_DIR.hdf5"

// HDF5 file to test structural orthotropic linear visco-elastic material
#define HDF5_STR_LIN_VIS_ORT_DIR			"plate_shell4_2elem_STR_LIN_VIS_ORT_DIR.hdf5"

// HDF5 file to test structural orthotropic linear visco-elastic laminate material
#define HDF5_STR_LIN_VIS_ORT_LAM			"plate_shell4_2elem_STR_LIN_VIS_ORT_LAM.hdf5"

// HDF5 file to test fluid linear viscoelastic material
#define HDF5_AF_LIN_EQF_ISO_DIR			"acousticmedium_fluid8_2elem_AF_LIN_EQF_ISO_DIR.hdf5"

// HDF5 file to test fluid linear elastic material
#define HDF5_AF_LIN_UAF_ISO_DIR			"acousticmedium_fluid8_2elem_AF_LIN_UAF_ISO_DIR.hdf5"

// HDF5 file to test beam -> Structural isotropic linear elastic material
#define HDF5_BEAM_STR_LIN_ELA_ISO_DIR		"beam_4elem.hdf5"

// HDF5 file to test structural isotropic linear elastic spring material
#define HDF5_STR_LIN_SPR_ORT_DIR		"plate_spring_dsg4_2elem_STR_LIN_SPR_ORT_DIR.hdf5"

// HDF5file to test structural boundary condition
#define HDF5_BC_STRUCTURE		"plate_shell4_uwbc.hdf5"

// HDF5 file to test structural orthotropic linear elastic material
#define HDF5_STR_LIN_ELA_ORT_DIR		"plate_shell4_2elem_STR_LIN_ELA_ORT_DIR.hdf5"