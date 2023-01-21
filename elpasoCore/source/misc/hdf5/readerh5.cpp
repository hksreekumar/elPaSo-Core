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

#include "inputh5.h"
#include <string>
#include <iostream>

#define RANK_VEC   1

/*BEGIN_NO_COVERAGE*/
cReaderH5::cReaderH5()
{
#ifdef HAVE_HDF5
    H5Eset_auto(NULL, NULL, NULL);
#endif // HAVE_HDF5
}

cReaderH5::~cReaderH5()
{
    // empty
}
/*END_NO_COVERAGE*/

#ifdef PETSC_USE_COMPLEX
void cReaderH5::readDenseComplexMatrix(std::vector<PetscComplex>& _vector, int& _m, int& _n, std::string _nameFirstLevel, std::string _nameSecondLevel, std::string _nameThirdLevel)
{
#ifdef HAVE_HDF5
    std::string datasetname = _nameFirstLevel + _nameSecondLevel + _nameThirdLevel;
    myDataSet = H5Dopen(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

    myDataSpace = H5Dget_space(myDataSet);
    hsize_t dims_out[2];
    int ndims = H5Sget_simple_extent_dims(myDataSpace, dims_out, NULL);

    _m = (int)dims_out[0];
    _n = (int)dims_out[1];

    std::vector<elpasoComplexDouble> convertPetscComplex;
    convertPetscComplex.resize(dims_out[0] * dims_out[1]);
    _vector.resize(dims_out[0] * dims_out[1]);

    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(elpasoComplexDouble));
    H5Tinsert(memtype, "real", HOFFSET(elpasoComplexDouble, real), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "imag", HOFFSET(elpasoComplexDouble, imag), H5T_NATIVE_DOUBLE);

    H5Dread(myDataSet, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &convertPetscComplex[0]);
    for (size_t i = 0; i < convertPetscComplex.size(); i++)
        _vector[i] = { convertPetscComplex[i].real, convertPetscComplex[i].imag };

    H5Sclose(myDataSpace);
    H5Dclose(myDataSet);
#endif
}
#endif

void cReaderH5::readIntegerVector(std::vector<int>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef HAVE_HDF5
    std::string datasetname = _nameMainGroup + _nameMatGroup + _nameVecGroup;
    myDataSet = H5Dopen(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

    myDataSpace = H5Dget_space(myDataSet);
    hsize_t dims_out[2];
    int ndims = H5Sget_simple_extent_dims(myDataSpace, dims_out, NULL);

    _vector.resize(dims_out[0]);
    H5Dread(myDataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_vector[0]);

    H5Sclose(myDataSpace);
    H5Dclose(myDataSet);
#endif
}

void cReaderH5::readIntegerVectorFromAttributeDataset(std::vector<int>& _vector, std::string _path, std::string _dataset, std::string _attribute)
{
#ifdef HAVE_HDF5
    std::string datasetname = _path + _dataset;
    myDataSet = H5Dopen(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

    hid_t myatt_out = H5Aopen(myDataSet, _attribute.c_str(), H5P_DEFAULT);

    myDataSpace = H5Aget_space(myatt_out);

    hsize_t dims_out[2];
    int ndims = H5Sget_simple_extent_dims(myDataSpace, dims_out, NULL);

    _vector.resize(dims_out[0]);
    H5Aread(myatt_out, H5T_NATIVE_INT, &_vector[0]);

    H5Aclose(myatt_out);
    H5Sclose(myDataSpace);
    H5Dclose(myDataSet);
#endif
}

void cReaderH5::readIntegerVectorFromAttributeGroup(std::vector<int>& _vector, std::string _groupname, std::string _attribute)
{
#ifdef HAVE_HDF5
    try{
        hid_t group_id = H5Gopen(cFileH5::getFileInstance(), _groupname.c_str(), H5P_DEFAULT);

        hid_t myatt_out = H5Aopen(group_id, _attribute.c_str(), H5P_DEFAULT);

        myDataSpace = H5Aget_space(myatt_out);

        hsize_t dims_out[2];
        int ndims = H5Sget_simple_extent_dims(myDataSpace, dims_out, NULL);

        _vector.resize(dims_out[0]);
        H5Aread(myatt_out, H5T_NATIVE_INT, &_vector[0]);

        H5Aclose(myatt_out);
        H5Sclose(myDataSpace);
        H5Gclose(group_id);
    }
    catch(...)
    {
        // do nothing
    }
#endif
}


void cReaderH5::readDenseDoubleMatrix(std::vector<double>& _vector, int& _m, int& _n, std::string _nameFirstLevel, std::string _nameSecondLevel, std::string _nameThirdLevel)
{
#ifdef HAVE_HDF5
    std::string datasetname = _nameFirstLevel + _nameSecondLevel + _nameThirdLevel;
    myDataSet = H5Dopen(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

    myDataSpace = H5Dget_space(myDataSet);
    hsize_t dims_out[2];
    int ndims = H5Sget_simple_extent_dims(myDataSpace, dims_out, NULL);

    _vector.resize(dims_out[0] * dims_out[1]);
    H5Dread(myDataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_vector[0]);

    _m = (int)dims_out[0];
    _n = (int)dims_out[1];

    H5Sclose(myDataSpace);
    H5Dclose(myDataSet);
#endif 
}

void cReaderH5::readDenseIntegerMatrix(std::vector<int>& _vector, int& _m, int& _n, std::string _nameFirstLevel, std::string _nameSecondLevel, std::string _nameThirdLevel)
{
#ifdef HAVE_HDF5
    std::string datasetname = _nameFirstLevel + _nameSecondLevel + _nameThirdLevel;
    myDataSet = H5Dopen(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

    myDataSpace = H5Dget_space(myDataSet);
    hsize_t dims_out[2];
    int ndims = H5Sget_simple_extent_dims(myDataSpace, dims_out, NULL);

    _vector.resize(dims_out[0] * dims_out[1]);
    H5Dread(myDataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_vector[0]);

    _m = (int)dims_out[0];
    _n = (int)dims_out[1];

    H5Sclose(myDataSpace);
    H5Dclose(myDataSet);
#endif 
}

void cReaderH5::readDoubleVector(std::vector<double>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef HAVE_HDF5
    try{
        std::string datasetname = _nameMainGroup + _nameMatGroup + _nameVecGroup;
        myDataSet = H5Dopen(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

        myDataSpace = H5Dget_space(myDataSet);
        hsize_t dims_out[2];
        int ndims = H5Sget_simple_extent_dims(myDataSpace, dims_out, NULL);

        _vector.resize(dims_out[0]);
        H5Dread(myDataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &_vector[0]);

        H5Sclose(myDataSpace);
        H5Dclose(myDataSet);
    }
    catch(...)
    {
        // do nothing
    }
#endif
}

#ifdef PETSC_USE_COMPLEX
void cReaderH5::readComplexVector(std::vector<PetscComplex>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef HAVE_HDF5
    std::string datasetname = _nameMainGroup + _nameMatGroup + _nameVecGroup;
    myDataSet = H5Dopen(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

    myDataSpace = H5Dget_space(myDataSet);
    hsize_t dims_out[2];
    int ndims = H5Sget_simple_extent_dims(myDataSpace, dims_out, NULL);

    std::vector<elpasoComplexDouble> convertPetscComplex;
    convertPetscComplex.resize(dims_out[0]);
    _vector.resize(dims_out[0]);

    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(elpasoComplexDouble));
    H5Tinsert(memtype, "real", HOFFSET(elpasoComplexDouble, real), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "imag", HOFFSET(elpasoComplexDouble, imag), H5T_NATIVE_DOUBLE);

    H5Dread(myDataSet, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &convertPetscComplex[0]);
    for (size_t i = 0; i < convertPetscComplex.size(); i++) {
        _vector[i] = { convertPetscComplex[i].real, convertPetscComplex[i].imag };
    }

    H5Sclose(myDataSpace);
    H5Dclose(myDataSet);
#endif
}
#endif

void cReaderH5::readNodeCompoundData(std::vector<PetscInt>& _ids, std::vector<double>& _coord, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef HAVE_HDF5
    std::string datasetname = _nameMainGroup + _nameMatGroup + _nameVecGroup;
    myDataSet = H5Dopen(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

    myDataSpace = H5Dget_space(myDataSet);
    hsize_t dims_out[2];
    int ndims = H5Sget_simple_extent_dims(myDataSpace, dims_out, NULL);

    _ids.resize(dims_out[0]);
    _coord.reserve(dims_out[0] * 3);

    std::vector<elpasoNodal> convertNodalStruct;
    convertNodalStruct.resize(dims_out[0]);

    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(elpasoNodal));
    H5Tinsert(memtype, "Ids", HOFFSET(elpasoNodal, Ids), H5T_NATIVE_INT);
    H5Tinsert(memtype, "xCoords", HOFFSET(elpasoNodal, xCoords), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "yCoords", HOFFSET(elpasoNodal, yCoords), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "zCoords", HOFFSET(elpasoNodal, zCoords), H5T_NATIVE_DOUBLE);

    H5Dread(myDataSet, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &convertNodalStruct[0]);

    for (size_t i = 0; i < convertNodalStruct.size(); i++) {
        _ids[i] = convertNodalStruct[i].Ids;
        _coord.push_back(convertNodalStruct[i].xCoords);
        _coord.push_back(convertNodalStruct[i].yCoords);
        _coord.push_back(convertNodalStruct[i].zCoords);
    }

    H5Sclose(myDataSpace);
    H5Dclose(myDataSet);
#endif
}

int cReaderH5::readIntegerAttributeFromGroup(std::string _name, std::string _path)
{
#ifdef HAVE_HDF5
    try {
        // Save old error handler 
        hid_t old_id;
        H5E_auto2_t old_func;
        //herr_t(*old_func)(void*);
        void* old_client_data;
        H5Eget_auto(old_id, &old_func, &old_client_data);

        // Turn off error handline 
        H5Eset_auto(NULL, NULL, NULL);

        int readBuffer = 0;
        hid_t group_id = H5Gopen(cFileH5::getFileInstance(), _path.c_str(), H5P_DEFAULT);

        hid_t l_attr = H5Aopen(group_id, _name.c_str(), H5P_DEFAULT);
        H5Aread(l_attr, H5T_NATIVE_INT, &readBuffer);

        H5Aclose(l_attr);
        H5Gclose(group_id);

        // Restore previous error handler 
        H5Eset_auto(old_id, old_func, old_client_data);

        return readBuffer;
    }
    catch (...) {
        return 0;
    }
#endif
}

double cReaderH5::readDoubleAttributeFromGroup(std::string _name, std::string _path)
{
#ifdef HAVE_HDF5
    try {
        double readBuffer = 0;
        hid_t group_id = H5Gopen(cFileH5::getFileInstance(), _path.c_str(), H5P_DEFAULT);

        hid_t l_attr = H5Aopen(group_id, _name.c_str(), H5P_DEFAULT);
        H5Aread(l_attr, H5T_NATIVE_DOUBLE, &readBuffer);

        H5Aclose(l_attr);
        H5Gclose(group_id);

        return readBuffer;
    }
    catch (...) {
        std::cout << "Attribute not found" << std::endl;
        return -1000;
    }
#endif
}

double cReaderH5::readDoubleAttributeFromDataset(std::string _name, std::string _path, std::string _dataset)
{
#ifdef HAVE_HDF5
    try {
        // Save old error handler 
        hid_t old_id;
        H5E_auto2_t old_func;
        //herr_t(*old_func)(void*);
        void* old_client_data;
        H5Eget_auto(old_id, &old_func, &old_client_data);

        // Turn off error handline 
        H5Eset_auto(NULL, NULL, NULL);

        double readBuffer = 0;
        std::string datasetname = _path + _dataset;
        myDataSet = H5Dopen(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

        hid_t l_attr = H5Aopen(myDataSet, _name.c_str(), H5P_DEFAULT);
        H5Aread(l_attr, H5T_NATIVE_DOUBLE, &readBuffer);

        H5Aclose(l_attr);
        H5Dclose(myDataSet);

        // Restore previous error handler 
        H5Eset_auto(old_id, old_func, old_client_data);

        return readBuffer;
    }
    catch (...) {
        // ignore
    }
#endif
}

std::string cReaderH5::readStringAttributeFromGroup(std::string _name, std::string _path)
{
#ifdef HAVE_HDF5
    try {
        char* readBuffer;
        hid_t group_id = H5Gopen(cFileH5::getFileInstance(), _path.c_str(), H5P_DEFAULT);

        hid_t l_attr = H5Aopen(group_id, _name.c_str(), H5P_DEFAULT);

        hid_t filetype = H5Aget_type(l_attr);
        size_t sdim = H5Tget_size(filetype);
        sdim++;

        readBuffer = (char*)malloc(sdim * sizeof(char));
        H5Aread(l_attr, filetype, &readBuffer);

        H5Aclose(l_attr);
        H5Gclose(group_id);
        std::string returnstring = readBuffer;
        return returnstring;
    }
    catch (...) {
        std::cout << "Attribute not found" << std::endl;
        return "Not found";
    }
#endif
}

std::string cReaderH5::readStringAttributeFromDataset(std::string _name, std::string _path, std::string _dataset)
{
#ifdef HAVE_HDF5
    // Save old error handler 
    hid_t old_id;
    H5E_auto2_t old_func;
    //herr_t(*old_func)(void*);
    void* old_client_data;
    H5Eget_auto(old_id, &old_func, &old_client_data);

    // Turn off error handline 
    H5Eset_auto(NULL, NULL, NULL);

    char* readBuffer;
    std::string datasetname = _path + _dataset;
    myDataSet = H5Dopen(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

    hid_t l_attr = H5Aopen(myDataSet, _name.c_str(), H5P_DEFAULT);
    hid_t filetype = H5Aget_type(l_attr);
    size_t sdim = H5Tget_size(filetype);
    sdim++;

    readBuffer = (char*)malloc(sdim * sizeof(char));

    H5Aread(l_attr, filetype, &readBuffer);

    H5Aclose(l_attr);
    H5Dclose(myDataSet);

    // Restore previous error handler 
    H5Eset_auto(old_id, old_func, old_client_data);

    std::string strreadbuf = readBuffer;
    return strreadbuf;
#endif
}

int cReaderH5::readIntegerAttributeFromDataset(std::string _name, std::string _path, std::string _dataset)
{
#ifdef HAVE_HDF5
    try {
        // Save old error handler 
        hid_t old_id;
        H5E_auto2_t old_func;
        //herr_t(*old_func)(void*);
        void* old_client_data;
        H5Eget_auto(old_id, &old_func, &old_client_data);

        // Turn off error handline 
        H5Eset_auto(NULL, NULL, NULL);

        int readBuffer = 0;
        std::string datasetname = _path + _dataset;
        myDataSet = H5Dopen(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

        hid_t l_attr = H5Aopen(myDataSet, _name.c_str(), H5P_DEFAULT);
        H5Aread(l_attr, H5T_NATIVE_INT, &readBuffer);

        H5Aclose(l_attr);
        H5Dclose(myDataSet);

        // Restore previous error handler 
        H5Eset_auto(old_id, old_func, old_client_data);

        return readBuffer;
    }
    catch (...) {
        // ignore
    }
#endif
}

int cReaderH5::getNumberOfMembersInGroup(std::string _path)
{
#ifdef HAVE_HDF5
    try {
        hid_t group_id = H5Gopen(cFileH5::getFileInstance(), _path.c_str(), H5P_DEFAULT);
        hsize_t num_members = 0;
        H5Gget_num_objs(group_id, &num_members);

        H5Gclose(group_id);
        return num_members;
    }
    catch (...) {
        std::cout << "Group not found" << std::endl;
        return -1;
    }
#endif
}

PetscInt cReaderH5::getNumberOfLinksInGroup(std::string _path)
{
#ifdef HAVE_HDF5
    hid_t group_id = H5Gopen(cFileH5::getFileInstance(), _path.c_str(), H5P_DEFAULT);
    H5G_info_t ginfo;
    herr_t status = H5Gget_info(group_id, &ginfo);

    H5Gclose(group_id);
    return ginfo.nlinks;
#endif
}

std::string cReaderH5::getMemberName(std::string _path, int _index)
{
#ifdef HAVE_HDF5
    try {
        hid_t group_id = H5Gopen(cFileH5::getFileInstance(), _path.c_str(), H5P_DEFAULT);
        char* name;
        name = new char[600];
        H5Gget_objname_by_idx(group_id, _index, name, 600);

        std::string returnstring(name);

        H5Gclose(group_id);
        return "/" + returnstring;
    }
    catch (...) {
        std::cout << "Group not found" << std::endl;
        return "";
    }
#endif
}

