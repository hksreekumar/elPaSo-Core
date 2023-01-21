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

#include "outputh5.h"
#include "../../basics/mpi/mpitools.h"
#include<iostream>

#define RANK_VEC   1
#define RANK_MAT   2

cOutputH5::cOutputH5()
{
#ifdef HAVE_HDF5
    H5Eset_auto(NULL, NULL, NULL);
#endif // HAVE_HDF5

}

cOutputH5::~cOutputH5()
{
#ifdef HAVE_HDF5
    //H5Gclose(myGroup);
    //H5Dclose(myDataSet);
    //H5Sclose(myDataSpace);
#endif
}

void cOutputH5::createGroup(std::string _nameMainGroup)
{
#ifdef HAVE_HDF5
    herr_t status = H5Gget_objinfo(cFileH5::getFileInstance(), _nameMainGroup.c_str(), 0, NULL);
        
    if (status == 0)
        myGroup = H5Gopen(cFileH5::getFileInstance(), _nameMainGroup.c_str(), H5P_DEFAULT);
    else 
        myGroup = H5Gcreate(cFileH5::getFileInstance(), _nameMainGroup.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(myGroup);
#endif
}

void cOutputH5::createSecoGroup(std::string _nameMainGroup, std::string _nameSecoGroup)
{
#ifdef HAVE_HDF5
    std::string groupname = _nameMainGroup + _nameSecoGroup;
    herr_t status = H5Gget_objinfo(cFileH5::getFileInstance(), groupname.c_str(), 0, NULL);

    if (status == 0) 
        myGroup = H5Gopen(cFileH5::getFileInstance(), groupname.c_str(), H5P_DEFAULT);
    else 
        myGroup = H5Gcreate(cFileH5::getFileInstance(), groupname.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Gclose(myGroup);
#endif
}

void cOutputH5::appendIntegerVector(std::vector<int>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef HAVE_HDF5
        hsize_t dimsf[2];
        dimsf[0] = _vector.size();
        myDataSpace = H5Screate_simple(RANK_VEC, dimsf, NULL);
        std::string datasetname = _nameMainGroup + _nameMatGroup + _nameVecGroup;
                
        myDataSet = H5Dcreate2(cFileH5::getFileInstance(), datasetname.c_str(), H5T_NATIVE_INT, myDataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        
        H5Dwrite(myDataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _vector.data());

        H5Dclose(myDataSet);
        H5Sclose(myDataSpace);
#endif
}

void cOutputH5::writeIntegerVectorAttributeToGroup(std::vector<int>& _vector, std::string _pathGroup, std::string _attribute)
{
#ifdef HAVE_HDF5
    myGroup = H5Gopen(cFileH5::getFileInstance(), _pathGroup.c_str(), H5P_DEFAULT);

    // delete existing attribute
    htri_t status = H5Aexists(myGroup, _attribute.c_str());
    if (status != 0)
        H5Adelete(myGroup, _attribute.c_str());

    hsize_t dimsf[2];
    dimsf[0] = _vector.size();
    myDataSpace = H5Screate_simple(RANK_VEC, dimsf, NULL);

    hid_t l_attr = H5Acreate2(myGroup, _attribute.c_str(), H5T_NATIVE_INT, myDataSpace, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(l_attr, H5T_NATIVE_INT, _vector.data());

    H5Aclose(l_attr);
    H5Sclose(myDataSpace);
    H5Gclose(myGroup);
#endif
}

#ifdef PETSC_USE_COMPLEX
void cOutputH5::appendComplexVector(std::vector<PetscComplex>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef HAVE_HDF5
    std::string datasetname = _nameMainGroup + _nameMatGroup + _nameVecGroup;

    hsize_t dimsf[2];
    dimsf[0] = _vector.size();
    myDataSpace = H5Screate_simple(RANK_VEC, dimsf, NULL);

    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(elpasoComplexDouble));

    std::vector<elpasoComplexDouble> convertPetscComplex;
    convertPetscComplex.resize(_vector.size());
    for (size_t i = 0; i < _vector.size(); i++)
        convertPetscComplex[i] = {_vector[i].real(), _vector[i].imag()};

    H5Tinsert(memtype, "real", HOFFSET(elpasoComplexDouble, real), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "imag", HOFFSET(elpasoComplexDouble, imag), H5T_NATIVE_DOUBLE);

    htri_t status = H5Lexists(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

    if (status != 0) {
        // delete dataset
        H5Ldelete(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);
    }
    //if (status == 0)
        myDataSet = H5Dcreate2(cFileH5::getFileInstance(), datasetname.c_str(), memtype, myDataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //else 
    //    myDataSet = H5Dopen(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);    

    H5Dwrite(myDataSet, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, convertPetscComplex.data());

    H5Sclose(myDataSpace);
    H5Dclose(myDataSet);
#endif
}
#endif


#ifdef PETSC_USE_COMPLEX
void cOutputH5::appendComplexVector(std::vector<PetscComplex>& _vector, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup, PetscInt _start, PetscInt _end)
{
#ifdef HAVE_HDF5
    std::string datasetname = _nameMainGroup + _nameMatGroup + _nameVecGroup;

    hid_t       filespace, memspace;      /* file and memory dataspace identifiers */
    hsize_t	count[2];	                  /* hyperslab selection parameters */
    hsize_t	offset[2];
    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    count[0] = _end - _start;
    offset[0] = _start;

    hsize_t dimsf[2];
    dimsf[0] = _vector.size();
    myDataSpace = H5Screate_simple(RANK_VEC, dimsf, NULL);

    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(elpasoComplexDouble));

    std::vector<elpasoComplexDouble> convertPetscComplex;
    convertPetscComplex.resize(count[0]);
    for (size_t i = 0; i < count[0]; i++)
        convertPetscComplex[i] = { _vector[i].real(), _vector[i].imag() };

    H5Tinsert(memtype, "real", HOFFSET(elpasoComplexDouble, real), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "imag", HOFFSET(elpasoComplexDouble, imag), H5T_NATIVE_DOUBLE);

    htri_t status = H5Lexists(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

    if (status != 0) {
        // delete dataset
        H5Ldelete(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);
    }
    myDataSet = H5Dcreate2(cFileH5::getFileInstance(), datasetname.c_str(), memtype, myDataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /*
     * Each process defines dataset in memory and writes it to the hyperslab
     * in the file.
     */
    memspace = H5Screate_simple(RANK_VEC, count, NULL);

    /*
     * Select hyperslab in the file.
     */
    filespace = H5Dget_space(myDataSet);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);

    hid_t plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    H5Dwrite(myDataSet, memtype, memspace, filespace, plist_id, convertPetscComplex.data());

    H5Pclose(plist_id);
    H5Sclose(filespace);
    H5Sclose(memspace);
    H5Sclose(myDataSpace);
    H5Dclose(myDataSet);
#endif
}
#endif

void cOutputH5::appendComplexMatrix(std::vector<PetscScalar>& _vectorMatrix, size_t _m, size_t _n, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup, PetscBool _isRowMajor)
{
#ifdef PETSC_USE_COMPLEX
#ifdef HAVE_HDF5
    std::string datasetname = _nameMainGroup + _nameMatGroup + _nameVecGroup;

    hsize_t dimsf[2];
    dimsf[0] = _m;
    dimsf[1] = _n;

    myDataSpace = H5Screate_simple(RANK_MAT, dimsf, NULL);

    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(elpasoComplexDouble));

    std::vector<elpasoComplexDouble> convertPetscComplex;
    convertPetscComplex.reserve(_vectorMatrix.size());

    if (_isRowMajor)
    {
        for (size_t i = 0; i < _vectorMatrix.size(); i++)
            convertPetscComplex.push_back({ _vectorMatrix[i].real(), _vectorMatrix[i].imag() });
    }
    else {
        // Colum-major storing
        for (size_t iRow = 0; iRow < _m; iRow++)
            for (size_t iCol = 0; iCol < _n; iCol++)
            {
                PetscInt index = iRow + iCol*_m;
                convertPetscComplex.push_back({ _vectorMatrix[index].real(), _vectorMatrix[index].imag() });
            }
    }
    
    H5Tinsert(memtype, "real", HOFFSET(elpasoComplexDouble, real), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "imag", HOFFSET(elpasoComplexDouble, imag), H5T_NATIVE_DOUBLE);

    htri_t status = H5Lexists(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

    if (status != 0) {
        herr_t delstatus = H5Ldelete(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);
    }
    myDataSet = H5Dcreate2(cFileH5::getFileInstance(), datasetname.c_str(), memtype, myDataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    H5Dwrite(myDataSet, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, convertPetscComplex.data());

    H5Tclose(memtype);
    H5Sclose(myDataSpace);
    H5Dclose(myDataSet);

    // write dimensions
    writeIntegerAttributeToDataset("numrows", _nameMainGroup + _nameMatGroup, _m, _nameVecGroup);
    writeIntegerAttributeToDataset("numcols", _nameMainGroup + _nameMatGroup, _n, _nameVecGroup);
#endif
#endif
}

void cOutputH5::writeStringAttributeToGroup(std::string _name, std::string _path, std::string _value)
{
#ifdef HAVE_HDF5
    myGroup = H5Gopen(cFileH5::getFileInstance(), _path.c_str(), H5P_DEFAULT);

    // open existing attribute in case if exists
    htri_t status = H5Aexists(myGroup, _name.c_str());
    hid_t l_attr;
    hid_t attr_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(attr_type, H5T_VARIABLE);
    H5Tset_cset(attr_type, H5T_CSET_UTF8);

    hid_t attr_space = H5Screate(H5S_SCALAR);

    if (status != 0) {
        l_attr = H5Aopen(myGroup, _name.c_str(), H5P_DEFAULT);
    }
    else {
        l_attr = H5Acreate2(myGroup, _name.c_str(), attr_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);
    }

    H5Awrite(l_attr, attr_type, &_value);

    H5Aclose(l_attr);
    H5Sclose(attr_space);
    H5Gclose(myGroup);
#endif
}

void cOutputH5::writeIntegerAttributeToGroup(std::string _name, std::string _path, int _value)
{
#ifdef HAVE_HDF5
    myGroup = H5Gopen(cFileH5::getFileInstance(), _path.c_str(), H5P_DEFAULT);

    // delete existing attribute
    htri_t status = H5Aexists(myGroup, _name.c_str());
    if (status != 0)
        H5Adelete(myGroup, _name.c_str());

    hid_t attr_type = H5Tcopy(H5T_NATIVE_INT);
    H5Tset_cset(attr_type, H5T_CSET_UTF8);

    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t l_attr = H5Acreate2(myGroup, _name.c_str(), attr_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(l_attr, attr_type, &_value);

    H5Tclose(attr_type);
    H5Sclose(attr_space);
    H5Aclose(l_attr);
    H5Gclose(myGroup);
#endif
}

void cOutputH5::writeIntegerAttributeToDataset(std::string _name, std::string _path, int _value, std::string _dataset)
{
#ifdef HAVE_HDF5
    std::string datasetname = _path + _dataset;
    myDataSet = H5Dopen(cFileH5::getFileInstance(), datasetname.c_str(), H5P_DEFAULT);

    // delete existing attribute
    htri_t status = H5Aexists(myDataSet, _name.c_str());
    if (status != 0)
        H5Adelete(myDataSet, _name.c_str());

    hid_t attr_type = H5Tcopy(H5T_NATIVE_INT);
    H5Tset_cset(attr_type, H5T_CSET_UTF8);

    hid_t attr_space = H5Screate(H5S_SCALAR);
    hid_t l_attr = H5Acreate2(myDataSet, _name.c_str(), attr_type, attr_space, H5P_DEFAULT, H5P_DEFAULT);

    H5Awrite(l_attr, attr_type, &_value);

    H5Tclose(attr_type);
    H5Sclose(attr_space);
    H5Aclose(l_attr);
    H5Dclose(myDataSet);
#endif
}

void cOutputH5::appendRealDenseMatrix(std::vector<PetscReal>& _matrix, int _m, int _n, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef HAVE_HDF5
    std::string datasetname = _nameMainGroup + _nameMatGroup + _nameVecGroup;
    hsize_t dimsf[2];
    dimsf[0] = _m;
    dimsf[1] = _n;
    myDataSpace = H5Screate_simple(RANK_MAT, dimsf, NULL);

    myDataSet = H5Dcreate2(cFileH5::getFileInstance(), datasetname.c_str(), H5T_NATIVE_DOUBLE, myDataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(myDataSet, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, _matrix.data());

    H5Sclose(myDataSpace);
    H5Dclose(myDataSet);
#endif
}

void cOutputH5::appendIntegerDenseMatrix(std::vector<PetscInt>& _matrix, int _m, int _n, std::string _nameMainGroup, std::string _nameMatGroup, std::string _nameVecGroup)
{
#ifdef HAVE_HDF5
    std::string datasetname = _nameMainGroup + _nameMatGroup + _nameVecGroup;
    hsize_t dimsf[2];
    dimsf[0] = _m;
    dimsf[1] = _n;
    myDataSpace = H5Screate_simple(RANK_MAT, dimsf, NULL);

    myDataSet = H5Dcreate2(cFileH5::getFileInstance(), datasetname.c_str(), H5T_NATIVE_INT, myDataSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    H5Dwrite(myDataSet, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, _matrix.data());

    H5Sclose(myDataSpace);
    H5Dclose(myDataSet);
#endif
}