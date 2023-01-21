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

#include "../../source/shape/shapefunctions.h"
#include "../../source/misc/mytypes.h"

#ifdef HAVE_GTEST
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#endif // HAVE_GTEST

//! @brief Mock class for cShapeFunctions
//! @author Harikrishnan Sreekumar
//! @date 08.01.2023
class cTestShapeFunctions: private cShapeFunctions, public testing::Test
{
public:
    void testShapeFunctionLagrangeQuad4()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        // act
        double retVal4Fun_0 = evaluateShapeFunction(Quadrilateral, 4, 0, N_fun, gpPoint);
        double retVal4Fun_1 = evaluateShapeFunction(Quadrilateral, 4, 1, N_fun, gpPoint);
        double retVal4Fun_2 = evaluateShapeFunction(Quadrilateral, 4, 2, N_fun, gpPoint);
        double retVal4Fun_3 = evaluateShapeFunction(Quadrilateral, 4, 3, N_fun, gpPoint);

        double retVal4Xi_0  = evaluateShapeFunction(Quadrilateral, 4, 0, N_xi, gpPoint);
        double retVal4Xi_1  = evaluateShapeFunction(Quadrilateral, 4, 1, N_xi, gpPoint);
        double retVal4Xi_2  = evaluateShapeFunction(Quadrilateral, 4, 2, N_xi, gpPoint);
        double retVal4Xi_3  = evaluateShapeFunction(Quadrilateral, 4, 3, N_xi, gpPoint);

        double retVal4Eta_0 = evaluateShapeFunction(Quadrilateral, 4, 0, N_eta, gpPoint);
        double retVal4Eta_1 = evaluateShapeFunction(Quadrilateral, 4, 1, N_eta, gpPoint);
        double retVal4Eta_2 = evaluateShapeFunction(Quadrilateral, 4, 2, N_eta, gpPoint);
        double retVal4Eta_3 = evaluateShapeFunction(Quadrilateral, 4, 3, N_eta, gpPoint);

        double retVal4NFun_0 = evaluateShapeFunction(Quadrilateral, 4, 0, N_fun, gpPoint);
        double retError      = evaluateShapeFunction(Quadrilateral, 10, 0, N_fun, gpPoint);
        // assert
        //std::cout<< " "<< retVal4Eta_0 << " " << retVal4Eta_1 << " " << retVal4Eta_2 << " " << retVal4Eta_3 << " " ;
        ASSERT_NEAR(retVal4Fun_0, 0.00482078, 1e-6);
        ASSERT_NEAR(retVal4Fun_1, 0.06461106, 1e-6);
        ASSERT_NEAR(retVal4Fun_2, 0.865957  , 1e-6);
        ASSERT_NEAR(retVal4Fun_3, 0.0646111 , 1e-6);

        ASSERT_NEAR(retVal4Xi_0, -0.03471592, 1e-6);
        ASSERT_NEAR(retVal4Xi_1,  0.03471592, 1e-6);
        ASSERT_NEAR(retVal4Xi_2,  0.465284  , 1e-6);
        ASSERT_NEAR(retVal4Xi_3, -0.465284  , 1e-6);

        ASSERT_NEAR(retVal4Eta_0, -0.03471592, 1e-6);
        ASSERT_NEAR(retVal4Eta_1, -0.465284  , 1e-6);
        ASSERT_NEAR(retVal4Eta_2,  0.465284  , 1e-6);
        ASSERT_NEAR(retVal4Eta_3,  0.0347159 , 1e-6);

        ASSERT_NEAR(retVal4NFun_0, 0.00482078, 1e-6);
        ASSERT_EQ(retError, 0.); 
    }
    
    void testShapeFunctionLagrangeQuad9()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        std::vector<double> retValFun(9), retValXi(9), retValEta(9);
        std::vector<double> refValFun   = { 0.00357488, -0.0479127, 0.642155, -0.0479127, -0.0154525,  0.207103, 0.207103, -0.0154525, 0.0667934};
        std::vector<double> refValXi    = {-0.0215924 , -0.0813827, 1.09074 ,  0.289395 ,  0.102975 ,  0.351778, -1.38014, 0.0933336, -0.445111};
        std::vector<double>  refValEta  = {-0.0215924 ,  0.289395 , 1.09074 , -0.0813827,  0.0933336 , -1.38014, 0.351778, 0.102975, -0.445111};
        // act
        for (size_t i = 0; i < 9; i++)
        {
            retValFun[i]    = evaluateShapeFunction(Quadrilateral, 9, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(Quadrilateral, 9, i, N_xi , gpPoint);
            retValEta[i]    = evaluateShapeFunction(Quadrilateral, 9, i, N_eta, gpPoint);
        }
        // assert
        for (size_t i = 0; i < 9; i++)
        {
            ASSERT_NEAR(retValFun[i], refValFun[i], 1e-4);
            ASSERT_NEAR(retValXi[i] , refValXi[i] , 1e-4);
            ASSERT_NEAR(retValEta[i], refValEta[i], 1e-4);
        }

        /*for (size_t i = 0; i < 9; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < 9; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < 9; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;*/
    }
    
    void testShapeFunctionHermiteQuad4()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        std::vector<double> retValFun(4), retValXi(4), retValEta(4), retValXi2(4), retValEta2(4);
        std::vector<double> refValFun   = { 0.000190244, 0.0136027, 0.972604, 0.0136027};
        std::vector<double> refValXi    = {-0.00267352, 0.0180158, 0.19116, -0.19116};
        std::vector<double> refValEta   = {-0.00267352, -0.19116, 0.19116, 0.00267352};
        std::vector<double> refValXi2   = {0.0178164, 0.0178164, -1.27389, -1.27389};
        std::vector<double> refValEta2  = {0.0178164, 1.27389, -1.27389, -0.0178164};
        // act
        for (size_t i = 0; i < 4; i++)
        {
            retValFun[i]    = HermiteQuadN(4, i, N_fun, gpPoint);
            retValXi[i]     = HermiteQuadN(4, i, N_xi , gpPoint);
            retValEta[i]    = HermiteQuadN(4, i, N_eta, gpPoint);
            retValXi2[i]    = HermiteQuadN(4, i, N_xi2, gpPoint);
            retValEta2[i]   = HermiteQuadN(4, i, N_eta2, gpPoint);
        }
        // assert
        for (size_t i = 0; i < 4; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
            ASSERT_NEAR(retValEta[i]    , refValEta[i]  , 1e-4);
            ASSERT_NEAR(retValXi2[i]    , refValXi2[i]  , 1e-4);
            ASSERT_NEAR(retValEta2[i]   , refValEta2[i] , 1e-4);
        }

        /*for (size_t i = 0; i < 4; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < 4; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < 4; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;
        for (size_t i = 0; i < 4; i++)
            std::cout << ", "<< retValXi2[i];
        std::cout << std::endl;
        for (size_t i = 0; i < 4; i++)
            std::cout << ", "<< retValEta2[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }
    
    void testShapeFunctionSerendipityQuad8()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        std::vector<double> retValFun(8), retValXi(8), retValEta(8);
        std::vector<double> refValFun   = { -0.0131235, -0.0646111, 0.625457, -0.0646111, 0.0179443, 0.2405, 0.2405, 0.0179443};
        std::vector<double> refValXi    = {0.0896854, 0.0298951, 1.20202, 0.400673, -0.119581, 0.129222, -1.60269, -0.129222};
        std::vector<double> refValEta   = {0.0896854, 0.400673, 1.20202, 0.0298951, -0.129222, -1.60269, 0.129222, -0.119581};
        // act
        for (size_t i = 0; i < 8; i++)
        {
            retValFun[i]    = evaluateShapeFunction(QuadSerendipity, 8, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(QuadSerendipity, 8, i, N_xi , gpPoint);
            retValEta[i]    = evaluateShapeFunction(QuadSerendipity, 8, i, N_eta, gpPoint);
        }
        // assert
        for (size_t i = 0; i < 8; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
            ASSERT_NEAR(retValEta[i]    , refValEta[i]  , 1e-4);
        }

        /*for (size_t i = 0; i < 8; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < 8; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < 8; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }
    
    void testShapeFunctionLagrangeHex20()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        std::vector<double> retValFun(20), retValXi(20), retValEta(20), retValZeta(20);
        std::vector<double> refValFun   = { -0.00897213, -0.0646111, -0.12025, -0.0646111, -0.00897213, -0.0646111, -0.12025, -0.0646111, 0.00897213, 0.12025, 0.12025, 0.00897213, 0.00482078, 0.0646111, 0.865957, 0.0646111, 0.00897213, 0.12025, 0.12025, 0.00897213};
        std::vector<double> refValXi    = {0.0622007, 0.432979, 0.368367, -0.00241039, 0.0622007, 0.432979, 0.368367, -0.00241039, -0.0646111, -0.801346, 0.0646111, -0.0597903, -0.0347159, -0.465284, 0.465284, 0.0347159, -0.0646111, -0.801346, 0.0646111, -0.0597903};
        std::vector<double> refValEta   = { 0.0622007, -0.00241039, 0.368367, 0.432979, 0.0622007, -0.00241039, 0.368367, 0.432979, -0.0597903, 0.0646111, -0.801346, -0.0646111, -0.0347159, 0.0347159, 0.465284, -0.465284, -0.0597903, 0.0646111, -0.801346, -0.0646111};
        std::vector<double> refValZeta  = { -0.00656174, -0.0323055, 0.312729, -0.0323055, 0.00656174, 0.0323055, -0.312729, 0.0323055, 0.00897213, 0.12025, 0.12025, 0.00897213, -0, -0, -0, -0, -0.00897213, -0.12025, -0.12025, -0.00897213};
        // act
        for (size_t i = 0; i < 20; i++)
        {
            retValFun[i]    = evaluateShapeFunction(Hexahedron, 20, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(Hexahedron, 20, i, N_xi , gpPoint);
            retValEta[i]    = evaluateShapeFunction(Hexahedron, 20, i, N_eta, gpPoint);
            retValZeta[i]   = evaluateShapeFunction(Hexahedron, 20, i, N_zeta, gpPoint);
        }
        // assert
        for (size_t i = 0; i < 20; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
            ASSERT_NEAR(retValEta[i]    , refValEta[i]  , 1e-4);
            ASSERT_NEAR(retValZeta[i]   , refValZeta[i]  , 1e-4);
        }

        /*for (size_t i = 0; i < 20; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < 20; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < 20; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;
        for (size_t i = 0; i < 20; i++)
            std::cout << ", "<< retValZeta[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }
    
    void testShapeFunctionLagrangeHex8()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        int nNode = 8;
        std::vector<double> retValFun(nNode), retValXi(nNode), retValEta(nNode), retValZeta(nNode);
        std::vector<double> refValFun   = { 0.00241039, 0.0323055, 0.432979, 0.0323055, 0.00241039, 0.0323055, 0.432979, 0.0323055};
        std::vector<double> refValXi    = {-0.017358, 0.017358, 0.232642, -0.232642, -0.017358, 0.017358, 0.232642, -0.232642};
        std::vector<double> refValEta   = { -0.017358, -0.232642, 0.232642, 0.017358, -0.017358, -0.232642, 0.232642, 0.017358};
        std::vector<double> refValZeta  = { -0.00241039, -0.0323055, -0.432979, -0.0323055, 0.00241039, 0.0323055, 0.432979, 0.0323055};
        // act
        for (size_t i = 0; i < nNode; i++)
        {
            retValFun[i]    = evaluateShapeFunction(Hexahedron, nNode, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(Hexahedron, nNode, i, N_xi , gpPoint);
            retValEta[i]    = evaluateShapeFunction(Hexahedron, nNode, i, N_eta, gpPoint);
            retValZeta[i]   = evaluateShapeFunction(Hexahedron, nNode, i, N_zeta, gpPoint);
        }
        // assert
        for (size_t i = 0; i < nNode; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
            ASSERT_NEAR(retValEta[i]    , refValEta[i]  , 1e-4);
            ASSERT_NEAR(retValZeta[i]   , refValZeta[i]  , 1e-4);
        }

        /*for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValZeta[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }
    
    void testShapeFunctionLagrangeHex27()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        int nNode = 27;
        std::vector<double> retValFun(nNode), retValXi(nNode), retValEta(nNode), retValZeta(nNode);
        std::vector<double> refValFun   = { -0, 0, -0, 0, 0, -0, 0, -0, 0, -0, -0, 0, 0.00357488, -0.0479127, 0.642155, -0.0479127, -0, 0, 0, -0, 0.0667934, -0, 0, -0.0154525, 0.207103, -0.0154525, 0.207103};
        std::vector<double> refValXi    = { 0, 0, -0, -0, -0, -0, 0, 0, -0, -0, 0, -0, -0.0215924, -0.0813827, 1.09074, 0.289395, 0, 0, -0, 0, -0.445111, 0, -0, 0.0933336, 0.351778, 0.102975, -1.38014};
        std::vector<double> refValEta   = { 0, -0, -0, 0, -0, 0, 0, -0, -0, 0, -0, -0, -0.0215924, 0.289395, 1.09074, -0.0813827, 0, -0, 0, 0, -0.445111, 0, -0, 0.102975, -1.38014, 0.0933336, 0.351778};
        std::vector<double> refValZeta  = { -0.00178744, 0.0239564, -0.321078, 0.0239564, 0.00178744, -0.0239564, 0.321078, -0.0239564, 0.00772623, -0.103552, -0.103552, 0.00772623, -0, 0, -0, 0, -0.00772623, 0.103552, 0.103552, -0.00772623, -0, -0.0333967, 0.0333967, 0, -0, 0, -0};
        // act
        for (size_t i = 0; i < nNode; i++)
        {
            retValFun[i]    = evaluateShapeFunction(Hexahedron, nNode, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(Hexahedron, nNode, i, N_xi , gpPoint);
            retValEta[i]    = evaluateShapeFunction(Hexahedron, nNode, i, N_eta, gpPoint);
            retValZeta[i]   = evaluateShapeFunction(Hexahedron, nNode, i, N_zeta, gpPoint);
        }
        // assert
        for (size_t i = 0; i < nNode; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
            ASSERT_NEAR(retValEta[i]    , refValEta[i]  , 1e-4);
            ASSERT_NEAR(retValZeta[i]   , refValZeta[i]  , 1e-4);
        }

        /*for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValZeta[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }

    void testShapeFunctionTria3()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        int nNode = 3;
        std::vector<double> retValFun(nNode), retValXi(nNode), retValEta(nNode);
        std::vector<double> refValFun   = { -0.722273, 0.861136, 0.861136};
        std::vector<double> refValXi    = { -1, 1, 0};
        std::vector<double> refValEta   = { -1, 0, 1};
        // act
        for (size_t i = 0; i < nNode; i++)
        {
            retValFun[i]    = evaluateShapeFunction(Tria, nNode, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(Tria, nNode, i, N_xi , gpPoint);
            retValEta[i]    = evaluateShapeFunction(Tria, nNode, i, N_eta, gpPoint);
        }
        // assert
        for (size_t i = 0; i < nNode; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
            ASSERT_NEAR(retValEta[i]    , refValEta[i]  , 1e-4);
        }

        /*for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }
    
    void testShapeFunctionTria6()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        int nNode = 6;
        std::vector<double> retValFun(nNode), retValXi(nNode), retValEta(nNode);
        std::vector<double> refValFun   = { 0.621975, 0.621975, -0, 2.96622, 0, 0};
        std::vector<double> refValXi    = { -2.44455, 2.44455, 0, 0, 0, -0};
        std::vector<double> refValEta   = { -2.44455, 0, -1, -3.44455, 3.44455, 3.44455};
        // act
        for (size_t i = 0; i < nNode; i++)
        {
            retValFun[i]    = evaluateShapeFunction(Tria, nNode, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(Tria, nNode, i, N_xi , gpPoint);
            retValEta[i]    = evaluateShapeFunction(Tria, nNode, i, N_eta, gpPoint);
        }
        // assert
        for (size_t i = 0; i < nNode; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
            ASSERT_NEAR(retValEta[i]    , refValEta[i]  , 1e-4);
        }

        /*for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }
    
    void testShapeFunctionTria9()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        int nNode = 9;
        std::vector<double> retValFun(nNode), retValXi(nNode), retValEta(nNode);
        std::vector<double> refValFun   = { 0.397748, 0.397748, 0, 5.28384, 5.28384, 0, -0, -0, 0};
        std::vector<double> refValXi    = { -3.26078, 3.26078, 0, -10.011, 10.011, 0, 0, -0, -0};
        std::vector<double> refValEta   = { -3.26078, 0, 1, -16.1469, -6.13589, 6.13589, -3.87511, -3.87511, 6.13589};
        // act
        for (size_t i = 0; i < nNode; i++)
        {
            retValFun[i]    = evaluateShapeFunction(Tria, nNode, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(Tria, nNode, i, N_xi , gpPoint);
            retValEta[i]    = evaluateShapeFunction(Tria, nNode, i, N_eta, gpPoint);
        }
        // assert
        for (size_t i = 0; i < nNode; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
            ASSERT_NEAR(retValEta[i]    , refValEta[i]  , 1e-4);
        }

        /*for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }
    
    void testShapeFunctionTet4()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        int nNode = 4;
        std::vector<double> retValFun(nNode), retValXi(nNode), retValEta(nNode), retValZeta(nNode);
        std::vector<double> refValFun   = { -0.722273, 0.861136, 0.861136, 0};
        std::vector<double> refValXi    = { -1, 1, 0, 0};
        std::vector<double> refValEta   = { -1, 0, 1, 0};
        std::vector<double> refValZeta   = { -1, 0, 0, 1};
        // act
        for (size_t i = 0; i < nNode; i++)
        {
            retValFun[i]    = evaluateShapeFunction(Tetrahedron,  nNode, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(Tetrahedron, nNode, i, N_xi , gpPoint);
            retValEta[i]    = evaluateShapeFunction(Tetrahedron, nNode, i, N_eta, gpPoint);
            retValZeta[i]    = evaluateShapeFunction(Tetrahedron, nNode, i, N_zeta, gpPoint);
        }
        // assert
        for (size_t i = 0; i < nNode; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
            ASSERT_NEAR(retValEta[i]    , refValEta[i]  , 1e-4);
            ASSERT_NEAR(retValZeta[i]    , refValZeta[i]  , 1e-4);
        }

        /*for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValZeta[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }
    
    void testShapeFunctionTet10()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        int nNode = 10;
        std::vector<double> retValFun(nNode), retValXi(nNode), retValEta(nNode), retValZeta(nNode);
        std::vector<double> refValFun   = { 1.76563, 0.621975, 0.621975, -0, -2.4879, 2.96622, -2.4879, -0, 0, 0, 1.76563, 0.621975, 0.621975, -0, -2.4879, 2.96622, -2.4879, -0, 0, 0};
        std::vector<double> refValXi    = { 3.88909, 2.44455, 0, 0, -6.33364, 3.44455, -3.44455, -0, 0, 0};
        std::vector<double> refValEta   = { 3.88909, 0, 2.44455, 0, -3.44455, 3.44455, -6.33364, -0, 0, 0};
        std::vector<double> refValZeta   = { 3.88909, 0, 0, -1, -3.44455, 0, -3.44455, -2.88909, 3.44455, 3.44455};
        // act
        for (size_t i = 0; i < nNode; i++)
        {
            retValFun[i]    = evaluateShapeFunction(Tetrahedron, nNode, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(Tetrahedron, nNode, i, N_xi , gpPoint);
            retValEta[i]    = evaluateShapeFunction(Tetrahedron, nNode, i, N_eta, gpPoint);
            retValZeta[i]    = evaluateShapeFunction(Tetrahedron, nNode, i, N_zeta, gpPoint);
        }
        // assert
        for (size_t i = 0; i < nNode; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
            ASSERT_NEAR(retValEta[i]    , refValEta[i]  , 1e-4);
            ASSERT_NEAR(retValZeta[i]    , refValZeta[i]  , 1e-4);
        }

        /*for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValZeta[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }
    
    void testShapeFunctionTet16()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        int nNode = 16;
        std::vector<double> retValFun(nNode), retValXi(nNode), retValEta(nNode), retValZeta(nNode);
        std::vector<double> refValFun   = { -4.76539, 0.397748, 0.397748, 0, 8.86357, -4.43178, 5.28384, 5.28384, -4.43178, 8.86357, 0, 0, 0, 0, -0, -0, -4.76539, 0.397748, 0.397748, 0, 8.86357, -4.43178, 5.28384, 5.28384, -4.43178, 8.86357, 0, 0, 0, 0, -0, -0};
        std::vector<double> refValXi    = { -14.5431, 3.26078, 0, 0, 30.9613, -19.679, 16.1469, 6.13589, -6.13589, 20.6684, 0, 0, 0, -0, 0, 0};
        std::vector<double> refValEta   = { -14.5431, 0, 3.26078, 0, 20.6684, -6.13589, 6.13589, 16.1469, -19.679, 30.9613, 0, 0, 0, -0, 0, 0};
        std::vector<double> refValZeta   = { -14.5431, 0, 0, 1, 20.6684, -6.13589, 0, 0, -6.13589, 20.6684, 10.2929, 6.13589, 6.13589, 3.25023, -3.87511, -3.87511};
        // act
        for (size_t i = 0; i < nNode; i++)
        {
            retValFun[i]    = evaluateShapeFunction(Tetrahedron, nNode, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(Tetrahedron, nNode, i, N_xi , gpPoint);
            retValEta[i]    = evaluateShapeFunction(Tetrahedron, nNode, i, N_eta, gpPoint);
            retValZeta[i]    = evaluateShapeFunction(Tetrahedron, nNode, i, N_zeta, gpPoint);
        }
        // assert
        for (size_t i = 0; i < nNode; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
            ASSERT_NEAR(retValEta[i]    , refValEta[i]  , 1e-4);
            ASSERT_NEAR(retValZeta[i]    , refValZeta[i]  , 1e-4);
        }

        /*for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValZeta[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }
    
    void testShapeFunctionTet4L()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        int nNode = 4;
        std::vector<double> retValFun(nNode), retValXi(nNode), retValEta(nNode), retValZeta(nNode);
        std::vector<double> refValFun   = { 0.861136, 0.861136, 0, 3.79474e-316};
        std::vector<double> refValXi    = { -1, 1, 0, 0};
        std::vector<double> refValEta   = { -1, 0, 1, 0};
        std::vector<double> refValZeta   = { -1, 0, 0, 1};
        // act
        for (size_t i = 0; i < nNode; i++)
        {
            retValFun[i]    = evaluateShapeFunction(Tetrahedron, 41, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(Tetrahedron, 41, i, N_xi , gpPoint);
            retValEta[i]    = evaluateShapeFunction(Tetrahedron, 41, i, N_eta, gpPoint);
            retValZeta[i]    = evaluateShapeFunction(Tetrahedron, 41, i, N_zeta, gpPoint);
        }
        // assert
        for (size_t i = 0; i < nNode; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
            ASSERT_NEAR(retValEta[i]    , refValEta[i]  , 1e-4);
            ASSERT_NEAR(retValZeta[i]    , refValZeta[i]  , 1e-4);
        }

        /*for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValZeta[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }
    
    void testShapeFunctionTet10L()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        int nNode = 10;
        std::vector<double> retValFun(nNode), retValXi(nNode), retValEta(nNode), retValZeta(nNode);
        std::vector<double> refValFun   = { 0.621975, 0.621975, -0, -0, 2.96622, 0, 0, 0, 0, 0};
        std::vector<double> refValXi    = { -2.44455, 2.44455, 0, 0, 0, 0, -0, -0, 0, 0};
        std::vector<double> refValEta   = { -2.44455, 0, -1, 0, -3.44455, 3.44455, 3.44455, -0, 0, 0};
        std::vector<double> refValZeta   = { -2.44455, 0, 0, -1, -3.44455, 0, -0, 3.44455, 3.44455, 0};
        // act
        for (size_t i = 0; i < nNode; i++)
        {
            retValFun[i]    = evaluateShapeFunction(Tetrahedron, 101, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(Tetrahedron, 101, i, N_xi , gpPoint);
            retValEta[i]    = evaluateShapeFunction(Tetrahedron, 101, i, N_eta, gpPoint);
            retValZeta[i]    = evaluateShapeFunction(Tetrahedron, 101, i, N_zeta, gpPoint);
        }
        // assert
        for (size_t i = 0; i < nNode; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
            ASSERT_NEAR(retValEta[i]    , refValEta[i]  , 1e-4);
            ASSERT_NEAR(retValZeta[i]    , refValZeta[i]  , 1e-4);
        }

        /*for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValZeta[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }
    
    void testShapeFunctionTet16L()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        int nNode = 16;
        std::vector<double> retValFun(nNode), retValXi(nNode), retValEta(nNode), retValZeta(nNode);
        std::vector<double> refValFun   = { 0.397748, 0.397748, 0, 0, 5.28384, 5.28384, 0, -0, -0, 0, 0, 0, -0, -0, -0, -0};
        std::vector<double> refValXi    = { -3.26078, 3.26078, 0, 0, -10.011, 10.011, 0, 0, -0, -0, -0, 0, 0, -0, 0, 0};
        std::vector<double> refValEta   = { -3.26078, 0, 1, 0, -16.1469, -6.13589, 6.13589, -3.87511, -3.87511, 6.13589, -0, 0, 0, -0, 0, 0};
        std::vector<double> refValZeta   = { -3.26078, 0, 0, 1, -16.1469, -6.13589, 0, 0, -0, -0, 6.13589, 6.13589, 0, -3.87511, -3.87511, 0};
        // act
        for (size_t i = 0; i < nNode; i++)
        {
            retValFun[i]    = evaluateShapeFunction(Tetrahedron, 161, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(Tetrahedron, 161, i, N_xi , gpPoint);
            retValEta[i]    = evaluateShapeFunction(Tetrahedron, 161, i, N_eta, gpPoint);
            retValZeta[i]    = evaluateShapeFunction(Tetrahedron, 161, i, N_zeta, gpPoint);
        }
        // assert
        for (size_t i = 0; i < nNode; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
            ASSERT_NEAR(retValEta[i]    , refValEta[i]  , 1e-4);
            ASSERT_NEAR(retValZeta[i]    , refValZeta[i]  , 1e-4);
        }

        /*for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValEta[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValZeta[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }
    
    void testShapeFunctionBeam2()
    {
        // arrange
        cGaussPoints gp;
        cPoint gpPoint( gp.getGaussPoint2D(4, 0) );
        int nNode = 2;
        std::vector<double> retValFun(nNode), retValXi(nNode);
        std::vector<double> refValFun   = { 0.0694318, 0.930568};
        std::vector<double> refValXi    = { -0.5, 0.5};
        // act
        for (size_t i = 0; i < nNode; i++)
        {
            retValFun[i]    = evaluateShapeFunction(Beam, nNode, i, N_fun, gpPoint);
            retValXi[i]     = evaluateShapeFunction(Beam, nNode, i, N_xi , gpPoint);
        }
        // assert
        for (size_t i = 0; i < nNode; i++)
        {
            ASSERT_NEAR(retValFun[i]    , refValFun[i]  , 1e-4);
            ASSERT_NEAR(retValXi[i]     , refValXi[i]   , 1e-4);
        }

        /*for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValFun[i];
        std::cout << std::endl;
        for (size_t i = 0; i < nNode; i++)
            std::cout << ", "<< retValXi[i];
        std::cout << std::endl;
        ASSERT_EQ(1,2);*/
    }
};

TEST_F(cTestShapeFunctions, correctShapeFunctionLagrangeQuad4)
{
    testShapeFunctionLagrangeQuad4();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionLagrangeQuad9)
{
    testShapeFunctionLagrangeQuad9();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionHermiteQuad4)
{
    testShapeFunctionHermiteQuad4();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionSerendipityQuad8)
{
    testShapeFunctionSerendipityQuad8();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionLagrangeHex8)
{
    testShapeFunctionLagrangeHex8();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionLagrangeHex27)
{
    testShapeFunctionLagrangeHex27();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionLagrangeHex20)
{
    testShapeFunctionLagrangeHex20();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionTria3)
{
    testShapeFunctionTria3();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionTria6)
{
    testShapeFunctionTria6();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionTria9)
{
    testShapeFunctionTria9();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionTet4)
{
    testShapeFunctionTet4();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionTet10)
{
    testShapeFunctionTet10();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionTet16)
{
    testShapeFunctionTet16();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionTet4L)
{
    testShapeFunctionTet4L();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionTet10L)
{
    testShapeFunctionTet10L();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionTet16L)
{
    testShapeFunctionTet16L();
}

TEST_F(cTestShapeFunctions, correctShapeFunctionBeam2)
{
    testShapeFunctionBeam2();
}