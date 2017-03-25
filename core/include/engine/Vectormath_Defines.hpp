#pragma once
#ifndef VECTORMATH_DEFINES_H
#define VECTORMATH_DEFINES_H

#include <Eigen/Core>

#include <vector>
#include <array>

#include "Spirit_Defines.h"

// Dynamic Eigen typedefs
typedef Eigen::Matrix<scalar, -1,  1> VectorX;
typedef Eigen::Matrix<scalar,  1, -1> RowVectorX;
typedef Eigen::Matrix<scalar, -1, -1> MatrixX;

// 3D Eigen typedefs
typedef Eigen::Matrix<scalar, 3, 1> Vector3;
typedef Eigen::Matrix<scalar, 1, 3> RowVector3;
typedef Eigen::Matrix<scalar, 3, 3> Matrix3;


// Vectorfield and Scalarfield typedefs
#ifdef USE_CUDA
    #include "Managed_Allocator.hpp"
    typedef std::vector<bool,       managed_allocator<bool>>       boolfield;
    typedef std::vector<int,        managed_allocator<int>>        intfield;
    typedef std::vector<scalar,     managed_allocator<scalar>>     scalarfield;
    typedef std::vector<Vector3,    managed_allocator<Vector3>>    vectorfield;
    struct Pair
    {
        int i, j;
        intfield translations;
    };
    struct Triplet
    {
        int i, j, k;
        int da_j, db_j, dc_j;
        int da_k, db_k, dc_k;
    };
    struct Quadruplet
    {
        int i, j, k, l;
        int da_j, db_j, dc_j;
        int da_k, db_k, dc_k;
        int da_l, db_l, dc_l;
    };
    typedef std::vector<Pair,       managed_allocator<Pair>>       pairfield;
    typedef std::vector<Triplet,    managed_allocator<Triplet>>    tripletfield;
    typedef std::vector<Quadruplet, managed_allocator<Quadruplet>> quadrupletfield;
#else
    typedef std::vector<bool>       boolfield;
    typedef std::vector<int>        intfield;
    typedef std::vector<scalar>     scalarfield;
    typedef std::vector<Vector3>    vectorfield;
    struct Pair
    {
        int i, j;
        std::array<int,3> translations;
    };
    struct Triplet
    {
        int i, j, k;
        int da_j, db_j, dc_j;
        int da_k, db_k, dc_k;
    };
    struct Quadruplet
    {
        int i, j, k, l;
        std::array<int,3> d_j, d_k, d_l;
    };
    typedef std::vector<Pair>       pairfield;
    typedef std::vector<Triplet>    tripletfield;
    typedef std::vector<Quadruplet> quadrupletfield;
#endif

#endif