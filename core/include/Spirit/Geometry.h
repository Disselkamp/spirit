#pragma once
#ifndef INTERFACE_GEOMETRY_H
#define INTERFACE_GEOMETRY_H
#include "DLL_Define_Export.h"
struct State;

#include "Spirit_Defines.h"

// Get Positions of Spins
DLLEXPORT scalar * Geometry_Get_Spin_Positions(State * state, int idx_image=-1, int idx_chain=-1);

// Get Bounds as array (x,y,z)
DLLEXPORT void Geometry_Get_Bounds(State *state, float min[3], float max[3], int idx_image=-1, int idx_chain=-1);
// Get Center as array (x,y,z)
DLLEXPORT void Geometry_Get_Center(State *state, float center[3], int idx_image=-1, int idx_chain=-1);
// Get Cell Bounds as array (x,y,z)
DLLEXPORT void Geometry_Get_Cell_Bounds(State *state, float min[3], float max[3], int idx_image=-1, int idx_chain=-1);

// Get basis vectors ta, tb, tc
DLLEXPORT void Geometry_Get_Basis_Vectors(State *state, float a[3], float b[3], float c[3], int idx_image=-1, int idx_chain=-1);
// Get number of atoms in a basis cell
DLLEXPORT int Geometry_Get_N_Basis_Atoms(State *state, int idx_image=-1, int idx_chain=-1);
// TODO: Get basis atoms
// DLLEXPORT void Geometry_Get_Basis_Atoms(State *state, float ** atoms);

// Get number of basis cells in the three translation directions
DLLEXPORT void Geometry_Get_N_Cells(State *state, int n_cells[3], int idx_image=-1, int idx_chain=-1);
// Get translation vectors ta, tb, tc
DLLEXPORT void Geometry_Get_Translation_Vectors(State *state, float ta[3], float tb[3], float tc[3], int idx_image=-1, int idx_chain=-1);

// Retrieve dimensionality of the system (0, 1, 2, 3)
DLLEXPORT int Geometry_Get_Dimensionality(State * state, int idx_image=-1, int idx_chain=-1);

// Get the 3D Delaunay triangulation. Returns the number of tetrahedrons and
// sets *indices_ptr to point to a list of index 4-tuples.
DLLEXPORT int Geometry_Get_Triangulation(State * state, const int **indices_ptr, int n_cell_step=1, int idx_image=-1, int idx_chain=-1);

#include "DLL_Undefine_Export.h"
#endif
