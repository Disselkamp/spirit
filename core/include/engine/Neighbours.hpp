#pragma once
#ifndef NEIGHBOURS_H
#define NEIGHBOURS_H

#include "Spirit_Defines.h"
#include <data/Geometry.hpp>
#include <engine/Vectormath_Defines.hpp>

#include <vector>
#include <memory>

namespace Engine
{
	namespace Neighbours
	{
		std::vector<scalar> Get_Shell_Radius(const Data::Geometry & geometry, const int n_shells);
		void Pairs_from_Neighbour_Shells(const Data::Geometry & geometry, int nShells, std::vector<int> & shellIndex, pairfield & pairs);
		Vector3 DMI_Normal_from_Pair(const Data::Geometry & geometry, Pair pair, int chirality=1);
		void DDI_Pairs_from_Neighbours(const Data::Geometry & geometry, scalar radius, pairfield & pairs, scalarfield & ddi_magnitude, vectorfield & ddi_normal);

		// void Triplets_from_Neighbours(const Data::Geometry & geometry, tripletfield & triplets);
		// void Quadruplets_from_Neighbours(const Data::Geometry & geometry, tripletfield & quadruplets);
		// void FSC_Quadruplets_from_Neighbours(const Data::Geometry & geometry, tripletfield & quadruplets);


		// // calculates the Bulk DMI vectors
		// void Create_DM_Norm_Vectors_Bulk(const int nos, const vectorfield &spin_pos, const int number_b_vectors,
		// 	const std::vector<Vector3> &boundary_vectors, const int n_shells, const std::vector<std::vector<int>> &n_spins_in_shell,
		// 	const std::vector<std::vector<std::vector<int>>> & neigh, std::vector<std::vector<std::vector<Vector3>>> & neigh_pos,
		// 	const int max_ndm, std::vector<vectorfield> &dm_normal, int dm_chirality=1);

		// // calculates the surface DMI vectors
		// void Create_DM_Norm_Vectors_Surface(const int nos, const vectorfield &spin_pos, const int number_b_vectors,
		// 	const std::vector<Vector3> &boundary_vectors, const int n_shells, const std::vector<std::vector<int>> &n_spins_in_shell,
		// 	const std::vector<std::vector<std::vector<int>>> & neigh, std::vector<std::vector<std::vector<Vector3>>> & neigh_pos,
		// 	const int max_ndm, std::vector<vectorfield> &dm_normal, int dm_chirality=0);

	};// end namespace Neighbours
}// end namespace Engine
#endif