#include <engine/Neighbours.hpp>
#include <engine/Vectormath.hpp>
#include <engine/Manifoldmath.hpp>
#include <utility/IO.hpp>
#include <utility/Logging.hpp>
#include <utility/Exception.hpp>

#include <Eigen/Dense>

#include <numeric>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <algorithm>

using namespace Utility;

namespace Engine
{
	namespace Neighbours
	{
		std::vector<scalar> Get_Shell_Radius(const Data::Geometry & geometry, const int n_shells)
		{
			auto shell_radius = std::vector<scalar>(n_shells);
			
			Vector3 a = geometry.basis[0];
			Vector3 b = geometry.basis[1];
			Vector3 c = geometry.basis[2];

			scalar current_radius=0, dx, min_distance=0;
			int i=0, j=0, k=0;
			int ii, jj, kk;

			// The 15 is a value that is big enough by experience to 
			// produce enough needed shells, but is small enough to run sufficiently fast
			int imax = 15, jmax = 15, kmax = 15;
			Vector3 x0={0,0,0}, x1={0,0,0};

			// Abort condidions for all 3 vectors
			if (a.norm() == 0.0) imax = 0;
			if (b.norm() == 0.0) jmax = 0;
			if (c.norm() == 0.0) kmax = 0;

			for (int n = 0; n < n_shells; ++n)
			{
				current_radius = min_distance;
				min_distance = 1e10;
				for (int iatom = 0; iatom < geometry.n_spins_basic_domain[0]; ++iatom)
				{
					x0 = geometry.basis_atoms[iatom];
					for (ii = imax; ii >= -imax; --ii)
					{
						for (jj = jmax; jj >= -jmax; --jj)
						{
							for (kk = kmax; kk >= -kmax; --kk)
							{
								for (int jatom = 0; jatom < geometry.n_spins_basic_domain[0]; ++jatom)
								{
									if ( !( iatom==jatom && ii==0 && jj==0 && kk==0 ) )
									{
										x1 = geometry.basis_atoms[jatom] + ii*a + jj*b + kk*c;
										dx = (x0-x1).norm();

										if (dx - current_radius > 1e-6 && dx < min_distance)
										{
											min_distance = dx;
											shell_radius[n] = dx;
										}
									}
								}//endfor jatom
							}//endfor kk
						}//endfor jj
					}//endfor ii
				}//endfor iatom
			}
			
			return shell_radius;
		}

		void Pairs_from_Neighbour_Shells(const Data::Geometry & geometry, int nShells, std::vector<int> & shellIndex, pairfield & pairs)
		{
			shellIndex = std::vector<int>(0);
			pairs = pairfield(0);

			auto shell_radius = Get_Shell_Radius(geometry, nShells);
			
			Vector3 a = geometry.basis[0];
			Vector3 b = geometry.basis[1];
			Vector3 c = geometry.basis[2];

			// The nShells + 10 is a value that is big enough by experience to 
			// produce enough needed shells, but is small enough to run sufficiently fast
			int tMax = nShells + 10;
			int imax = tMax, jmax = tMax, kmax = tMax;
			int i,j,k;
			scalar dx, delta, radius;
			Vector3 x0={0,0,0}, x1={0,0,0};

			// Abort condidions for all 3 vectors
			if (a.norm() == 0.0) imax = 0;
			if (b.norm() == 0.0) jmax = 0;
			if (c.norm() == 0.0) kmax = 0;

			for (int iatom = 0; iatom < geometry.n_spins_basic_domain[0]; ++iatom)
			{
				x0 = geometry.basis_atoms[iatom];
				for (int ishell = 0; ishell < nShells; ++ishell)
				{
					radius = shell_radius[ishell];
					for (i = imax; i >= -imax; --i)
					{
						for (j = jmax; j >= -jmax; --j)
						{
							for (k = kmax; k >= -kmax; --k)
							{
								for (int jatom = 0; jatom < geometry.n_spins_basic_domain[0]; ++jatom)
								{
									x1 = geometry.basis_atoms[jatom] + i*a + j*b + k*c;
									dx = (x0-x1).norm();
									delta = std::abs(dx - radius);
									if (delta < 1e-6)
									{
										shellIndex.push_back(ishell);
										pairs.push_back( {iatom, jatom, {i, j, k} } );
									}
								}//endfor jatom
							}//endfor k
						}//endfor j
					}//endfor i
				}//endfor ishell
			}//endfor iatom
		}

		Vector3 DMI_Normal_from_Pair(const Data::Geometry & geometry, Pair pair, int chirality)
		{
			int N = geometry.n_spins_basic_domain[0];
			int Na = geometry.n_cells[0];
			int Nb = geometry.n_cells[1];
			int Nc = geometry.n_cells[2];

			int ispin = 0;
			int jspin = pair.translations[0]*N + pair.translations[1]*N*Na + pair.translations[2]*N*Na*Nb;

			auto ipos = geometry.spin_pos[ispin];
			auto jpos = geometry.spin_pos[jspin];

			if (chirality == 1)
			{
				return (jpos - ipos).normalized();
			}
			else if (chirality == -1)
			{
				return (ipos - jpos).normalized();
			}
			else
			{
				return Vector3{0,0,0};
			}
		}

		void DDI_Pairs_from_Neighbours(const Data::Geometry & geometry, scalar radius, pairfield & pairs)
		{
			auto diagonal = geometry.bounds_max - geometry.bounds_min;
			scalar maxradius = std::min(diagonal[0], std::min(diagonal[1], diagonal[2]));
			// scalar ratio = radius/diagonal.norm();

			// Check for too large DDI radius
			if (radius > maxradius)
			{
				radius = maxradius;
				Log(Log_Level::Warning, Log_Sender::All, "DDI radius is larger than your system! Setting to minimum of system bounds: " + std::to_string(radius));
			}

			if (radius > 1e-6)
			{
				Vector3 a = geometry.basis[0];
				Vector3 b = geometry.basis[1];
				Vector3 c = geometry.basis[2];

				Vector3 ratio{radius/diagonal[0], radius/diagonal[1], radius/diagonal[2]};

				// This should give enough translations to contain all DDI pairs
				int imax = std::min(geometry.n_cells[0], (int)(1.5 * ratio[0] * geometry.n_cells[0]) + 1);
				int jmax = std::min(geometry.n_cells[1], (int)(1.5 * ratio[1] * geometry.n_cells[1]) + 1);
				int kmax = std::min(geometry.n_cells[1], (int)(1.5 * ratio[2] * geometry.n_cells[2]) + 1);

				int i,j,k;
				scalar dx, delta, radius;
				Vector3 x0={0,0,0}, x1={0,0,0};

				// Abort condidions for all 3 vectors
				if (a.norm() == 0.0) imax = 0;
				if (b.norm() == 0.0) jmax = 0;
				if (c.norm() == 0.0) kmax = 0;

				for (int iatom = 0; iatom < geometry.n_spins_basic_domain[0]; ++iatom)
				{
					x0 = geometry.basis_atoms[iatom];
					for (i = imax; i >= -imax; --i)
					{
						for (j = jmax; j >= -jmax; --j)
						{
							for (k = kmax; k >= -kmax; --k)
							{
								for (int jatom = 0; jatom < geometry.n_spins_basic_domain[0]; ++jatom)
								{
									x1 = geometry.basis_atoms[jatom] + i*a + j*b + k*c;
									dx = (x0-x1).norm();
									if (dx < radius)
									{
										pairs.push_back( {iatom, jatom, {i, j, k} } );
									}
								}//endfor jatom
							}//endfor k
						}//endfor j
					}//endfor i
				}//endfor iatom
			}
		}


		// void Triplets_from_Neighbours(const Data::Geometry & geometry, tripletfield & triplets);
		// void Quadruplets_from_Neighbours(const Data::Geometry & geometry, tripletfield & quadruplets);
		// void FSC_Quadruplets_from_Neighbours(const Data::Geometry & geometry, tripletfield & quadruplets);


		// void Triplets_from_Neighbours(const Data::Geometry & geometry, const int n_shells,
		// 	const int nos, const vectorfield &spin_pos, const std::vector<std::vector<int>> &n_spins_in_shell,
		// 	const std::vector<std::vector<std::vector<int>>> & neigh, std::vector<std::vector<std::vector<Vector3>>> & neigh_pos,
		// 	std::vector<std::vector<int>> &segments, std::vector<std::vector<Vector3>> &segments_pos)
		// {
		// 	//========================= Init local vars ================================
		// 	// Retrieve basis vectors from geometry.basis
		// 	Vector3 a = geometry.basis[0];
		// 	Vector3 b = geometry.basis[1];
		// 	Vector3 c = geometry.basis[2];

		// 	Vector3 ipos = { 0, 0, 0 };
		// 	Vector3 jpos = { 0, 0, 0 };
		// 	scalar error = 1.0E-5;
		// 	int ispin, jneigh, jspin, shell, shell_fin;
		// 	Vector3 build_array = { 0, 0, 0 };

		// 	std::string output_to_file = "List of spins in segments\n";
		// 	const int buffer_length = 45;
		// 	output_to_file.reserve(n_shells * buffer_length);	//memory allocation for fast append
		// 	char buffer_string_conversion[buffer_length];
		// 	//------------------------ End Init ----------------------------------------

		// 	if (n_shells >= 2) {
		// 		shell_fin = 2;
		// 	}
		// 	else {
		// 		shell_fin = 1;
		// 	}

		// 	for (ispin = 0; ispin < nos; ++ispin)
		// 	{
		// 		ipos = spin_pos[ispin];
		// 		for (shell = 0; shell < shell_fin; ++shell)
		// 		{
		// 			for (jneigh = 0; jneigh < n_spins_in_shell[ispin][shell]; ++jneigh)
		// 			{
		// 				jspin = neigh[ispin][shell][jneigh];
		// 				jpos = neigh_pos[ispin][shell][jneigh];

		// 				build_array = ipos + a - jpos;
		// 				if (build_array.norm() < error) {
		// 					segments[ispin][0] = jspin;
		// 					segments_pos[ispin][0] = jpos;
		// 				}
		// 				build_array = ipos + b - jpos;
		// 				if (build_array.norm() < error) {
		// 					segments[ispin][1] = jspin;
		// 					segments_pos[ispin][1] = jpos;
		// 				}
		// 				build_array = ipos - a - jpos;
		// 				if (build_array.norm() < error) {
		// 					segments[ispin][2] = jspin;
		// 					segments_pos[ispin][2] = jpos;
		// 				}
		// 				build_array = ipos - b - jpos;
		// 				if (build_array.norm() < error) {
		// 					segments[ispin][3] = jspin;
		// 					segments_pos[ispin][3] = jpos;
		// 				}
		// 			}//endfor jneigh
		// 		}//endfor shell
		// 		snprintf(buffer_string_conversion, buffer_length, "%+06i  %+06i  %+06i  %+06i  %+06i\n", ispin, segments[ispin][0], segments[ispin][1], segments[ispin][2], segments[ispin][3]);
		// 		output_to_file.append(buffer_string_conversion);
		// 	}//endfor ispin
		// 	IO::Dump_to_File(output_to_file, "output/segments.dat");
		// }//end Neighbours::Create_Segments

	};



	// void Neighbours::Create_DM_Norm_Vectors_Bulk(const int nos, const vectorfield &spin_pos, const int number_b_vectors,
	// 	const std::vector<Vector3> &boundary_vectors, const int n_shells, const std::vector<std::vector<int>> &n_spins_in_shell,
	// 	const std::vector<std::vector<std::vector<int>>> & neigh, std::vector<std::vector<std::vector<Vector3>>> & neigh_pos,
	// 	const int max_ndm, std::vector<vectorfield> &dm_normal, int chirality)
	// {
	// 	//========================= Init local vars ================================
	// 	int ispin, jspin, jneigh;
	// 	Vector3 ispin_pos = { 0, 0, 0 };
	// 	Vector3 jspin_pos = { 0, 0, 0 };
	// 	Vector3 r_a = { 0, 0, 0 };
	// 	Vector3 r_b = { 0, 0, 1 };
	// 	Vector3 build_array = { 0, 0, 0 };
	// 	//------------------------ End Init ----------------------------------------

	// 	Log(Log_Level::Debug, Log_Sender::All, "Calculating Bulk DMI norm vectors...");
	// 	for (ispin = 0; ispin < nos; ++ispin)
	// 	{								// loop over all spins
	// 		//Vectormath::Vector_Copy_o(ispin_pos, spin_pos, 0, ispin);
	// 		ispin_pos = spin_pos[ispin];
	// 		for (jneigh = 0; jneigh < n_spins_in_shell[ispin][0]; ++jneigh)
	// 		{	// loop over all neighbours of that spin
	// 			jspin = neigh[ispin][0][jneigh];
	// 			jspin_pos = neigh_pos[ispin][0][jneigh];
	// 			if (chirality == -1)
	// 			{
	// 				r_a = ispin_pos - jspin_pos; //get DMI vec with chirality "-"
	// 			}
	// 			else if (chirality == 2)
	// 			{
	// 				// This should only be used in 2D case, not in bulk
	// 				r_a = (jspin_pos - ispin_pos).cross(Vector3{ 0,0,1 });
	// 			}
	// 			else if (chirality == -2)
	// 			{
	// 				// This should only be used in 2D case, not in bulk
	// 				r_a = (ispin_pos - jspin_pos).cross(Vector3{ 0,0,1 });
	// 			}
	// 			else
	// 			{
	// 				r_a = jspin_pos - ispin_pos; // get DMI vec with chirality "+"
	// 			}
	// 			r_a.normalize();
	// 			dm_normal[ispin][jneigh] = r_a;
	// 		}//endfor jneigh
	// 	}//endfor ispin
	// 	Log(Log_Level::Debug, Log_Sender::All, "Done calculating Bulk DMI norm vectors.");
	// }//end Neighbours::Create_DM_Norm_Vectors_Bulk

	// void Neighbours::Create_DM_Norm_Vectors_Surface(const int nos, const vectorfield &spin_pos, const int number_b_vectors,
	// 	const std::vector<Vector3> &boundary_vectors, const int n_shells, const std::vector<std::vector<int>> &n_spins_in_shell,
	// 	const std::vector<std::vector<std::vector<int>>> & neigh, std::vector<std::vector<std::vector<Vector3>>> & neigh_pos,
	// 	const int max_ndm, std::vector<vectorfield> &dm_normal, int chirality)
	// {
	// 	//========================= Init local vars ================================
	// 	int ispin, jneigh;
	// 	Vector3 unit_vec_z = { 0, 0, 1 };
	// 	Vector3 build_array_1 = { 0, 0, 0 };
	// 	Vector3 build_array_2 = { 0, 0, 0 };
	// 	//------------------------ End Init ----------------------------------------
	// 	Create_DM_Norm_Vectors_Bulk(nos, spin_pos, number_b_vectors, boundary_vectors, n_shells, n_spins_in_shell, neigh, neigh_pos, max_ndm, dm_normal);

	// 	for (ispin = 0; ispin < nos; ++ispin)
	// 	{
	// 		for (jneigh = 0; jneigh < max_ndm; ++jneigh)
	// 		{
	// 			if (chirality == 2)
	// 			{
	// 				// This should only be used in 2D case, not in bulk
	// 				dm_normal[ispin][jneigh] = dm_normal[ispin][jneigh].cross(Vector3{ 0,0,1 });
	// 			}
	// 			else if (chirality == -2)
	// 			{
	// 				// This should only be used in 2D case, not in bulk
	// 				dm_normal[ispin][jneigh] = -dm_normal[ispin][jneigh].cross(Vector3{ 0,0,1 });
	// 			}
	// 			else if (chirality == -1)
	// 			{
	// 				dm_normal[ispin][jneigh] *= -1;
	// 			}
	// 		}
	// 	}
	// 	Log(Log_Level::Debug, Log_Sender::All, "Done calculating Surface DMI norm vectors.");
	// }//end Neighbours::Create_DM_Norm_Vectors_Surface


}// end Namespace Engine