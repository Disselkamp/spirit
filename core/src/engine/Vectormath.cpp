#ifndef USE_CUDA

#include <engine/Neighbours.hpp>
#include <engine/Vectormath.hpp>
#include <engine/Manifoldmath.hpp>
#include <utility/Logging.hpp>
#include <utility/Exception.hpp>

#include <Eigen/Dense>

#include <array>
#include <algorithm>

namespace Engine
{
	namespace Vectormath
	{

		void rotate(const Vector3 & v, const Vector3 & axis, const scalar & angle, Vector3 & v_out)
		{
			v_out = v * std::cos(angle) + axis.cross(v) * std::sin(angle);
		}

		Vector3 decompose(const Vector3 & v, const std::vector<Vector3> & basis)
		{
			Eigen::Ref<const Matrix3> A = Eigen::Map<const Matrix3>(basis[0].data());
			return A.colPivHouseholderQr().solve(v);
		}

		/////////////////////////////////////////////////////////////////
		

		void Build_Spins(vectorfield & spin_pos, const std::vector<Vector3> & basis_atoms, const std::vector<Vector3> & translation_vectors, const std::vector<int> & n_cells)
		{
			// Check for erronous input placing two spins on the same location
			Vector3 sp;
			for (unsigned int i = 0; i < basis_atoms.size(); ++i)
			{
				for (unsigned int j = 0; j < basis_atoms.size(); ++j)
				{
					for (int k1 = -2; k1 <= 2; ++k1)
					{
						for (int k2 = -2; k2 <= 2; ++k2)
						{
							for (int k3 = -2; k3 <= 2; ++k3)
							{ 
								// Norm is zero if translated basis atom is at position of another basis atom
								sp = basis_atoms[i] - (basis_atoms[j]
									+ k1*translation_vectors[0] + k2*translation_vectors[1] + k3*translation_vectors[2]);
								if ((i != j || k1 != 0 || k2 != 0 || k3 != 0) && std::abs(sp[0]) < 1e-9 && std::abs(sp[1]) < 1e-9 && std::abs(sp[2]) < 1e-9)
								{
									Log(Utility::Log_Level::Severe, Utility::Log_Sender::All, "Unable to initialize Spin-System, since 2 spins occupy the same space.\nPlease check the config file!");
									Log.Append_to_File();
									throw Utility::Exception::System_not_Initialized;
								}
							}
						}
					}
				}
			}

			// Build up the spins array
			int i, j, k, s, pos;
			int nos_basic = basis_atoms.size();
			//int nos = nos_basic * n_cells[0] * n_cells[1] * n_cells[2];
			Vector3 build_array;
			for (k = 0; k < n_cells[2]; ++k) {
				for (j = 0; j < n_cells[1]; ++j) {
					for (i = 0; i < n_cells[0]; ++i) {
						for (s = 0; s < nos_basic; ++s) {
							pos = k*n_cells[1] * n_cells[0] * nos_basic + j*n_cells[0] * nos_basic + i*nos_basic + s;
							build_array = i*translation_vectors[0] + j*translation_vectors[1] + k*translation_vectors[2];
							// paste initial spin orientations across the lattice translations
							//spins[dim*nos + pos] = spins[dim*nos + s];
							// calculate the spin positions
							spin_pos[pos] = basis_atoms[s] + build_array;
						}// endfor s
					}// endfor k
				}// endfor j
			}// endfor dim

		};// end Build_Spins


		std::array<scalar,3> Magnetization(const vectorfield & vf)
		{
			std::array<scalar, 3> M{0, 0, 0};
			int nos = vf.size();
			scalar scale = 1/(scalar)nos;
			for (int i=0; i<nos; ++i)
			{
				M[0] += vf[i][0]*scale;
				M[1] += vf[i][1]*scale;
				M[2] += vf[i][2]*scale;
			}
			return M;
		}

		scalar TopologicalCharge(const vectorfield & vf)
		{
        	Log(Utility::Log_Level::Warning, Utility::Log_Sender::All, std::string("Calculating the topological charge is not yet implemented"));
			return 0;
		}

		// Utility function for the SIB Optimizer
		void transform(const vectorfield & spins, const vectorfield & force, vectorfield & out)
		{
			Vector3 e1, a2, A;
			scalar detAi;
			for (unsigned int i = 0; i < spins.size(); ++i)
			{
				e1 = spins[i];
				A = force[i];

				// 1/determinant(A)
				detAi = 1.0 / (1 + pow(A.norm(), 2.0));

				// calculate equation without the predictor?
				a2 = e1 + e1.cross(A);

				out[i][0] = (a2[0] * (1 + A[0] * A[0])    + a2[1] * (A[0] * A[1] + A[2]) + a2[2] * (A[0] * A[2] - A[1]))*detAi;
				out[i][1] = (a2[0] * (A[1] * A[0] - A[2]) + a2[1] * (1 + A[1] * A[1])    + a2[2] * (A[1] * A[2] + A[0]))*detAi;
				out[i][2] = (a2[0] * (A[2] * A[0] + A[1]) + a2[1] * (A[2] * A[1] - A[0]) + a2[2] * (1 + A[2] * A[2]))*detAi;
			}
		}
		void get_random_vectorfield(const Data::Spin_System & sys, scalar epsilon, vectorfield & xi)
		{
			for (int i = 0; i < sys.nos; ++i)
			{
				for (int dim = 0; dim < 3; ++dim)
				{
					// PRNG gives RN int [0,1] -> [-1,1] -> multiply with epsilon
					xi[i][dim] = epsilon*(sys.llg_parameters->distribution_int(sys.llg_parameters->prng) * 2 - 1); 
				}
			}
		}



		/////////////////////////////////////////////////////////////////

		void assign(scalarfield & sf_dest, const scalarfield& sf_source)
		{
			sf_dest = sf_source;
		}

		void fill(scalarfield & sf, scalar s)
		{
			for (unsigned int i = 0; i<sf.size(); ++i)
			{
				sf[i] = s;
			}
		}

		void scale(scalarfield & sf, scalar s)
		{
			for (unsigned int i = 0; i<sf.size(); ++i)
			{
				sf[i] *= s;
			}
		}

		scalar sum(const scalarfield & sf)
		{
			scalar ret = 0;
			for (unsigned int i = 0; i<sf.size(); ++i)
			{
				ret += sf[i];
			}
			return ret;
		}

		scalar mean(const scalarfield & sf)
		{
			scalar ret = 0;
			for (unsigned int i = 0; i<sf.size(); ++i)
			{
				ret = (i - 1) / i * ret + sf[i] / i;
			}
			return ret;
		}

		void fill(vectorfield & vf, const Vector3 & v)
		{
			for (unsigned int i=0; i<vf.size(); ++i)
			{
				vf[i] = v;
			}
		}

		void normalize_vectors(vectorfield & vf)
		{
			for (unsigned int i=0; i<vf.size(); ++i)
			{
				vf[i].normalize();
			}
		}
		
		std::pair<scalar, scalar> minmax_component(const vectorfield & v1)
		{
			scalar min=1e6, max=-1e6;
			std::pair<scalar, scalar> minmax;
			for (unsigned int i = 0; i < v1.size(); ++i)
			{
				for (int dim = 0; dim < 3; ++dim)
				{
					if (v1[i][dim] < min) min = v1[i][dim];
					if (v1[i][dim] > max) max = v1[i][dim];
				}
			}
			minmax.first = min;
			minmax.second = max;
			return minmax;
		}
    	scalar  max_abs_component(const vectorfield & vf)
		{
			// We want the Maximum of Absolute Values of all force components on all images
			scalar absmax = 0;
			// Find minimum and maximum values
			std::pair<scalar,scalar> minmax = minmax_component(vf);
			// Mamimum of absolute values
			absmax = std::max(absmax, std::abs(minmax.first));
			absmax = std::max(absmax, std::abs(minmax.second));
			// Return
			return absmax;
		}

		void scale(vectorfield & vf, const scalar & sc)
		{
			for (unsigned int i=0; i<vf.size(); ++i)
			{
				vf[i] *= sc;
			}
		}

		Vector3 sum(const vectorfield & vf)
		{
			Vector3 ret = { 0,0,0 };
			for (unsigned int i = 0; i<vf.size(); ++i)
			{
				ret += vf[i];
			}
			return ret;
		}

		Vector3 mean(const vectorfield & vf)
		{
			Vector3 ret = { 0,0,0 };
			for (unsigned int i = 0; i<vf.size(); ++i)
			{
				ret = (i-1)/i * ret + vf[i]/i;
			}
			return ret;
		}



		// computes the inner product of two vectorfields v1 and v2
		scalar dot(const vectorfield & v1, const vectorfield & v2)
		{
			scalar x = 0;
			for (unsigned int i = 0; i<v1.size(); ++i)
			{
				x += v1[i].dot(v2[i]);
			}
			return x;
		}

		// computes the inner products of vectors in v1 and v2
		// v1 and v2 are vectorfields
		void dot(const vectorfield & v1, const vectorfield & v2, scalarfield & out)
		{
			for (unsigned int i=0; i<v1.size(); ++i)
			{
				out[i] = v1[i].dot(v2[i]);
			}
		}

		// computes the vector (cross) products of vectors in v1 and v2
		// v1 and v2 are vector fields
		void cross(const vectorfield & v1, const vectorfield & v2, vectorfield & out)
		{
			for (unsigned int i=0; i<v1.size(); ++i)
			{
				out[i] = v1[i].cross(v2[i]);
			}
		}



		// out[i] += c*a
		void add_c_a(const scalar & c, const Vector3 & a, vectorfield & out)
		{
			for(unsigned int idx = 0; idx < out.size(); ++idx)
			{
				out[idx] += c*a;
			}
		}
		// out[i] += c*a[i]
		void add_c_a(const scalar & c, const vectorfield & a, vectorfield & out)
		{
			for(unsigned int idx = 0; idx < out.size(); ++idx)
			{
				out[idx] += c*a[idx];
			}
		}
		
		// out[i] = c*a
		void set_c_a(const scalar & c, const Vector3 & a, vectorfield & out)
		{
			for(unsigned int idx = 0; idx < out.size(); ++idx)
			{
				out[idx] = c*a;
			}
		}
		// out[i] = c*a[i]
		void set_c_a(const scalar & c, const vectorfield & a, vectorfield & out)
		{
			for(unsigned int idx = 0; idx < out.size(); ++idx)
			{
				out[idx] = c*a[idx];
			}
		}


		// out[i] += c * a*b[i]
		void add_c_dot(const scalar & c, const Vector3 & a, const vectorfield & b, scalarfield & out)
		{
			for(unsigned int idx = 0; idx < out.size(); ++idx)
			{
				out[idx] += c*a.dot(b[idx]);
			}
		}
		// out[i] += c * a[i]*b[i]
		void add_c_dot(const scalar & c, const vectorfield & a, const vectorfield & b, scalarfield & out)
		{
			for(unsigned int idx = 0; idx < out.size(); ++idx)
			{
				out[idx] += c*a[idx].dot(b[idx]);
			}
		}

		// out[i] = c * a*b[i]
		void set_c_dot(const scalar & c, const Vector3 & a, const vectorfield & b, scalarfield & out)
		{
			for(unsigned int idx = 0; idx < out.size(); ++idx)
			{
				out[idx] = c*a.dot(b[idx]);
			}
		}
		// out[i] = c * a[i]*b[i]
		void set_c_dot(const scalar & c, const vectorfield & a, const vectorfield & b, scalarfield & out)
		{
			for(unsigned int idx = 0; idx < out.size(); ++idx)
			{
				out[idx] = c*a[idx].dot(b[idx]);
			}
		}


		// out[i] += c * a x b[i]
		void add_c_cross(const scalar & c, const Vector3 & a, const vectorfield & b, vectorfield & out)
		{
			for(unsigned int idx = 0; idx < out.size(); ++idx)
			{
				out[idx] += c*a.cross(b[idx]);
			}
		}
		// out[i] += c * a[i] x b[i]
		void add_c_cross(const scalar & c, const vectorfield & a, const vectorfield & b, vectorfield & out)
		{
			for(unsigned int idx = 0; idx < out.size(); ++idx)
			{
				out[idx] += c*a[idx].cross(b[idx]);
			}
		}
		
		// out[i] = c * a x b[i]
		void set_c_cross(const scalar & c, const Vector3 & a, const vectorfield & b, vectorfield & out)
		{
			for(unsigned int idx = 0; idx < out.size(); ++idx)
			{
				out[idx] = c*a.cross(b[idx]);
			}
		}
		// out[i] = c * a[i] x b[i]
		void set_c_cross(const scalar & c, const vectorfield & a, const vectorfield & b, vectorfield & out)
		{
			for(unsigned int idx = 0; idx < out.size(); ++idx)
			{
				out[idx] = c*a[idx].cross(b[idx]);
			}
		}
		// SST -- gradient (j_e*grad)*S
		void gradient(const vectorfield & spins, const Engine::Hamiltonian & hamiltonian, const Data::Geometry & geometry, const Vector3 & je, vectorfield & s_c_grad)
		{
			// std::cout << "start gradient" << std::endl;
			vectorfield translations = { { 0,0,0 }, { 0,0,0 }, { 0,0,0 } };
			std::vector <int> n_cells = geometry.n_cells;

			std::vector<bool> boundary_conditions = hamiltonian.boundary_conditions;

			Vector3 a = geometry.translation_vectors[0];
			Vector3 b = geometry.translation_vectors[1];
			Vector3 c = geometry.translation_vectors[2];

			Neighbourfield neigh;
			Neighbours::get_Neighbours(geometry, neigh);
			// std::cout << "Neighbours: " << neigh[k][j].jatom << " " << neigh[k][j].translations[0]

			Vector3  diffq, diffqx, diffqy, diffqz;

			for(unsigned int i = 0; i < spins.size(); ++i)
			{
				auto translations_i = translations_from_idx(n_cells, geometry.n_spins_basic_domain, i);
				int k = i%geometry.n_spins_basic_domain;
				double n = 0;
				for(unsigned int j = 0; j < neigh[k].size(); ++j)
				{
					std::array<int,3> neightranslation;
					neightranslation[0] = neigh[k][j].translations[0]; neightranslation[1] = neigh[k][j].translations[1]; neightranslation[2] = neigh[k][j].translations[2];
					if ( boundary_conditions_fulfilled(geometry.n_cells, boundary_conditions, translations_i, neightranslation) )
					{
						Vector3 translationVec3 = geometry.basis_atoms[neigh[k][j].jatom] - geometry.basis_atoms[k] + neigh[k][j].translations[0]*a + neigh[k][j].translations[1]*b + neigh[k][j].translations[2]*c; // -geometry.basis_atoms[k]
						// std::array<int,3> translation;
						// translation[0] = (int)translationVec3[0]; translation[1] = (int)translationVec3[1]; translation[2] = (int)translationVec3[2];
						int idx = idx_from_translations(n_cells, geometry.n_spins_basic_domain, translations_i, neightranslation);
						
						// std::cout << "i: " << i << "  j: " << j << "  idx: " << idx << std::endl;
						// std::cout << "neighbours   " << " : ";
						// std::cout << neigh[k][j].jatom << " ," << neigh[k][j].translations[0] << " ," << neigh[k][j].translations[1] << " ," << neigh[k][j].translations[2] << std::endl; 
						// std::cout << "translations " << i << j << " : ";
						// std::cout << neightranslation[0] << " ," << neightranslation[1] << " ," << neightranslation[2] << std::endl;

						diffq = (spins[idx]-spins[i])/translationVec3.norm();

						// std::cout << "Spin[i]: " << spins[i][0] << " " << spins[i][1] << " " << spins[i][2] 
						// << " Spin[idx]: " << spins[idx][0] << " " << spins[idx][1] << " " << spins[idx][2] << std::endl;

						// std::cout << "translvec: " << translationVec3[0] << std::endl;

						diffqx += translationVec3[0]*diffq;
						diffqy += translationVec3[1]*diffq;
						diffqz += translationVec3[2]*diffq;
						n += 1.0; // 1 fuer Aussen-Spins entsprechend boundary conditions sonst 2
					}
				}

				diffqx = diffqx/n; diffqy = diffqy/n; diffqz = diffqz/n;
				s_c_grad.push_back(je[0]*diffqx+je[1]*diffqy+je[2]*diffqz); // dot(je, diffqxyz, scalarfield & out)
				
				diffqx = { 0,0,0 }; diffqy = { 0,0,0 }; diffqz = { 0,0,0 };
			}
		}

		inline int idx_from_translations(const intfield & n_cells, const int & n_spins_basic_domain, const std::array<int,3> & translations_i, const std::array<int,3> & translations)
		{
			int Na = n_cells[0];
			int Nb = n_cells[1];
			int Nc = n_cells[2];
			int N  = n_spins_basic_domain;
			
			int da = translations_i[0]+translations[0];
			int db = translations_i[1]+translations[1];
			int dc = translations_i[2]+translations[2];

			if (translations[0] < 0)
				da += N*Na;
			if (translations[1] < 0)
				db += N*Na*Nb;
			if (translations[2] < 0)
				dc += N*Na*Nb*Nc;
				
			int idx = (da%Na)*N + (db%Nb)*N*Na + (dc%Nc)*N*Na*Nb;
			
			return idx;
		}

		inline std::array<int,3> translations_from_idx(const intfield & n_cells, const int & n_spins_basic_domain, int idx)
		{
			std::array<int,3> ret;
			int Na = n_cells[0];
			int Nb = n_cells[1];
			int Nc = n_cells[2];
			int N  = n_spins_basic_domain;
			ret[2] = idx/(Na*Nb);
			ret[1] = (idx-ret[2]*Na*Nb)/Na;
			ret[0] = idx-ret[2]*Na*Nb-ret[1]*Na;
			return ret;
		}

		inline bool boundary_conditions_fulfilled(const intfield & n_cells, const std::vector<bool> & boundary_conditions, const std::array<int, 3> & translations_i, const std::array<int, 3> & translations_j)
		{
			int da = translations_i[0] + translations_j[0];
			int db = translations_i[1] + translations_j[1];
			int dc = translations_i[2] + translations_j[2];
			return ((boundary_conditions[0] || (0 <= da && da < n_cells[0])) &&
				(boundary_conditions[1] || (0 <= db && db < n_cells[1])) &&
				(boundary_conditions[2] || (0 <= dc && dc < n_cells[2])));
		}

		// index from translationvector and vise versa
		// inline int idx_from_translations(const std::vector<int> & n_cells, int n_spins_basic_domain, const Vector3 & translations_i, const Vector3 & translations)
		// {
		// 	int Na = n_cells[0];
		// 	int Nb = n_cells[1];
		// 	int Nc = n_cells[2];
		// 	int N  = n_spins_basic_domain;

		// 	int da = translations_i[0]+translations[0];
		// 	int db = translations_i[1]+translations[1];
		// 	int dc = translations_i[2]+translations[2];

		// 	if (translations[0] < 0)
		// 		da += N*Na;
		// 	if (translations[1] < 0)
		// 		db += N*Na*Nb;
		// 	if (translations[2] < 0)
		// 		dc += N*Na*Nb*Nc;
				
		// 	int idx = (da%Na)*N + (db%Nb)*N*Na + (dc%Nc)*N*Na*Nb;
			
		// 	return idx;
		// }
		// inline intfield translations_from_idx(const std::vector<int> & n_cells, int n_spins_basic_domain, int idx)
		// {
		// 	// spin_pos[idx]
		// 	std::cout << "idx: " << idx << std::endl;

		// 	intfield ret;
		// 	int Na = n_cells[0];
		// 	int Nb = n_cells[1];
		// 	int Nc = n_cells[2];
		// 	int N  = n_spins_basic_domain;
		// 	ret[2] = idx/(Na*Nb);
		// 	ret[1] = (idx-ret[2]*Na*Nb)/Na;
		// 	ret[0] = idx-ret[2]*Na*Nb-ret[1]*Na;
		
		// 	std::cout << "ret: " << ret[0] << " " << ret[1] << " " << ret[2] << std::endl;
		
		// 	return ret;
		// }
	}
}

#endif