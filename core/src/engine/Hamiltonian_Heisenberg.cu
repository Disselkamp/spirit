#ifdef USE_CUDA

#define _USE_MATH_DEFINES
#include <cmath>

#include <Eigen/Dense>

#include <engine/Hamiltonian_Heisenberg.hpp>
#include <engine/Vectormath.hpp>
#include <data/Spin_System.hpp>
#include <utility/Constants.hpp>

using std::vector;
using std::function;

using namespace Data;
using namespace Utility;

namespace Engine
{
	Hamiltonian_Heisenberg::Hamiltonian_Heisenberg(
			scalarfield mu_s,
			intfield external_field_index, scalarfield external_field_magnitude, vectorfield external_field_normal,
			intfield anisotropy_index, scalarfield anisotropy_magnitude, vectorfield anisotropy_normal,
			std::vector<indexPairs> Exchange_indices, std::vector<scalarfield> Exchange_magnitude,
			std::vector<indexPairs> DMI_indices, std::vector<scalarfield> DMI_magnitude, std::vector<vectorfield> DMI_normal,
			std::vector<indexPairs> DD_indices, std::vector<scalarfield> DD_magnitude, std::vector<vectorfield> DD_normal,
			std::vector<indexQuadruplets> quadruplet_indices, std::vector<scalarfield> quadruplet_magnitude,
			std::vector<bool> boundary_conditions
	) :
		Hamiltonian(boundary_conditions),
		mu_s(mu_s),
		external_field_index(external_field_index), external_field_magnitude(external_field_magnitude), external_field_normal(external_field_normal),
		anisotropy_index(anisotropy_index), anisotropy_magnitude(anisotropy_magnitude), anisotropy_normal(anisotropy_normal),
		Exchange_indices(Exchange_indices), Exchange_magnitude(Exchange_magnitude),
		DMI_indices(DMI_indices), DMI_magnitude(DMI_magnitude), DMI_normal(DMI_normal),
		DD_indices(DD_indices), DD_magnitude(DD_magnitude), DD_normal(DD_normal),
		Quadruplet_indices(quadruplet_indices), Quadruplet_magnitude(quadruplet_magnitude)
	{
		// Renormalize the external field from Tesla to whatever
		for (unsigned int i = 0; i < external_field_magnitude.size(); ++i)
		{
			this->external_field_magnitude[i] = this->external_field_magnitude[i] * Constants::mu_B * mu_s[i];
		}

		this->Update_Energy_Contributions();
	}

	void Hamiltonian_Heisenberg::Update_Energy_Contributions()
	{
		this->energy_contributions_per_spin = std::vector<std::pair<std::string, scalarfield>>(0);

		// External field
		if (this->external_field_index.size() > 0)
		{
			this->energy_contributions_per_spin.push_back({"Zeeman", scalarfield(0)});
			this->idx_zeeman = this->energy_contributions_per_spin.size()-1;
		}
		else this->idx_zeeman = -1;
		// Anisotropy
		if (this->anisotropy_index.size() > 0)
		{
			this->energy_contributions_per_spin.push_back({"Anisotropy", scalarfield(0) });
			this->idx_anisotropy = this->energy_contributions_per_spin.size()-1;
		}
		else this->idx_anisotropy = -1;
		// Exchange
		if (this->Exchange_indices[0].size() > 0)
		{
			this->energy_contributions_per_spin.push_back({"Exchange", scalarfield(0) });
			this->idx_exchange = this->energy_contributions_per_spin.size()-1;
		}
		else this->idx_exchange = -1;
		// DMI
		if (this->DMI_indices[0].size() > 0)
		{
			this->energy_contributions_per_spin.push_back({"DMI", scalarfield(0) });
			this->idx_dmi = this->energy_contributions_per_spin.size()-1;
		}
		else this->idx_dmi = -1;
		// Dipole-Dipole
		if (this->DD_indices[0].size() > 0)
		{
			this->energy_contributions_per_spin.push_back({"DD", scalarfield(0) });
			this->idx_dd = this->energy_contributions_per_spin.size()-1;
		}
		else this->idx_dd = -1;
		// Quadruplet
		if (this->Quadruplet_indices[0].size() > 0)
		{
			this->energy_contributions_per_spin.push_back({"Quadruplet", scalarfield(0) });
			this->idx_quadruplet = this->energy_contributions_per_spin.size()-1;
		}
		else this->idx_quadruplet = -1;
	}


	void Hamiltonian_Heisenberg::Energy_Contributions_per_Spin(const vectorfield & spins, std::vector<std::pair<std::string, scalarfield>> & contributions)
	{
		int nos = spins.size();
		for (auto& pair : contributions)
		{
			// Allocate if not already allocated
			if (pair.second.size() != nos) pair.second = scalarfield(nos, 0);
			// Otherwise set to zero
			else for (auto& pair : contributions) Vectormath::fill(pair.second, 0);
		}
		

		// External field
		if (this->idx_zeeman >=0 ) E_Zeeman(spins, contributions[idx_zeeman].second);

		// Anisotropy
		if (this->idx_anisotropy >=0 ) E_Anisotropy(spins, contributions[idx_anisotropy].second);

		// Pairs
		//		Loop over periodicity
		for (int i_periodicity = 0; i_periodicity < 8; ++i_periodicity)
		{
			// Check if boundary conditions contain this periodicity
			if ((i_periodicity == 0)
				|| (i_periodicity == 1 && this->boundary_conditions[0])
				|| (i_periodicity == 2 && this->boundary_conditions[1])
				|| (i_periodicity == 3 && this->boundary_conditions[2])
				|| (i_periodicity == 4 && this->boundary_conditions[0] && this->boundary_conditions[1])
				|| (i_periodicity == 5 && this->boundary_conditions[0] && this->boundary_conditions[2])
				|| (i_periodicity == 6 && this->boundary_conditions[1] && this->boundary_conditions[2])
				|| (i_periodicity == 7 && this->boundary_conditions[0] && this->boundary_conditions[1] && this->boundary_conditions[2]))
			{
				//		Energies of this periodicity
				// Exchange
				if (this->idx_exchange >=0 ) E_Exchange(spins, Exchange_indices[i_periodicity], Exchange_magnitude[i_periodicity], contributions[idx_exchange].second);
				// DMI
				if (this->idx_dmi >=0 ) E_DMI(spins, DMI_indices[i_periodicity], DMI_magnitude[i_periodicity], DMI_normal[i_periodicity], contributions[idx_dmi].second);
				// DD
				if (this->idx_dd >=0 ) E_DD(spins, DD_indices[i_periodicity], DD_magnitude[i_periodicity], DD_normal[i_periodicity], contributions[idx_dd].second);
				// Quadruplet
				if (this->idx_quadruplet >=0 ) E_Quadruplet(spins, Quadruplet_indices[i_periodicity], Quadruplet_magnitude[i_periodicity], contributions[idx_quadruplet].second);
			}
		}
		
		cudaDeviceSynchronize();

		// Return
		//return this->E;
	}

	
	__global__ void CU_E_Zeeman(const Vector3 * spins, const int * external_field_index, const scalar * external_field_magnitude, const Vector3 * external_field_normal, scalar * Energy, size_t size)
	{
		for(auto idx = blockIdx.x * blockDim.x + threadIdx.x;
			idx < size;
			idx +=  blockDim.x * gridDim.x)
		{
			atomicAdd(&Energy[external_field_index[idx]], - external_field_magnitude[idx] * external_field_normal[idx].dot(spins[external_field_index[idx]]));
		}
	}
	void Hamiltonian_Heisenberg::E_Zeeman(const vectorfield & spins, scalarfield & Energy)
	{
		int size = this->external_field_index.size();
		CU_E_Zeeman<<<(size+1023)/1024, 1024>>>(spins.data(), this->external_field_index.data(), this->external_field_magnitude.data(), this->external_field_normal.data(), Energy.data(), size);
	}


	__global__ void CU_E_Anisotropy(const Vector3 * spins, const int * anisotropy_index, const scalar * anisotropy_magnitude, const Vector3 * anisotropy_normal, scalar * Energy, size_t size)
	{
		for(auto idx = blockIdx.x * blockDim.x + threadIdx.x;
			idx < size;
			idx +=  blockDim.x * gridDim.x)
		{
			atomicAdd(&Energy[anisotropy_index[idx]], - anisotropy_magnitude[idx] * std::pow(anisotropy_normal[idx].dot(spins[anisotropy_index[idx]]), 2.0));
		}
	}
	void Hamiltonian_Heisenberg::E_Anisotropy(const vectorfield & spins, scalarfield & Energy)
	{
		int size = this->anisotropy_index.size();
		CU_E_Anisotropy<<<(size+1023)/1024, 1024>>>(spins.data(), this->anisotropy_index.data(), this->anisotropy_magnitude.data(), this->anisotropy_normal.data(), Energy.data(), size);
	}


	__global__ void CU_E_Exchange(const Vector3 * spins, const indexPair * pairs, const scalar * J_ij, scalar * Energy, size_t size)
	{
		for(auto iPair = blockIdx.x * blockDim.x + threadIdx.x;
			iPair < size;
			iPair +=  blockDim.x * gridDim.x)
		{
			int ispin = pairs[iPair][0];
			int jspin = pairs[iPair][1];
			scalar sc = - 0.5 * J_ij[iPair] * spins[ispin].dot(spins[jspin]);
			atomicAdd(&Energy[ispin], sc);
			atomicAdd(&Energy[jspin], sc);
		}
	}
	void Hamiltonian_Heisenberg::E_Exchange(const vectorfield & spins, indexPairs & indices, scalarfield & J_ij, scalarfield & Energy)
	{
		int size = indices.size();
		CU_E_Exchange<<<(size+1023)/1024, 1024>>>(spins.data(), indices.data(), J_ij.data(), Energy.data(), size);
	}


	__global__ void CU_E_DMI(const Vector3 * spins, const indexPair * pairs, const scalar * DMI_magnitude, const Vector3 * DMI_normal, scalar * Energy, size_t size)
	{
		for(auto iPair = blockIdx.x * blockDim.x + threadIdx.x;
			iPair < size;
			iPair +=  blockDim.x * gridDim.x)
		{
			int ispin = pairs[iPair][0];
			int jspin = pairs[iPair][1];
			scalar sc = - 0.5 *  DMI_magnitude[iPair] * DMI_normal[iPair].dot(spins[ispin].cross(spins[jspin]));
			atomicAdd(&Energy[ispin], sc);
			atomicAdd(&Energy[jspin], sc);
		}
	}
	void Hamiltonian_Heisenberg::E_DMI(const vectorfield & spins, indexPairs & indices, scalarfield & DMI_magnitude, vectorfield & DMI_normal, scalarfield & Energy)
	{
		int size = indices.size();
		CU_E_DMI<<<(size+1023)/1024, 1024>>>(spins.data(), indices.data(), DMI_magnitude.data(), DMI_normal.data(), Energy.data(), size);
	}


	void Hamiltonian_Heisenberg::E_DD(const vectorfield & spins, indexPairs & indices, scalarfield & DD_magnitude, vectorfield & DD_normal, scalarfield & Energy)
	{
		//scalar mult = -Utility::Constants::mu_B*Utility::Constants::mu_B*1.0 / 4.0 / M_PI; // multiply with mu_B^2
		scalar mult = 0.5*0.0536814951168; // mu_0*mu_B**2/(4pi*10**-30) -- the translations are in angstr�m, so the |r|[m] becomes |r|[m]*10^-10

		for (unsigned int i_pair = 0; i_pair < indices.size(); ++i_pair)
		{
			if (DD_magnitude[i_pair] > 0.0)
			{
				Energy[indices[i_pair][0]] -= mult / std::pow(DD_magnitude[i_pair], 3.0) *
					(3 * spins[indices[i_pair][1]].dot(DD_normal[i_pair]) * spins[indices[i_pair][0]].dot(DD_normal[i_pair]) - spins[indices[i_pair][0]].dot(spins[indices[i_pair][1]]));
				Energy[indices[i_pair][1]] -= mult / std::pow(DD_magnitude[i_pair], 3.0) *
					(3 * spins[indices[i_pair][1]].dot(DD_normal[i_pair]) * spins[indices[i_pair][0]].dot(DD_normal[i_pair]) - spins[indices[i_pair][0]].dot(spins[indices[i_pair][1]]));
			}

		}
	}// end DipoleDipole


	void Hamiltonian_Heisenberg::E_Quadruplet(const vectorfield & spins, indexQuadruplets & indices, scalarfield & magnitude, scalarfield & Energy)
	{
		for (unsigned int i_pair = 0; i_pair < indices.size(); ++i_pair)
		{
			Energy[indices[i_pair][0]] -= 0.25*magnitude[i_pair] * (spins[indices[i_pair][0]].dot(spins[indices[i_pair][1]])) * (spins[indices[i_pair][2]].dot(spins[indices[i_pair][3]]));
			Energy[indices[i_pair][1]] -= 0.25*magnitude[i_pair] * (spins[indices[i_pair][0]].dot(spins[indices[i_pair][1]])) * (spins[indices[i_pair][2]].dot(spins[indices[i_pair][3]]));
			Energy[indices[i_pair][2]] -= 0.25*magnitude[i_pair] * (spins[indices[i_pair][0]].dot(spins[indices[i_pair][1]])) * (spins[indices[i_pair][2]].dot(spins[indices[i_pair][3]]));
			Energy[indices[i_pair][3]] -= 0.25*magnitude[i_pair] * (spins[indices[i_pair][0]].dot(spins[indices[i_pair][1]])) * (spins[indices[i_pair][2]].dot(spins[indices[i_pair][3]]));
		}
	}



	void Hamiltonian_Heisenberg::Gradient(const vectorfield & spins, vectorfield & gradient)
	{
		// Set to zero
		Vectormath::fill(gradient, {0,0,0});

		// External field
		Gradient_Zeeman(gradient);

		// Anisotropy
		Gradient_Anisotropy(spins, gradient);

		// Pairs
		//		Loop over periodicity
		for (int i_periodicity = 0; i_periodicity < 8; ++i_periodicity)
		{
			// Check if boundary conditions contain this periodicity
			if ((i_periodicity == 0)
				|| (i_periodicity == 1 && this->boundary_conditions[0])
				|| (i_periodicity == 2 && this->boundary_conditions[1])
				|| (i_periodicity == 3 && this->boundary_conditions[2])
				|| (i_periodicity == 4 && this->boundary_conditions[0] && this->boundary_conditions[1])
				|| (i_periodicity == 5 && this->boundary_conditions[0] && this->boundary_conditions[2])
				|| (i_periodicity == 6 && this->boundary_conditions[1] && this->boundary_conditions[2])
				|| (i_periodicity == 7 && this->boundary_conditions[0] && this->boundary_conditions[1] && this->boundary_conditions[2]))
			{
				//		Fields of this periodicity
				// Exchange
				this->Gradient_Exchange(spins, Exchange_indices[i_periodicity], Exchange_magnitude[i_periodicity], gradient);
				// DMI
				this->Gradient_DMI(spins, DMI_indices[i_periodicity], DMI_magnitude[i_periodicity], DMI_normal[i_periodicity], gradient);
				// DD
				this->Gradient_DD(spins, DD_indices[i_periodicity], DD_magnitude[i_periodicity], DD_normal[i_periodicity], gradient);
				// Quadruplet
				this->Gradient_Quadruplet(spins, Quadruplet_indices[i_periodicity], Quadruplet_magnitude[i_periodicity], gradient);
			}
		}

		// Triplet Interactions

		// Quadruplet Interactions

		cudaDeviceSynchronize();
	}


	__global__ void CU_Gradient_Zeeman( const int * external_field_index, const scalar * external_field_magnitude, const Vector3 * external_field_normal, Vector3 * gradient, size_t size)
	{
		for(auto idx = blockIdx.x * blockDim.x + threadIdx.x;
			idx < size;
			idx +=  blockDim.x * gridDim.x)
		{
			int ispin = external_field_index[idx];
			for (int dim=0; dim<3 ; dim++)
			{
				atomicAdd(&gradient[ispin][dim], -external_field_magnitude[idx]*external_field_normal[idx][dim]);
			}
		}
	}
	void Hamiltonian_Heisenberg::Gradient_Zeeman(vectorfield & gradient)
	{
		int size = this->external_field_index.size();
		CU_Gradient_Zeeman<<<(size+1023)/1024, 1024>>>( this->external_field_index.data(), this->external_field_magnitude.data(), this->external_field_normal.data(), gradient.data(), size );
	}


	__global__ void CU_Gradient_Anisotropy(const Vector3 * spins, const int * anisotropy_index, const scalar * anisotropy_magnitude, const Vector3 * anisotropy_normal, Vector3 * gradient, size_t size)
	{
		for(auto idx = blockIdx.x * blockDim.x + threadIdx.x;
			idx < size;
			idx +=  blockDim.x * gridDim.x)
		{
			int ispin = anisotropy_index[idx];
			scalar sc = -2 * anisotropy_magnitude[idx] * anisotropy_normal[idx].dot(spins[ispin]);
			for (int dim=0; dim<3 ; dim++)
			{
				atomicAdd(&gradient[ispin][dim], sc*anisotropy_normal[idx][dim]);
			}
		}
	}
	void Hamiltonian_Heisenberg::Gradient_Anisotropy(const vectorfield & spins, vectorfield & gradient)
	{
		int size = this->anisotropy_index.size();
		CU_Gradient_Anisotropy<<<(size+1023)/1024, 1024>>>( spins.data(), this->anisotropy_index.data(), this->anisotropy_magnitude.data(), this->anisotropy_normal.data(), gradient.data(), size );
	}


	__global__ void CU_Gradient_Exchange(const Vector3 * spins, const indexPair * pairs, const scalar * J_ij, Vector3 * gradient, size_t size)
	{
		for(auto iPair = blockIdx.x * blockDim.x + threadIdx.x;
			iPair < size;
			iPair +=  blockDim.x * gridDim.x)
		{
			int ispin = pairs[iPair][0];
			int jspin = pairs[iPair][1];
			for (int dim=0; dim<3 ; dim++)
			{
				atomicAdd(&gradient[ispin][dim], -J_ij[iPair]*spins[jspin][dim]);
				atomicAdd(&gradient[jspin][dim], -J_ij[iPair]*spins[ispin][dim]);
			}
		}
	}
	void Hamiltonian_Heisenberg::Gradient_Exchange(const vectorfield & spins, indexPairs & indices, scalarfield & J_ij, vectorfield & gradient)
	{
		int size = indices.size();
		CU_Gradient_Exchange<<<(size+1023)/1024, 1024>>>( spins.data(), indices.data(), J_ij.data(), gradient.data(), size );
	}


	__global__ void CU_Gradient_DMI(const Vector3 * spins, const indexPair * pairs, const scalar * DMI_magnitude, const Vector3 * DMI_normal, Vector3 * gradient, size_t size)
	{
		for(auto iPair = blockIdx.x * blockDim.x + threadIdx.x;
			iPair < size;
			iPair +=  blockDim.x * gridDim.x)
		{
			int ispin = pairs[iPair][0];
			int jspin = pairs[iPair][1];
			Vector3 jcross = DMI_magnitude[iPair]*spins[jspin].cross(DMI_normal[iPair]);
			Vector3 icross = DMI_magnitude[iPair]*spins[ispin].cross(DMI_normal[iPair]);
			for (int dim=0; dim<3 ; dim++)
			{
				atomicAdd(&gradient[ispin][dim], -jcross[dim]);
				atomicAdd(&gradient[jspin][dim],  icross[dim]);
			}
		}
	}
	void Hamiltonian_Heisenberg::Gradient_DMI(const vectorfield & spins, indexPairs & indices, scalarfield & DMI_magnitude, vectorfield & DMI_normal, vectorfield & gradient)
	{
		int size = indices.size();
		CU_Gradient_DMI<<<(size+1023)/1024, 1024>>>( spins.data(), indices.data(), DMI_magnitude.data(), DMI_normal.data(), gradient.data(), size );
	}


	void Hamiltonian_Heisenberg::Gradient_DD(const vectorfield & spins, indexPairs & indices, scalarfield & DD_magnitude, vectorfield & DD_normal, vectorfield & gradient)
	{
		//scalar mult = Utility::Constants::mu_B*Utility::Constants::mu_B*1.0 / 4.0 / M_PI; // multiply with mu_B^2
		scalar mult = 0.0536814951168; // mu_0*mu_B**2/(4pi*10**-30) -- the translations are in angstr�m, so the |r|[m] becomes |r|[m]*10^-10
		
		for (unsigned int i_pair = 0; i_pair < indices.size(); ++i_pair)
		{
			if (DD_magnitude[i_pair] > 0.0)
			{
				scalar skalar_contrib = mult / std::pow(DD_magnitude[i_pair], 3.0);
				gradient[indices[i_pair][0]] -= skalar_contrib * (3 * DD_normal[i_pair] * spins[indices[i_pair][1]].dot(DD_normal[i_pair]) - spins[indices[i_pair][1]]);
				gradient[indices[i_pair][1]] -= skalar_contrib * (3 * DD_normal[i_pair] * spins[indices[i_pair][0]].dot(DD_normal[i_pair]) - spins[indices[i_pair][0]]);
			}
		}
	}//end Field_DipoleDipole


	void Hamiltonian_Heisenberg::Gradient_Quadruplet(const vectorfield & spins, indexQuadruplets & indices, scalarfield & magnitude, vectorfield & gradient)
	{
		for (unsigned int i_pair = 0; i_pair < indices.size(); ++i_pair)
		{
			gradient[indices[i_pair][0]] -= magnitude[i_pair] * spins[indices[i_pair][1]] * (spins[indices[i_pair][2]].dot(spins[indices[i_pair][3]]));
			gradient[indices[i_pair][1]] -= magnitude[i_pair] * spins[indices[i_pair][0]] *  (spins[indices[i_pair][2]].dot(spins[indices[i_pair][3]]));
			gradient[indices[i_pair][2]] -= magnitude[i_pair] * (spins[indices[i_pair][0]].dot(spins[indices[i_pair][1]])) * spins[indices[i_pair][3]];
			gradient[indices[i_pair][3]] -= magnitude[i_pair] * (spins[indices[i_pair][0]].dot(spins[indices[i_pair][1]])) * spins[indices[i_pair][2]];
		}
	}


	void Hamiltonian_Heisenberg::Hessian(const vectorfield & spins, MatrixX & hessian)
	{
		int nos = spins.size();

		// Set to zero
		// for (auto& h : hessian) h = 0;
		hessian.setZero();

		// Single Spin elements
		for (int alpha = 0; alpha < 3; ++alpha)
		{
			for (unsigned int i = 0; i < anisotropy_index.size(); ++i)
			{
				int idx = anisotropy_index[i];
				// scalar x = -2.0*this->anisotropy_magnitude[i] * std::pow(this->anisotropy_normal[i][alpha], 2);
				hessian(3*idx + alpha, 3*idx + alpha) += -2.0*this->anisotropy_magnitude[i]*std::pow(this->anisotropy_normal[i][alpha],2);
			}
		}

		// std::cerr << "calculated hessian" << std::endl;

		 // Spin Pair elements
		 for (int i_periodicity = 0; i_periodicity < 8; ++i_periodicity)
		 {
		 	//		Check if boundary conditions contain this periodicity
		 	if ((i_periodicity == 0)
		 		|| (i_periodicity == 1 && this->boundary_conditions[0])
		 		|| (i_periodicity == 2 && this->boundary_conditions[1])
		 		|| (i_periodicity == 3 && this->boundary_conditions[2])
		 		|| (i_periodicity == 4 && this->boundary_conditions[0] && this->boundary_conditions[1])
		 		|| (i_periodicity == 5 && this->boundary_conditions[0] && this->boundary_conditions[2])
		 		|| (i_periodicity == 6 && this->boundary_conditions[1] && this->boundary_conditions[2])
		 		|| (i_periodicity == 7 && this->boundary_conditions[0] && this->boundary_conditions[1] && this->boundary_conditions[2]))
		 	{
		 		//		Loop over pairs of this periodicity
		 		// Exchange
		 		for (unsigned int i_pair = 0; i_pair < this->Exchange_indices[i_periodicity].size(); ++i_pair)
		 		{
		 			for (int alpha = 0; alpha < 3; ++alpha)
		 			{
		 				int idx_i = 3*Exchange_indices[i_periodicity][i_pair][0] + alpha;
		 				int idx_j = 3*Exchange_indices[i_periodicity][i_pair][1] + alpha;
		 				hessian(idx_i,idx_j) += -Exchange_magnitude[i_periodicity][i_pair];
		 				hessian(idx_j,idx_i) += -Exchange_magnitude[i_periodicity][i_pair];
		 			}
		 		}
		 		// DMI
		 		for (unsigned int i_pair = 0; i_pair < this->DMI_indices[i_periodicity].size(); ++i_pair)
		 		{
		 			for (int alpha = 0; alpha < 3; ++alpha)
		 			{
		 				for (int beta = 0; beta < 3; ++beta)
		 				{
		 					int idx_i = 3*DMI_indices[i_periodicity][i_pair][0] + alpha;
		 					int idx_j = 3*DMI_indices[i_periodicity][i_pair][1] + beta;
		 					if ( (alpha == 0 && beta == 1) )
		 					{
		 						hessian(idx_i,idx_j) +=
		 							-DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][2];
		 						hessian(idx_j,idx_i) +=
		 							-DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][2];
		 					}
		 					else if ( (alpha == 1 && beta == 0) )
		 					{
		 						hessian(idx_i,idx_j) +=
		 							DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][2];
		 						hessian(idx_j,idx_i) +=
		 							DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][2];
		 					}
		 					else if ( (alpha == 0 && beta == 2) )
		 					{
		 						hessian(idx_i,idx_j) +=
		 							DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][1];
		 						hessian(idx_j,idx_i) +=
		 							DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][1];
		 					}
		 					else if ( (alpha == 2 && beta == 0) )
		 					{
		 						hessian(idx_i,idx_j) +=
		 							-DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][1];
		 						hessian(idx_j,idx_i) +=
		 							-DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][1];
		 					}
		 					else if ( (alpha == 1 && beta == 2) )
		 					{
		 						hessian(idx_i,idx_j) +=
		 							-DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][0];
		 						hessian(idx_j,idx_i) +=
		 							-DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][0];
		 					}
		 					else if ( (alpha == 2 && beta == 1) )
		 					{
		 						hessian(idx_i,idx_j) +=
		 							DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][0];
		 						hessian(idx_j,idx_i) +=
		 							DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][0];
		 					}
		 				}
		 			}
		 		}
		 //		// Dipole-Dipole
		 //		for (unsigned int i_pair = 0; i_pair < this->DD_indices[i_periodicity].size(); ++i_pair)
		 //		{
		 //			// indices
		 //			int idx_1 = DD_indices[i_periodicity][i_pair][0];
		 //			int idx_2 = DD_indices[i_periodicity][i_pair][1];
		 //			// prefactor
		 //			scalar prefactor = 0.0536814951168
		 //				* this->mu_s[idx_1] * this->mu_s[idx_2]
		 //				/ std::pow(DD_magnitude[i_periodicity][i_pair], 3);
		 //			// components
		 //			for (int alpha = 0; alpha < 3; ++alpha)
		 //			{
		 //				for (int beta = 0; beta < 3; ++beta)
		 //				{
		 //					int idx_h = idx_1 + alpha*nos + 3 * nos*(idx_2 + beta*nos);
		 //					if (alpha == beta)
		 //						hessian[idx_h] += prefactor;
		 //					hessian[idx_h] += -3.0*prefactor*DD_normal[i_periodicity][i_pair][alpha] * DD_normal[i_periodicity][i_pair][beta];
		 //				}
		 //			}
		 //		}
		 	}// end if periodicity
		 }// end for periodicity
	}

	// Hamiltonian name as string
	static const std::string name = "Heisenberg Heisenberg";
	const std::string& Hamiltonian_Heisenberg::Name() { return name; }
}

#endif