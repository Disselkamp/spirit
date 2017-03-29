#ifndef USE_CUDA

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
	inline bool boundary_conditions_fulfilled(const intfield & n_cells, const boolfield & boundary_conditions, const std::array<int,3> & translations_i, const std::array<int,3> & translations_j)
	{
		int da = translations_i[0]+translations_j[0];
		int db = translations_i[1]+translations_j[1];
		int dc = translations_i[2]+translations_j[2];
		return  ( ( boundary_conditions[0] || (0 <= da && da < n_cells[0]) ) &&
				  ( boundary_conditions[1] || (0 <= db && db < n_cells[1]) ) &&
				  ( boundary_conditions[2] || (0 <= dc && dc < n_cells[2]) ) );
	}

	inline int idx_from_translations(const intfield & n_cells, const intfield & n_spins_basic_domain, const std::array<int,3> & translations_i, const std::array<int,3> & translations)
	{
		int Na = n_cells[0];
		int Nb = n_cells[1];
		int Nc = n_cells[2];
		int N  = n_spins_basic_domain[0];
		
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

	inline std::array<int,3> translations_from_idx(const intfield & n_cells, const intfield & n_spins_basic_domain, int idx)
	{
		std::array<int,3> ret;
		int Na = n_cells[0];
		int Nb = n_cells[1];
		int Nc = n_cells[2];
		int N  = n_spins_basic_domain[0];
		ret[2] = idx/(Na*Nb);
		ret[1] = (idx-ret[2]*Na*Nb)/Na;
		ret[0] = idx-ret[2]*Na*Nb-ret[1]*Na;
		return ret;
	}

	Hamiltonian_Heisenberg::Hamiltonian_Heisenberg(
			scalarfield mu_s,
			intfield external_field_index, scalarfield external_field_magnitude, vectorfield external_field_normal,
			intfield anisotropy_index, scalarfield anisotropy_magnitude, vectorfield anisotropy_normal,
			pairfield Exchange_pairs, scalarfield Exchange_magnitude,
			pairfield DMI_pairs, scalarfield DMI_magnitude, vectorfield DMI_normal,
			pairfield DD_pairs, scalarfield DD_magnitude, vectorfield DD_normal,
			quadrupletfield quadruplets, scalarfield quadruplet_magnitude,
			std::shared_ptr<Data::Geometry> geometry,
			std::vector<bool> boundary_conditions
	) :
		Hamiltonian(boundary_conditions),
		geometry(geometry),
		mu_s(mu_s),
		external_field_index(external_field_index), external_field_magnitude(external_field_magnitude), external_field_normal(external_field_normal),
		anisotropy_index(anisotropy_index), anisotropy_magnitude(anisotropy_magnitude), anisotropy_normal(anisotropy_normal),
		Exchange_pairs(Exchange_pairs), Exchange_magnitude(Exchange_magnitude),
		DMI_pairs(DMI_pairs), DMI_magnitude(DMI_magnitude), DMI_normal(DMI_normal),
		DD_pairs(DD_pairs), DD_magnitude(DD_magnitude), DD_normal(DD_normal),
		quadruplets(quadruplets), quadruplet_magnitude(quadruplet_magnitude)
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
		if (this->Exchange_pairs.size() > 0)
		{
			this->energy_contributions_per_spin.push_back({"Exchange", scalarfield(0) });
			this->idx_exchange = this->energy_contributions_per_spin.size()-1;
		}
		else this->idx_exchange = -1;
		// DMI
		if (this->DMI_pairs.size() > 0)
		{
			this->energy_contributions_per_spin.push_back({"DMI", scalarfield(0) });
			this->idx_dmi = this->energy_contributions_per_spin.size()-1;
		}
		else this->idx_dmi = -1;
		// Dipole-Dipole
		if (this->DD_pairs.size() > 0)
		{
			this->energy_contributions_per_spin.push_back({"DD", scalarfield(0) });
			this->idx_dd = this->energy_contributions_per_spin.size()-1;
		}
		else this->idx_dd = -1;
		// Quadruplet
		if (this->quadruplets.size() > 0)
		{
			this->energy_contributions_per_spin.push_back({"Quadruplet", scalarfield(0) });
			this->idx_quadruplet = this->energy_contributions_per_spin.size()-1;
		}
		else this->idx_quadruplet = -1;
	}

	void Hamiltonian_Heisenberg::Energy_Contributions_per_Spin(const vectorfield & spins, std::vector<std::pair<std::string, scalarfield>> & contributions)
	{
		int nos = spins.size();
		for (auto& pair : energy_contributions_per_spin)
		{
			// Allocate if not already allocated
			if (pair.second.size() != nos) pair.second = scalarfield(nos, 0);
			// Otherwise set to zero
			else for (auto& pair : energy_contributions_per_spin) Vectormath::fill(pair.second, 0);
		}
		

		// External field
		if (this->idx_zeeman >=0 ) E_Zeeman(spins, energy_contributions_per_spin[idx_zeeman].second);

		// Anisotropy
		if (this->idx_anisotropy >=0 ) E_Anisotropy(spins, energy_contributions_per_spin[idx_anisotropy].second);

		// Pairs
		// Exchange
		if (this->idx_exchange >=0 )   E_Exchange(spins,energy_contributions_per_spin[idx_exchange].second);
		// DMI
		if (this->idx_dmi >=0 )        E_DMI(spins, energy_contributions_per_spin[idx_dmi].second);
		// DD
		if (this->idx_dd >=0 )         E_DD(spins, energy_contributions_per_spin[idx_dd].second);
		// Quadruplet
		if (this->idx_quadruplet >=0 ) E_Quadruplet(spins, energy_contributions_per_spin[idx_quadruplet].second);

		// Return
		//return this->E;
	}

	void Hamiltonian_Heisenberg::E_Zeeman(const vectorfield & spins, scalarfield & Energy)
	{
		for (unsigned int i = 0; i < this->external_field_index.size(); ++i)
		{
			Energy[external_field_index[i]] -= this->external_field_magnitude[i] * this->external_field_normal[i].dot(spins[external_field_index[i]]);
		}
	}

	void Hamiltonian_Heisenberg::E_Anisotropy(const vectorfield & spins, scalarfield & Energy)
	{
		for (unsigned int i = 0; i < this->anisotropy_index.size(); ++i)
		{
			Energy[anisotropy_index[i]] -= this->anisotropy_magnitude[i] * std::pow(anisotropy_normal[i].dot(spins[anisotropy_index[i]]), 2.0);
		}
	}

	void Hamiltonian_Heisenberg::E_Exchange(const vectorfield & spins, scalarfield & Energy)
	{
		for (unsigned int ispin = 0; ispin < spins.size(); ++ispin)
		{
			auto translations = translations_from_idx(geometry->n_cells, geometry->n_spins_basic_domain, ispin);
			for (unsigned int i_pair = 0; i_pair < Exchange_pairs.size(); ++i_pair)
			{
				if ( boundary_conditions_fulfilled(geometry->n_cells, boundary_conditions, translations, Exchange_pairs[i_pair].translations) )
				{
					int jspin = idx_from_translations(geometry->n_cells, geometry->n_spins_basic_domain, translations, Exchange_pairs[i_pair].translations);
					Energy[ispin] -= 0.5 * Exchange_magnitude[i_pair] * spins[ispin].dot(spins[jspin]);
				}
			}
		}
	}

	void Hamiltonian_Heisenberg::E_DMI(const vectorfield & spins, scalarfield & Energy)
	{
		for (unsigned int ispin = 0; ispin < spins.size(); ++ispin)
		{
			auto translations = translations_from_idx(geometry->n_cells, geometry->n_spins_basic_domain, ispin);
			for (unsigned int i_pair = 0; i_pair < DMI_pairs.size(); ++i_pair)
			{
				if ( boundary_conditions_fulfilled(geometry->n_cells, boundary_conditions, translations, DMI_pairs[i_pair].translations) )
				{
					int jspin = idx_from_translations(geometry->n_cells, geometry->n_spins_basic_domain, translations, DMI_pairs[i_pair].translations);
					Energy[ispin] -= 0.5 * DMI_magnitude[i_pair] * DMI_normal[i_pair].dot(spins[ispin].cross(spins[jspin]));
				}
			}
		}
	}

	void Hamiltonian_Heisenberg::E_DD(const vectorfield & spins, scalarfield & Energy)
	{
		//scalar mult = -Constants::mu_B*Constants::mu_B*1.0 / 4.0 / M_PI; // multiply with mu_B^2
		scalar mult = 0.5*0.0536814951168; // mu_0*mu_B**2/(4pi*10**-30) -- the translations are in angstr�m, so the |r|[m] becomes |r|[m]*10^-10
		scalar result = 0.0;

		for (unsigned int i_pair = 0; i_pair < DD_pairs.size(); ++i_pair)
		{
			if (DD_magnitude[i_pair] > 0.0)
			{
				// Energy[pairs[i_pair][0]] -= mult / std::pow(DD_magnitude[i_pair], 3.0) *
				// 	(3 * spins[pairs[i_pair][1]].dot(DD_normal[i_pair]) * spins[pairs[i_pair][0]].dot(DD_normal[i_pair]) - spins[pairs[i_pair][0]].dot(spins[pairs[i_pair][1]]));
				// Energy[pairs[i_pair][1]] -= mult / std::pow(DD_magnitude[i_pair], 3.0) *
				// 	(3 * spins[pairs[i_pair][1]].dot(DD_normal[i_pair]) * spins[pairs[i_pair][0]].dot(DD_normal[i_pair]) - spins[pairs[i_pair][0]].dot(spins[pairs[i_pair][1]]));
			}

		}
	}// end DipoleDipole


	void Hamiltonian_Heisenberg::E_Quadruplet(const vectorfield & spins, scalarfield & Energy)
	{
		for (unsigned int i_pair = 0; i_pair < quadruplets.size(); ++i_pair)
		{
			// Energy[pairs[i_pair][0]] -= 0.25*magnitude[i_pair] * (spins[pairs[i_pair][0]].dot(spins[pairs[i_pair][1]])) * (spins[pairs[i_pair][2]].dot(spins[pairs[i_pair][3]]));
			// Energy[pairs[i_pair][1]] -= 0.25*magnitude[i_pair] * (spins[pairs[i_pair][0]].dot(spins[pairs[i_pair][1]])) * (spins[pairs[i_pair][2]].dot(spins[pairs[i_pair][3]]));
			// Energy[pairs[i_pair][2]] -= 0.25*magnitude[i_pair] * (spins[pairs[i_pair][0]].dot(spins[pairs[i_pair][1]])) * (spins[pairs[i_pair][2]].dot(spins[pairs[i_pair][3]]));
			// Energy[pairs[i_pair][3]] -= 0.25*magnitude[i_pair] * (spins[pairs[i_pair][0]].dot(spins[pairs[i_pair][1]])) * (spins[pairs[i_pair][2]].dot(spins[pairs[i_pair][3]]));
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
		// Exchange
		this->Gradient_Exchange(spins, gradient);
		// DMI
		this->Gradient_DMI(spins, gradient);
		// DD
		this->Gradient_DD(spins, gradient);

		// Triplet

		// Quadruplet
		this->Gradient_Quadruplet(spins, gradient);
	}

	void Hamiltonian_Heisenberg::Gradient_Zeeman(vectorfield & gradient)
	{
		for (unsigned int i = 0; i < this->external_field_index.size(); ++i)
		{
			gradient[external_field_index[i]] -= this->external_field_magnitude[i] * this->external_field_normal[i];
		}
	}

	void Hamiltonian_Heisenberg::Gradient_Anisotropy(const vectorfield & spins, vectorfield & gradient)
	{
		for (unsigned int i = 0; i < this->anisotropy_index.size(); ++i)
		{
			gradient[anisotropy_index[i]] -= 2.0 * this->anisotropy_magnitude[i] * this->anisotropy_normal[i] * anisotropy_normal[i].dot(spins[anisotropy_index[i]]);
		}
	}

	void Hamiltonian_Heisenberg::Gradient_Exchange(const vectorfield & spins, vectorfield & gradient)
	{
		for (unsigned int ispin = 0; ispin < spins.size(); ++ispin)
		{
			auto translations = translations_from_idx(geometry->n_cells, geometry->n_spins_basic_domain, ispin);
			for (unsigned int i_pair = 0; i_pair < Exchange_pairs.size(); ++i_pair)
			{
				if ( boundary_conditions_fulfilled(geometry->n_cells, boundary_conditions, translations, Exchange_pairs[i_pair].translations) )
				{
					int jspin = idx_from_translations(geometry->n_cells, geometry->n_spins_basic_domain, translations, Exchange_pairs[i_pair].translations);
					gradient[ispin] -= Exchange_magnitude[i_pair] * spins[jspin];
				}
			}
		}
	}

	void Hamiltonian_Heisenberg::Gradient_DMI(const vectorfield & spins, vectorfield & gradient)
	{
		for (unsigned int ispin = 0; ispin < spins.size(); ++ispin)
		{
			auto translations = translations_from_idx(geometry->n_cells, geometry->n_spins_basic_domain, ispin);
			for (unsigned int i_pair = 0; i_pair < DMI_pairs.size(); ++i_pair)
			{
				if ( boundary_conditions_fulfilled(geometry->n_cells, boundary_conditions, translations, DMI_pairs[i_pair].translations) )
				{
					int jspin = idx_from_translations(geometry->n_cells, geometry->n_spins_basic_domain, translations, DMI_pairs[i_pair].translations);
					gradient[ispin] -= DMI_magnitude[i_pair] * spins[jspin].cross(DMI_normal[i_pair]);
				}
			}
		}
	}

	void Hamiltonian_Heisenberg::Gradient_DD(const vectorfield & spins, vectorfield & gradient)
	{
		//scalar mult = Constants::mu_B*Constants::mu_B*1.0 / 4.0 / M_PI; // multiply with mu_B^2
		scalar mult = 0.0536814951168; // mu_0*mu_B**2/(4pi*10**-30) -- the translations are in angstr�m, so the |r|[m] becomes |r|[m]*10^-10
		
		for (unsigned int i_pair = 0; i_pair < DD_pairs.size(); ++i_pair)
		{
			if (DD_magnitude[i_pair] > 0.0)
			{
				scalar skalar_contrib = mult / std::pow(DD_magnitude[i_pair], 3.0);
				// gradient[indices[i_pair][0]] -= skalar_contrib * (3 * DD_normal[i_pair] * spins[indices[i_pair][1]].dot(DD_normal[i_pair]) - spins[indices[i_pair][1]]);
				// gradient[indices[i_pair][1]] -= skalar_contrib * (3 * DD_normal[i_pair] * spins[indices[i_pair][0]].dot(DD_normal[i_pair]) - spins[indices[i_pair][0]]);
			}
		}
	}//end Field_DipoleDipole


	void Hamiltonian_Heisenberg::Gradient_Quadruplet(const vectorfield & spins, vectorfield & gradient)
	{
		for (unsigned int i_pair = 0; i_pair < quadruplets.size(); ++i_pair)
		{
			// gradient[indices[i_pair][0]] -= magnitude[i_pair] * spins[indices[i_pair][1]] * (spins[indices[i_pair][2]].dot(spins[indices[i_pair][3]]));
			// gradient[indices[i_pair][1]] -= magnitude[i_pair] * spins[indices[i_pair][0]] *  (spins[indices[i_pair][2]].dot(spins[indices[i_pair][3]]));
			// gradient[indices[i_pair][2]] -= magnitude[i_pair] * (spins[indices[i_pair][0]].dot(spins[indices[i_pair][1]])) * spins[indices[i_pair][3]];
			// gradient[indices[i_pair][3]] -= magnitude[i_pair] * (spins[indices[i_pair][0]].dot(spins[indices[i_pair][1]])) * spins[indices[i_pair][2]];
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

		//  // Spin Pair elements
		//  for (int i_periodicity = 0; i_periodicity < 8; ++i_periodicity)
		//  {
		//  	//		Check if boundary conditions contain this periodicity
		//  	if ((i_periodicity == 0)
		//  		|| (i_periodicity == 1 && this->boundary_conditions[0])
		//  		|| (i_periodicity == 2 && this->boundary_conditions[1])
		//  		|| (i_periodicity == 3 && this->boundary_conditions[2])
		//  		|| (i_periodicity == 4 && this->boundary_conditions[0] && this->boundary_conditions[1])
		//  		|| (i_periodicity == 5 && this->boundary_conditions[0] && this->boundary_conditions[2])
		//  		|| (i_periodicity == 6 && this->boundary_conditions[1] && this->boundary_conditions[2])
		//  		|| (i_periodicity == 7 && this->boundary_conditions[0] && this->boundary_conditions[1] && this->boundary_conditions[2]))
		//  	{
		//  		//		Loop over pairs of this periodicity
		//  		// Exchange
		//  		for (unsigned int i_pair = 0; i_pair < this->Exchange_pairs.size(); ++i_pair)
		//  		{
		//  			for (int alpha = 0; alpha < 3; ++alpha)
		//  			{
		//  				int idx_i = 3*Exchange_pairs[i_pair][0] + alpha;
		//  				int idx_j = 3*Exchange_pairs[i_pair][1] + alpha;
		//  				hessian(idx_i,idx_j) += -Exchange_magnitude[i_pair];
		//  				hessian(idx_j,idx_i) += -Exchange_magnitude[i_pair];
		//  			}
		//  		}
		//  		// DMI
		//  		for (unsigned int i_pair = 0; i_pair < this->DMI_pairs[i_periodicity].size(); ++i_pair)
		//  		{
		//  			for (int alpha = 0; alpha < 3; ++alpha)
		//  			{
		//  				for (int beta = 0; beta < 3; ++beta)
		//  				{
		//  					int idx_i = 3*DMI_pairs[i_periodicity][i_pair][0] + alpha;
		//  					int idx_j = 3*DMI_pairs[i_periodicity][i_pair][1] + beta;
		//  					if ( (alpha == 0 && beta == 1) )
		//  					{
		//  						hessian(idx_i,idx_j) +=
		//  							-DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][2];
		//  						hessian(idx_j,idx_i) +=
		//  							-DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][2];
		//  					}
		//  					else if ( (alpha == 1 && beta == 0) )
		//  					{
		//  						hessian(idx_i,idx_j) +=
		//  							DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][2];
		//  						hessian(idx_j,idx_i) +=
		//  							DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][2];
		//  					}
		//  					else if ( (alpha == 0 && beta == 2) )
		//  					{
		//  						hessian(idx_i,idx_j) +=
		//  							DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][1];
		//  						hessian(idx_j,idx_i) +=
		//  							DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][1];
		//  					}
		//  					else if ( (alpha == 2 && beta == 0) )
		//  					{
		//  						hessian(idx_i,idx_j) +=
		//  							-DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][1];
		//  						hessian(idx_j,idx_i) +=
		//  							-DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][1];
		//  					}
		//  					else if ( (alpha == 1 && beta == 2) )
		//  					{
		//  						hessian(idx_i,idx_j) +=
		//  							-DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][0];
		//  						hessian(idx_j,idx_i) +=
		//  							-DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][0];
		//  					}
		//  					else if ( (alpha == 2 && beta == 1) )
		//  					{
		//  						hessian(idx_i,idx_j) +=
		//  							DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][0];
		//  						hessian(idx_j,idx_i) +=
		//  							DMI_magnitude[i_periodicity][i_pair] * DMI_normal[i_periodicity][i_pair][0];
		//  					}
		//  				}
		//  			}
		//  		}
		//  //		// Dipole-Dipole
		//  //		for (unsigned int i_pair = 0; i_pair < this->DD_pairs[i_periodicity].size(); ++i_pair)
		//  //		{
		//  //			// indices
		//  //			int idx_1 = DD_pairs[i_periodicity][i_pair][0];
		//  //			int idx_2 = DD_pairs[i_periodicity][i_pair][1];
		//  //			// prefactor
		//  //			scalar prefactor = 0.0536814951168
		//  //				* this->mu_s[idx_1] * this->mu_s[idx_2]
		//  //				/ std::pow(DD_magnitude[i_periodicity][i_pair], 3);
		//  //			// components
		//  //			for (int alpha = 0; alpha < 3; ++alpha)
		//  //			{
		//  //				for (int beta = 0; beta < 3; ++beta)
		//  //				{
		//  //					int idx_h = idx_1 + alpha*nos + 3 * nos*(idx_2 + beta*nos);
		//  //					if (alpha == beta)
		//  //						hessian[idx_h] += prefactor;
		//  //					hessian[idx_h] += -3.0*prefactor*DD_normal[i_periodicity][i_pair][alpha] * DD_normal[i_periodicity][i_pair][beta];
		//  //				}
		//  //			}
		//  //		}
		//  	}// end if periodicity
		//  }// end for periodicity
	}

	// Hamiltonian name as string
	static const std::string name = "Heisenberg";
	const std::string& Hamiltonian_Heisenberg::Name() { return name; }
}

#endif