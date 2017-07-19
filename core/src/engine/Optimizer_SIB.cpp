#include <engine/Optimizer_SIB.hpp>
#include <engine/Vectormath.hpp>
#include <utility/Constants.hpp>

using namespace Utility;

namespace Engine
{
	Optimizer_SIB::Optimizer_SIB(std::shared_ptr<Engine::Method> method) :
        Optimizer(method)
    {
		this->xi = vectorfield(this->nos, {0,0,0});
		this->virtualforce = std::vector<vectorfield>(this->noi, vectorfield(this->nos));	// [noi][nos]
		
		this->spins_temp = std::vector<std::shared_ptr<vectorfield>>(this->noi);
		for (int i=0; i<this->noi; ++i) spins_temp[i] = std::shared_ptr<vectorfield>(new vectorfield(this->nos)); // [noi][nos]
    }

    void Optimizer_SIB::Iteration()
    {
		std::shared_ptr<Data::Spin_System> s;

		// Random Numbers
		for (int i = 0; i < this->noi; ++i)
		{
			s = method->systems[i];
			if (s->llg_parameters->temperature > 0)
			{
				this->epsilon = std::sqrt(2.0*s->llg_parameters->damping / (1.0 + std::pow(s->llg_parameters->damping, 2))*s->llg_parameters->temperature*Constants::k_B);
				// Precalculate RNs --> move this up into Iterate and add array dimension n for no of iterations?
				Vectormath::get_random_vectorfield(*s, epsilon, xi);
			}
		}

		// First part of the step
		this->method->Calculate_Force(configurations, force);
		for (int i = 0; i < this->noi; ++i)
		{
			s = method->systems[i];
			this->VirtualForce(*s->spins, *s->llg_parameters, force[i], xi, virtualforce[i]);
			Vectormath::transform(*s->spins, virtualforce[i], *spins_temp[i]);
			Vectormath::add_c_a(1, *s->spins, *spins_temp[i]);
			Vectormath::scale(*spins_temp[i], 0.5);
		}

		// Second part of the step
		this->method->Calculate_Force(this->spins_temp, force); ////// I cannot see a difference if this step is included or not...
		for (int i = 0; i < this->noi; ++i)
		{
			s = method->systems[i];
			this->VirtualForce(*spins_temp[i], *s->llg_parameters, force[i], xi, virtualforce[i]);
			Vectormath::transform(*s->spins, virtualforce[i], *s->spins);
		}
	}


	void Optimizer_SIB::VirtualForce(vectorfield & spins, Data::Parameters_Method_LLG & llg_params, vectorfield & gradient,  vectorfield & xi, vectorfield & force)
	{

		//========================= Init local vars ================================
		// time steps
		scalar damping = llg_params.damping;
		scalar sqrtdt = std::sqrt(llg_params.dt)/(1+std::pow(damping, 2)); // is this correct for the thermal term? 
		scalar dtg = llg_params.dt/(1+std::pow(damping, 2));
		scalar sqrtdtg = sqrtdt;
		// STT
		scalar a_j = llg_params.stt_magnitude;
		Vector3 s_c_vec = llg_params.stt_polarisation_normal; 
		//STT gradient
        // hardcoded parameters
		bool gradON = true; // switch between STT methods
		scalar b_j = 1.0;   // pre-factor b_j = u*mu_s/gamma (see bachelorthesis Constantin)
		scalar beta = llg_params.beta;  // non-adiabatic parameter of correction term
		Vector3 je = { 1,0,0 }; // direction of current
		vectorfield s_c_grad;
		//------------------------ End Init ----------------------------------------

		Vectormath::fill       (force, {0,0,0});
		Vectormath::add_c_a    (-0.5 * dtg, gradient, force);
		Vectormath::add_c_cross(-0.5 * dtg * damping, spins, gradient, force);

		// STT
		if (a_j > 0 && gradON == false)
		{
			Vectormath::add_c_a    ( 0.5 * dtg * a_j * damping, s_c_vec, force); // Vorzeichen!
			Vectormath::add_c_cross( 0.5 * dtg * a_j, s_c_vec, spins, force);
		}

		// STT gradient
		auto & geometry = *this->method->systems[0]->geometry;
		auto & hamiltonian = *this->method->systems[0]->hamiltonian;
		if (a_j > 0 && gradON == true)
		{
			Vectormath::gradient   (spins, hamiltonian, geometry, je, s_c_grad); // s_c_grad = (j_e*grad)*S
			Vectormath::add_c_a    ( -0.5 * dtg * a_j * ( damping - beta ), s_c_grad, force); // TODO: a_j durch b_j ersetzen 
			Vectormath::add_c_cross( -0.5 * dtg * a_j * ( 1 + beta * damping ), s_c_grad, spins, force); // TODO: a_j durch b_j ersetzen 
			// Gradient in current richtung, daher => *(-1)
		}

		// Temperature
		if (llg_params.temperature > 0)
		{
			Vectormath::add_c_a    (-0.5 * sqrtdtg, xi, force);
			Vectormath::add_c_cross(-0.5 * sqrtdtg * damping, spins, xi, force);
		}
	}


    // Optimizer name as string
    std::string Optimizer_SIB::Name() { return "SIB"; }
    std::string Optimizer_SIB::FullName() { return "Semi-implicit B"; }
}
