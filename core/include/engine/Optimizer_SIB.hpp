#pragma once
#ifndef OPTIMIZER_SIB_H
#define OPTIMIZER_SIB_H

#include <vector>

#include <Eigen/Dense>

#include "Core_Defines.h"
#include <engine/Optimizer.hpp>
#include <data/Spin_System_Chain.hpp>

namespace Engine
{
	/*
		Semi-Implicit Method B (SIB) Optimizer:
			The Spin System is optimized w.r.t. the force while following the physical LLG trajectory.
			Note that this means this is not a direct optimizer and the system posesses "momentum".
			Method taken from: Mentink et. al., Stable and fast semi-implicit integration of the stochastic Landau�Lifshitz equation, J. Phys.: Condens. Matter 22 (2010) 176001 (12pp)
	*/
	class Optimizer_SIB : public Optimizer
	{

	public:
		Optimizer_SIB(std::shared_ptr<Engine::Method> method);
		
		// One step in the optimization
		void Iteration() override;
		
		// Optimizer name as string
		std::string Name() override;
		std::string FullName() override;

    private:
		// Temporary Spins arrays
		std::vector<std::shared_ptr<std::vector<Vector3>>> spins_temp;

		// Virtual Heun Forces used in the Steps
		std::vector<std::vector<Vector3>> virtualforce;

		// Random vector array
		std::vector<Vector3> xi;
		// Some variable
		scalar epsilon;

		// Generate an array of random spins?
		void Gen_Xi(Data::Spin_System & s, std::vector<Vector3> & xi, scalar eps);
		// Calculate the virtual Heun force to be used in the Steps
		void VirtualForce(const int nos, std::vector<Vector3> & spins, Data::Parameters_Method_LLG & llg_params, std::vector<Vector3> & eff_field, std::vector<Vector3> & xi, std::vector<Vector3> & force);
		// First Part of one Optimization step
		void FirstStep(const int nos, std::vector<Vector3> & spins, std::vector<Vector3> & force, std::vector<Vector3> & spins_temp);
		// Second Part of one Optimization step
		void SecondStep(const int nos, std::vector<Vector3> & force, std::vector<Vector3> & spins);

    };
}

#endif