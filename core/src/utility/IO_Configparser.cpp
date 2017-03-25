﻿#include <utility/IO.hpp>
#include <utility/IO_Filter_File_Handle.hpp>
#include <engine/Vectormath.hpp>
#include <engine/Neighbours.hpp>
#include <utility/Constants.hpp>
#include <utility/Logging.hpp>
#include <utility/Exception.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>

namespace Utility
{
	namespace IO
	{
		void Log_from_Config(const std::string configFile)
		{
			// Verbosity and Reject Level are read as integers
			int i_print_level = 5, i_accept_level = 5;
			std::string output_folder = ".";
			bool save_output = true, save_input = false;

			//------------------------------- Parser --------------------------------
			if (configFile != "")
			{
				try
				{
					Log(Log_Level::Info, Log_Sender::IO, "Building Log");
					IO::Filter_File_Handle myfile(configFile);

					// Accept Level
					if (myfile.Find("log_accept")) myfile.iss >> i_accept_level;
					else Log(Log_Level::Error, Log_Sender::IO, "Keyword 'log_accept' not found. Using default: " + std::to_string(i_accept_level));

					// Print level
					if (myfile.Find("log_print")) myfile.iss >> i_print_level;
					else Log(Log_Level::Error, Log_Sender::IO, "Keyword 'log_print' not found. Using default: " + std::to_string(i_print_level));

					// Output folder
					if (myfile.Find("log_output_folder")) myfile.iss >> output_folder;
					else Log(Log_Level::Error, Log_Sender::IO, "Keyword 'log_output_folder' not found. Using default: '" + output_folder + "'");
					
					// Save Output (Log Messages)
					if (myfile.Find("log_output_save")) myfile.iss >> save_output;
					else Log(Log_Level::Error, Log_Sender::IO, "Keyword 'log_output_save' not found. Using default: " + std::to_string(save_output));
					
					// Save Input (parameters from config file and defaults)
					if (myfile.Find("log_input_save")) myfile.iss >> save_input;
					else Log(Log_Level::Error, Log_Sender::IO, "Keyword 'log_input_save' not found. Using default: " + std::to_string(save_input));

				}// end try
				catch (Exception ex) {
					if (ex == Exception::File_not_Found) {
						Log(Log_Level::Error, Log_Sender::IO, "Log_Levels: Unable to open Config File " + configFile + " Leaving values at default.");
					}
					else throw ex;
				}// end catch
			}
			// Log the parameters
			Log(Log_Level::Parameter, Log_Sender::IO, "Log accept level  = " + std::to_string(i_accept_level));
			Log(Log_Level::Parameter, Log_Sender::IO, "Log print level   = " + std::to_string(i_print_level));
			Log(Log_Level::Parameter, Log_Sender::IO, "Log output folder = " + output_folder);
			Log(Log_Level::Parameter, Log_Sender::IO, "Log output save   = " + std::to_string(save_output));
			Log(Log_Level::Parameter, Log_Sender::IO, "Log input save    = " + std::to_string(save_input));
			// Update the Log
			Log.accept_level  = Log_Level(i_accept_level);
			Log.print_level   = Log_Level(i_print_level);
			Log.output_folder = output_folder;
			Log.save_output   = save_output;
			Log.save_input    = save_input;
		}// End Log_Levels_from_Config


		std::unique_ptr<Data::Spin_System> Spin_System_from_Config(const std::string configFile)
		{
			Log(Log_Level::Info, Log_Sender::IO, "-------------- Initialising Spin System ------------");
			// ----------------------------------------------------------------------------------------------
			// Geometry
			auto geometry = Geometry_from_Config(configFile);
			// LLG Parameters
			auto llg_params = Parameters_Method_LLG_from_Config(configFile);
			// Hamiltonian
			auto hamiltonian = std::move(Hamiltonian_from_Config(configFile, geometry));
			// Spin System
			auto system = std::unique_ptr<Data::Spin_System>(new Data::Spin_System(std::move(hamiltonian), geometry, std::move(llg_params), false));
			// ----------------------------------------------------------------------------------------------
			Log(Log_Level::Info, Log_Sender::IO, "-------------- Spin System Initialised -------------");

			// Return
			return system;
		}// End Spin_System_from_Config		


		void Basis_from_Config(const std::string configFile, std::vector<Vector3> & basis, std::vector<Vector3> & basis_atoms, scalar & lattice_constant)
		{
			// ---------- Default values
			// Lattice constant [Angtrom]
			lattice_constant = 1.0;
			// Basis: vector {a, b, c}
			basis = { Vector3{1,0,0}, Vector3{0,1,0}, Vector3{0,0,1} };
			// Atoms in the basis [dim][n_basis_atoms]
			basis_atoms = { Vector3{0,0,0} };
			// NoS in the basic domain (= unit cell for periodic lattices)
			int n_spins_basic_domain = 0;
			
			Log(Log_Level::Info, Log_Sender::IO, "Basis: building");

			if (configFile != "")
			{
				try {
					IO::Filter_File_Handle myfile(configFile);

					myfile.Read_Single(lattice_constant, "lattice_constant");

					// Utility 1D array to build vectors and use Vectormath
					Vector3 build_array = { 0, 0, 0 };

					if (myfile.Find("basis"))
					{
						// Read the basis vectors a, b, c
						myfile.GetLine();
						myfile.iss >> basis[0][0] >> basis[0][1] >> basis[0][2];
						myfile.GetLine();
						myfile.iss >> basis[1][0] >> basis[1][1] >> basis[1][2];
						myfile.GetLine();
						myfile.iss >> basis[2][0] >> basis[2][1] >> basis[2][2];

						// Read no_spins_basic_domain and atoms in basis
						myfile.GetLine();
						myfile.iss >> n_spins_basic_domain;
						basis_atoms = std::vector<Vector3>(n_spins_basic_domain);

						// Read spins per basic domain
						for (int iatom = 0; iatom < n_spins_basic_domain; ++iatom)
						{
							myfile.GetLine();
							myfile.iss >> basis_atoms[iatom][0] >> basis_atoms[iatom][1] >> basis_atoms[iatom][2];
							// Get x,y,z of component of spin_pos in unit of length (instead of in units of a,b,c)
							build_array = basis[0] * basis_atoms[iatom][0] + basis[1] * basis_atoms[iatom][1] + basis[2] * basis_atoms[iatom][2];
							basis_atoms[iatom] = lattice_constant * build_array;
						}// endfor iatom

					}// end find "basis"
					else {
						Log(Log_Level::Error, Log_Sender::IO, "Keyword 'basis' not found. Using Default (sc)");
					}
				}// end try
				catch (Exception ex) {
					if (ex == Exception::File_not_Found)
					{
						Log(Log_Level::Error, Log_Sender::IO, "Basis: Unable to open Config File " + configFile + " Leaving values at default.");
						throw Exception::System_not_Initialized;
					}
					else throw ex;
				}// end catch
			}
			else Log(Log_Level::Warning, Log_Sender::IO, "Basis: No config file specified. Leaving values at default.");
			
			// Log the parameters
			Log(Log_Level::Parameter, Log_Sender::IO, "Lattice constant = " + std::to_string(lattice_constant) + " angstrom");
			Log(Log_Level::Debug, Log_Sender::IO, "Basis: vectors in units of lattice constant");
			Log(Log_Level::Debug, Log_Sender::IO, "        a = " + std::to_string(basis[0][0]/lattice_constant) + " " + std::to_string(basis[0][1]/lattice_constant) + " " + std::to_string(basis[0][2]/lattice_constant));
			Log(Log_Level::Debug, Log_Sender::IO, "        b = " + std::to_string(basis[1][0]/lattice_constant) + " " + std::to_string(basis[1][1]/lattice_constant) + " " + std::to_string(basis[1][2]/lattice_constant));
			Log(Log_Level::Debug, Log_Sender::IO, "        c = " + std::to_string(basis[2][0]/lattice_constant) + " " + std::to_string(basis[2][1]/lattice_constant) + " " + std::to_string(basis[2][2]/lattice_constant));
			Log(Log_Level::Parameter, Log_Sender::IO, "Basis: vectors");
			Log(Log_Level::Parameter, Log_Sender::IO, "        a = " + std::to_string(basis[0][0]) + " " + std::to_string(basis[0][1]) + " " + std::to_string(basis[0][2]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        b = " + std::to_string(basis[1][0]) + " " + std::to_string(basis[1][1]) + " " + std::to_string(basis[1][2]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        c = " + std::to_string(basis[2][0]) + " " + std::to_string(basis[2][1]) + " " + std::to_string(basis[2][2]));
			Log(Log_Level::Parameter, Log_Sender::IO, "Basis: " + std::to_string(n_spins_basic_domain) + " atom(s) at the following positions:");
			for (int iatom = 0; iatom < n_spins_basic_domain; ++iatom)
			{
				Log(Log_Level::Parameter, Log_Sender::IO, "            " + std::to_string(iatom) + " = " + std::to_string(basis_atoms[iatom][0]) + " " + std::to_string(basis_atoms[iatom][1]) + " " + std::to_string(basis_atoms[iatom][2]));
			}
			Log(Log_Level::Info, Log_Sender::IO, "Basis: built");
		}// End Basis_from_Config

		std::shared_ptr<Data::Geometry> Geometry_from_Config(const std::string configFile)
		{
			//-------------- Insert default values here -----------------------------
			// Basis from separate file?
			std::string basis_file = "";
			// Basis: vector {a, b, c}
			std::vector<Vector3> basis = { Vector3{1,0,0}, Vector3{0,1,0}, Vector3{0,0,1} };
			// Atoms in the basis [dim][n_basis_atoms]
			std::vector<Vector3> basis_atoms = { Vector3{0,0,0} };
			// Lattice Constant [Angstrom]
			scalar lattice_constant = 1;
			// Translation vectors [dim][nov]
			std::vector<Vector3> translation_vectors = { Vector3{1,0,0}, Vector3{0,1,0}, Vector3{0,0,1} };
			// Number of translations nT for each basis direction
			std::vector<int> n_cells = { 100, 100, 1 };
			// Number of Spins
			int nos;
			vectorfield spin_pos;

			// Utility 1D array to build vectors and use Vectormath
			Vector3 build_array = { 0, 0, 0 };

			Log(Log_Level::Info, Log_Sender::IO, "Geometry: building");
			//------------------------------- Parser --------------------------------
			// iteration variables
			int iatom = 0, dim = 0;
			if (configFile != "")
			{
				try {
					Log(Log_Level::Info, Log_Sender::IO, "Reading Geometry Parameters");
					IO::Filter_File_Handle myfile(configFile);

					// Read Shape of spins in term of the basis
					if (myfile.Find("translation_vectors"))
					{
						// Read translation vectors into translation_vectors & nTa, nTb, nTc
						myfile.GetLine();
						myfile.iss >> translation_vectors[0][0] >> translation_vectors[0][1] >> translation_vectors[0][2] >> n_cells[0];
						myfile.GetLine();
						myfile.iss >> translation_vectors[1][0] >> translation_vectors[1][1] >> translation_vectors[1][2] >> n_cells[1];
						myfile.GetLine();
						myfile.iss >> translation_vectors[2][0] >> translation_vectors[2][1] >> translation_vectors[2][2] >> n_cells[2];
					}// finish Reading Shape in terms of basis
					else {
						Log(Log_Level::Error, Log_Sender::IO, "Keyword 'translation_vectors' not found. Using default. (sc 30x30x0)");
					}
					// Read Basis
						
					if (myfile.Find("basis_from_config"))
					{
						myfile.iss >> basis_file;
						Basis_from_Config(basis_file, basis, basis_atoms, lattice_constant);
					}
					else if (myfile.Find("basis"))
					{
						Basis_from_Config(configFile, basis, basis_atoms, lattice_constant);
					}
					else {
						Log(Log_Level::Error, Log_Sender::IO, "Neither Keyword 'basis_from_config', nor Keyword 'basis' found. Using Default (sc)");
					}// end Basis
				}// end try
				catch (Exception ex)
				{
					if (ex == Exception::File_not_Found)
					{
						Log(Log_Level::Error, Log_Sender::IO, "Geometry: Unable to open Config File " + configFile + " Leaving values at default.");
					}
					else throw ex;
				}// end catch
			}// end if file=""
			else Log(Log_Level::Warning, Log_Sender::IO, "Geometry: Using default configuration!");

			Log(Log_Level::Parameter, Log_Sender::IO, "Translation: vectors");
			Log(Log_Level::Parameter, Log_Sender::IO, "        a = " + std::to_string(translation_vectors[0][0]) + " " + std::to_string(translation_vectors[1][0]) + " " + std::to_string(translation_vectors[2][0]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        b = " + std::to_string(translation_vectors[0][1]) + " " + std::to_string(translation_vectors[1][1]) + " " + std::to_string(translation_vectors[2][1]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        c = " + std::to_string(translation_vectors[0][2]) + " " + std::to_string(translation_vectors[1][2]) + " " + std::to_string(translation_vectors[2][2]));

			// Get x,y,z of component of translation_vectors in unit of length (instead of in units of a,b,c)
			for (dim = 0; dim < 3; ++dim)
			{
				for (int i = 0; i < 3; ++i)
				{
					build_array[i] = basis[0][i] * translation_vectors[dim][0] + basis[1][i] * translation_vectors[dim][1] + basis[2][i] * translation_vectors[dim][2];
				}
				translation_vectors[dim] = build_array;
			}
			// Calculate NOS
			nos = basis_atoms.size() * n_cells[0] * n_cells[1] * n_cells[2];

			// Spin Positions
			spin_pos = vectorfield(nos);
			Engine::Vectormath::Build_Spins(spin_pos, basis_atoms, translation_vectors, n_cells);
			
			// Log parameters
			Log(Log_Level::Parameter, Log_Sender::IO, "Translation: vectors transformed by basis");
			Log(Log_Level::Parameter, Log_Sender::IO, "        a = " + std::to_string(translation_vectors[0][0]) + " " + std::to_string(translation_vectors[0][1]) + " " + std::to_string(translation_vectors[0][2]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        b = " + std::to_string(translation_vectors[1][0]) + " " + std::to_string(translation_vectors[1][1]) + " " + std::to_string(translation_vectors[1][2]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        c = " + std::to_string(translation_vectors[2][0]) + " " + std::to_string(translation_vectors[2][1]) + " " + std::to_string(translation_vectors[2][2]));
			Log(Log_Level::Parameter, Log_Sender::IO, "Translation: n_cells");
			Log(Log_Level::Parameter, Log_Sender::IO, "        na = " + std::to_string(n_cells[0]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        nb = " + std::to_string(n_cells[1]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        nc = " + std::to_string(n_cells[2]));
			Log(Log_Level::Parameter, Log_Sender::IO, "Geometry: " + std::to_string(nos) + " spins");
			
			// Return geometry
			auto geometry = std::shared_ptr<Data::Geometry>(new Data::Geometry(basis, translation_vectors, n_cells, basis_atoms, lattice_constant, spin_pos));
			Log(Log_Level::Parameter, Log_Sender::IO, "Geometry is " + std::to_string(geometry->dimensionality) + "-dimensional"); 
			Log(Log_Level::Info, Log_Sender::IO, "Geometry: built");
			return geometry;
		}// end Geometry from Config

		std::unique_ptr<Data::Parameters_Method_LLG> Parameters_Method_LLG_from_Config(const std::string configFile)
		{
			//-------------- Insert default values here -----------------------------
			// Output folder for results
			std::string output_folder = "output_llg";
			// Save output when logging
			bool save_output_any = true, save_output_initial = false, save_output_final = true, save_output_energy = true;
			bool save_output_archive = false, save_output_single = false;
			// PRNG Seed
			std::srand(std::time(0));
			int seed = std::rand();
			// number of iterations carried out when pressing "play" or calling "iterate"
			int n_iterations = (int)2E+6;
			// Number of iterations after which the system is logged to file
			int n_iterations_log = 100;
			// Temperature in K
			scalar temperature = 0.0;
			// Damping constant
			scalar damping = 0.5;
			// iteration time step
			scalar dt = 1.0E-02;
			// Whether to renormalize spins after every SD iteration
			bool renorm_sd = 1;
			// spin transfer torque vector
			scalar stt_magnitude = 1.5;
			// spin_current polarisation normal vector
			Vector3 stt_polarisation_normal = { 1.0, -1.0, 0.0 };
			// Force convergence parameter
			scalar force_convergence = 10e-9;

			//------------------------------- Parser --------------------------------
			Log(Log_Level::Info, Log_Sender::IO, "Parameters LLG: building");
			if (configFile != "")
			{
				try {
					IO::Filter_File_Handle myfile(configFile);

					myfile.Read_Single(output_folder, "llg_output_folder");
					myfile.Read_Single(save_output_any, "llg_output_save_any");
					myfile.Read_Single(save_output_initial, "llg_output_save_initial");
					myfile.Read_Single(save_output_final, "llg_output_save_final");
					myfile.Read_Single(save_output_energy, "llg_output_save_energy");
					myfile.Read_Single(save_output_archive, "llg_output_save_archive");
					myfile.Read_Single(save_output_single, "llg_output_save_single");
					myfile.Read_Single(seed, "llg_seed");
					myfile.Read_Single(n_iterations, "llg_n_iterations");
					myfile.Read_Single(n_iterations_log, "llg_n_iterations_log");
					myfile.Read_Single(temperature, "llg_temperature");
					myfile.Read_Single(damping, "llg_damping");
					myfile.Read_Single(dt, "llg_dt");
					// dt = time_step [ps] * 10^-12 * gyromagnetic raio / mu_B  { / (1+damping^2)} <- not implemented
					dt = dt*std::pow(10, -12) / Constants::mu_B*1.760859644*std::pow(10, 11);
					myfile.Read_Single(renorm_sd, "llg_renorm");
					myfile.Read_Single(stt_magnitude, "llg_stt_magnitude");
					myfile.Read_Vector3(stt_polarisation_normal, "llg_stt_polarisation_normal");
					myfile.Read_Single(force_convergence, "llg_force_convergence");
				}// end try
				catch (Exception ex) {
					if (ex == Exception::File_not_Found)
					{
						Log(Log_Level::Error, Log_Sender::IO, "Parameters LLG: Unable to open Config File " + configFile + " Leaving values at default.");
					}
					else throw ex;
				}// end catch
			}
			else Log(Log_Level::Warning, Log_Sender::IO, "Parameters LLG: Using default configuration!");

			// Return
			Log(Log_Level::Parameter, Log_Sender::IO, "Parameters LLG:");
			Log(Log_Level::Parameter, Log_Sender::IO, "        seed                = " + std::to_string(seed));
			Log(Log_Level::Parameter, Log_Sender::IO, "        temperature         = " + std::to_string(temperature));
			Log(Log_Level::Parameter, Log_Sender::IO, "        damping             = " + std::to_string(damping));
			Log(Log_Level::Parameter, Log_Sender::IO, "        time step           = " + std::to_string(dt));
			Log(Log_Level::Parameter, Log_Sender::IO, "        stt magnitude       = " + std::to_string(stt_magnitude));
			Log(Log_Level::Parameter, Log_Sender::IO, "        stt normal          = " + std::to_string(stt_polarisation_normal[0]) + " " + std::to_string(stt_polarisation_normal[1]) + " " + std::to_string(stt_polarisation_normal[2]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        force convergence   = " + std::to_string(force_convergence));
			Log(Log_Level::Parameter, Log_Sender::IO, "        n_iterations        = " + std::to_string(n_iterations));
			Log(Log_Level::Parameter, Log_Sender::IO, "        n_iterations_log    = " + std::to_string(n_iterations_log));
			Log(Log_Level::Parameter, Log_Sender::IO, "        output_folder       = " + output_folder);
			Log(Log_Level::Parameter, Log_Sender::IO, "        save_output_any     = " + std::to_string(save_output_any));
			Log(Log_Level::Parameter, Log_Sender::IO, "        save_output_initial = " + std::to_string(save_output_initial));
			Log(Log_Level::Parameter, Log_Sender::IO, "        save_output_final   = " + std::to_string(save_output_final));
			Log(Log_Level::Parameter, Log_Sender::IO, "        save_output_energy  = " + std::to_string(save_output_energy));
			Log(Log_Level::Parameter, Log_Sender::IO, "        save_output_archive = " + std::to_string(save_output_archive));
			Log(Log_Level::Parameter, Log_Sender::IO, "        save_output_single  = " + std::to_string(save_output_single));
			auto llg_params = std::unique_ptr<Data::Parameters_Method_LLG>(new Data::Parameters_Method_LLG(output_folder, {save_output_any, save_output_initial, save_output_final, save_output_energy, save_output_archive, save_output_single}, force_convergence, n_iterations, n_iterations_log, seed, temperature, damping, dt, renorm_sd, stt_magnitude, stt_polarisation_normal));
			Log(Log_Level::Info, Log_Sender::IO, "Parameters LLG: built");
			return llg_params;
		}// end Parameters_Method_LLG_from_Config

		std::unique_ptr<Data::Parameters_Method_GNEB> Parameters_Method_GNEB_from_Config(const std::string configFile)
		{
			//-------------- Insert default values here -----------------------------
			// Output folder for results
			std::string output_folder = "output_gneb";
			// Save output when logging
			bool save_output_any = true, save_output_initial = false, save_output_final = true, save_output_energy = true;
			// Spring constant
			scalar spring_constant = 1.0;
			// Force convergence parameter
			scalar force_convergence = 10e-9;
			// number of iterations carried out when pressing "play" or calling "iterate"
			int n_iterations = (int)2E+6;
			// Number of iterations after which the system is logged to file
			int n_iterations_log = 100;
			// Number of Energy Interpolation points
			int n_E_interpolations = 10;
			//------------------------------- Parser --------------------------------
			Log(Log_Level::Info, Log_Sender::IO, "Parameters GNEB: building");
			if (configFile != "")
			{
				try {
					IO::Filter_File_Handle myfile(configFile);
					
					myfile.Read_Single(output_folder, "gneb_output_folder");
					myfile.Read_Single(save_output_any, "gneb_output_save_any");
					myfile.Read_Single(save_output_initial, "gneb_output_save_initial");
					myfile.Read_Single(save_output_final, "gneb_output_save_final");
					myfile.Read_Single(save_output_energy, "gneb_output_save_energy");
					myfile.Read_Single(spring_constant, "gneb_spring_constant");
					myfile.Read_Single(force_convergence, "gneb_force_convergence");
					myfile.Read_Single(n_iterations, "gneb_n_iterations");
					myfile.Read_Single(n_iterations_log, "gneb_n_iterations_log");
					myfile.Read_Single(n_E_interpolations, "gneb_n_energy_interpolations");
				}// end try
				catch (Exception ex) {
					if (ex == Exception::File_not_Found)
					{
						Log(Log_Level::Error, Log_Sender::IO, "Parameters GNEB: Unable to open Config File " + configFile + " Leaving values at default.");
					}
					else throw ex;
				}// end catch
			}
			else Log(Log_Level::Warning, Log_Sender::IO, "Parameters GNEB: Using default configuration!");

			// Return
			Log(Log_Level::Parameter, Log_Sender::IO, "Parameters GNEB:");
			Log(Log_Level::Parameter, Log_Sender::IO, "        spring_constant     = " + std::to_string(spring_constant));
			Log(Log_Level::Parameter, Log_Sender::IO, "        force_convergence   = " + std::to_string(force_convergence));
			Log(Log_Level::Parameter, Log_Sender::IO, "        n_E_interpolations  = " + std::to_string(n_E_interpolations));
			Log(Log_Level::Parameter, Log_Sender::IO, "        n_iterations        = " + std::to_string(n_iterations));
			Log(Log_Level::Parameter, Log_Sender::IO, "        n_iterations_log    = " + std::to_string(n_iterations_log));
			Log(Log_Level::Parameter, Log_Sender::IO, "        output_folder       = " + output_folder);
			Log(Log_Level::Parameter, Log_Sender::IO, "        save_output_any     = " + std::to_string(save_output_any));
			Log(Log_Level::Parameter, Log_Sender::IO, "        save_output_initial = " + std::to_string(save_output_initial));
			Log(Log_Level::Parameter, Log_Sender::IO, "        save_output_final   = " + std::to_string(save_output_final));
			Log(Log_Level::Parameter, Log_Sender::IO, "        save_output_energy  = " + std::to_string(save_output_energy));
			auto gneb_params = std::unique_ptr<Data::Parameters_Method_GNEB>(new Data::Parameters_Method_GNEB(output_folder, {save_output_any, save_output_initial, save_output_final, save_output_energy}, force_convergence, n_iterations, n_iterations_log, spring_constant, n_E_interpolations));
			Log(Log_Level::Info, Log_Sender::IO, "Parameters GNEB: built");
			return gneb_params;
		}// end Parameters_Method_LLG_from_Config

		std::unique_ptr<Data::Parameters_Method_MMF> Parameters_Method_MMF_from_Config(const std::string configFile)
		{
			//-------------- Insert default values here -----------------------------
			// Output folder for results
			std::string output_folder = "output_mmf";
			// Save output when logging
			bool save_output_any = true, save_output_initial = false, save_output_final = true, save_output_energy = true;
			// Force convergence parameter
			scalar force_convergence = 10e-9;
			// Number of iterations carried out when pressing "play" or calling "iterate"
			int n_iterations = (int)2E+6;
			// Number of iterations after which the system is logged to file
			int n_iterations_log = 100;
			
			//------------------------------- Parser --------------------------------
			Log(Log_Level::Info, Log_Sender::IO, "Parameters MMF: building");
			if (configFile != "")
			{
				try {
					IO::Filter_File_Handle myfile(configFile);
					
					myfile.Read_Single(output_folder, "mmf_output_folder");
					myfile.Read_Single(save_output_any, "mmf_output_save_any");
					myfile.Read_Single(save_output_initial, "mmf_output_save_initial");
					myfile.Read_Single(save_output_final, "mmf_output_save_final");
					myfile.Read_Single(save_output_energy, "mmf_output_save_energy");
					myfile.Read_Single(force_convergence, "mmf_force_convergence");
					myfile.Read_Single(n_iterations, "mmf_n_iterations");
					myfile.Read_Single(n_iterations_log, "mmf_n_iterations_log");
				}// end try
				catch (Exception ex) {
					if (ex == Exception::File_not_Found)
					{
						Log(Log_Level::Error, Log_Sender::IO, "Parameters MMF: Unable to open Config File " + configFile + " Leaving values at default.");
					}
					else throw ex;
				}// end catch
			}
			else Log(Log_Level::Warning, Log_Sender::IO, "Parameters MMF: Using default configuration!");

			// Return
			Log(Log_Level::Parameter, Log_Sender::IO, "Parameters MMF:");
			Log(Log_Level::Parameter, Log_Sender::IO, "        force_convergence   = " + std::to_string(force_convergence));
			Log(Log_Level::Parameter, Log_Sender::IO, "        n_iterations        = " + std::to_string(n_iterations));
			Log(Log_Level::Parameter, Log_Sender::IO, "        n_iterations_log    = " + std::to_string(n_iterations_log));
			Log(Log_Level::Parameter, Log_Sender::IO, "        output_folder       = " + output_folder);
			Log(Log_Level::Parameter, Log_Sender::IO, "        save_output_any     = " + std::to_string(save_output_any));
			Log(Log_Level::Parameter, Log_Sender::IO, "        save_output_initial = " + std::to_string(save_output_initial));
			Log(Log_Level::Parameter, Log_Sender::IO, "        save_output_final   = " + std::to_string(save_output_final));
			Log(Log_Level::Parameter, Log_Sender::IO, "        save_output_energy  = " + std::to_string(save_output_energy));
			auto mmf_params = std::unique_ptr<Data::Parameters_Method_MMF>(new Data::Parameters_Method_MMF(output_folder, {save_output_any, save_output_initial, save_output_final, save_output_energy}, force_convergence, n_iterations, n_iterations_log));
			Log(Log_Level::Info, Log_Sender::IO, "Parameters MMF: built");
			return mmf_params;
		}

		std::unique_ptr<Engine::Hamiltonian> Hamiltonian_from_Config(const std::string configFile, const std::shared_ptr<Data::Geometry> geometry)
		{
			//-------------- Insert default values here -----------------------------
			// The type of hamiltonian we will use
			std::string hamiltonian_type = "isotropic";

			//------------------------------- Parser --------------------------------
			Log(Log_Level::Info, Log_Sender::IO, "Hamiltonian: building");

			// Hamiltonian type
			if (configFile != "")
			{
				try {
					Log(Log_Level::Info, Log_Sender::IO, "Hamiltonian: deciding type");
					IO::Filter_File_Handle myfile(configFile);

					// What hamiltonian do we use?
					myfile.Read_Single(hamiltonian_type, "hamiltonian");
				}// end try
				catch (Exception ex) {
					if (ex == Exception::File_not_Found) {
						Log(Log_Level::Error, Log_Sender::IO, "Hamiltonian: Unable to open Config File " + configFile + " Using default Hamiltonian: " + hamiltonian_type);
					}
					else throw ex;
				}// end catch
			}
			else Log(Log_Level::Warning, Log_Sender::IO, "Hamiltonian: Using default Hamiltonian: " + hamiltonian_type);
			
			// Hamiltonian
			std::unique_ptr<Engine::Hamiltonian> hamiltonian;
			if (hamiltonian_type == "heisenberg")
			{
				// TODO: to std::move or not to std::move, that is the question...
				hamiltonian = std::move(Hamiltonian_Heisenberg_from_Config(configFile, geometry));
			}
			else if (hamiltonian_type == "gaussian")
			{
				hamiltonian = std::move(Hamiltonian_Gaussian_from_Config(configFile, geometry));
			}
			else
			{
				Log(Log_Level::Error, Log_Sender::IO, "Hamiltonian: Invalid type: " + hamiltonian_type);
			}
			
			// Return
			Log(Log_Level::Info, Log_Sender::IO, "Hamiltonian: built hamiltonian of type: " + hamiltonian_type);
			return hamiltonian;
		}

		
		std::unique_ptr<Engine::Hamiltonian_Heisenberg> Hamiltonian_Heisenberg_from_Config(const std::string configFile, const std::shared_ptr<Data::Geometry> geometry)
		{
			//-------------- Insert default values here -----------------------------
			// Boundary conditions (a, b, c)
			std::vector<int> boundary_conditions_i = { 0, 0, 0 };
			std::vector<bool> boundary_conditions = { false, false, false };
			// Spin moment
			scalarfield mu_s = scalarfield(geometry->nos, 2);	// [nos]
			// External Magnetic Field
			std::string external_field_file = "";
			scalar B = 0;
			Vector3 B_normal = { 0.0, 0.0, 1.0 };
			intfield    external_field_index(geometry->nos);				// [nos]
			scalarfield external_field_magnitude(geometry->nos, 0);	// [nos]
			vectorfield external_field_normal(geometry->nos, B_normal);	// [3][nos]
			
			// Anisotropy
			std::string anisotropy_file = "";
			scalar K = 0;
			Vector3 K_normal = { 0.0, 0.0, 1.0 };
			bool anisotropy_from_file = false;
			intfield    anisotropy_index(geometry->nos);				// [nos]
			scalarfield anisotropy_magnitude(geometry->nos, 0.0);	// [nos]
			vectorfield anisotropy_normal(geometry->nos, K_normal);	// [nos][3]

			// ------------ Pair Interactions ------------
			int n_pairs = 0;
			std::string interaction_pairs_file = "";
			bool interaction_pairs_from_file = false;
			pairfield Exchange_pairs; scalarfield Exchange_magnitude;
			pairfield DMI_pairs; scalarfield DMI_magnitude; vectorfield DMI_normal;
			pairfield DD_pairs; scalarfield DD_magnitude; vectorfield DD_normal;

			scalar dd_radius = 0.0;

			// ------------ Quadruplet Interactions ------------
			int n_quadruplets = 0;
			std::string quadruplets_file = "";
			bool quadruplets_from_file = false;
			quadrupletfield quadruplets; scalarfield quadruplet_magnitude;

			//------------------------------- Parser --------------------------------
			Log(Log_Level::Info, Log_Sender::IO, "Hamiltonian_Heisenberg: building");
			// iteration variables
			int iatom = 0;
			if (configFile != "")
			{
				try {
					IO::Filter_File_Handle myfile(configFile);

					// Boundary conditions
					myfile.Read_3Vector(boundary_conditions_i, "boundary_conditions");
					boundary_conditions[0] = (boundary_conditions_i[0] != 0);
					boundary_conditions[1] = (boundary_conditions_i[1] != 0);
					boundary_conditions[2] = (boundary_conditions_i[2] != 0);

					// Spin moment
					mu_s = scalarfield(geometry->nos, 2.0);
					int N = geometry->n_spins_basic_domain[0];
					if (myfile.Find("mu_s"))
					{
						for (iatom = 0; iatom < N; ++iatom)
						{
							myfile.iss >> mu_s[iatom];
							for (int ispin = 0; ispin < geometry->nos / N; ++ispin)
							{
								mu_s[ispin*N + iatom] = mu_s[iatom];
							}
						}
					}
					else Log(Log_Level::Error, Log_Sender::IO, "Keyword 'mu_s' not found. Using Default: 2.0");

					// External Field
					if (myfile.Find("external_field_file")) myfile.iss >> external_field_file;
					if (external_field_file.length() > 0)
					{
						Log(Log_Level::Warning, Log_Sender::IO, "Hamiltonian_Heisenberg: Read external field file has not been implemented yet. Using 0 field for now.");
						// The file name should be valid so we try to read it
						// Not yet implemented!

						B = external_field_magnitude[0];
						B_normal = external_field_normal[0];
					}
					else 
					{
						// Read parameters from config if available
						myfile.Read_Single(B, "external_field_magnitude");
						myfile.Read_Vector3(B_normal, "external_field_normal");
						B_normal.normalize();

						if (B != 0)
						{
							// Fill the arrays
							for (int i = 0; i < geometry->nos; ++i)
							{
								external_field_index[i] = i;
								external_field_magnitude[i] = B;
								external_field_normal[i] = B_normal;
							}
						}
						else
						{
							external_field_index = intfield(0);
							external_field_magnitude = scalarfield(0);
							external_field_normal = vectorfield(0);
						}
					}

					// Anisotropy
					if (myfile.Find("anisotropy_file")) myfile.iss >> anisotropy_file;
					if (anisotropy_file.length() > 0)
					{
						// The file name should be valid so we try to read it
						Anisotropy_from_File(anisotropy_file, geometry, n_pairs,
							anisotropy_index, anisotropy_magnitude, anisotropy_normal);
						K = anisotropy_magnitude[0];
						K_normal = anisotropy_normal[0];
					}
					else
					{
						// Read parameters from config
						myfile.Read_Single(K, "anisotropy_magnitude");
						myfile.Read_Vector3(K_normal, "anisotropy_normal");
						K_normal.normalize();

						if (K != 0)
						{
							// Fill the arrays
							for (int i = 0; i < geometry->nos; ++i)
							{
								anisotropy_index[i] = i;
								anisotropy_magnitude[i] = K;
								anisotropy_normal[i] = K_normal;
							}
						}
						else
						{
							anisotropy_index = intfield(0);
							anisotropy_magnitude = scalarfield(0);
							anisotropy_normal = vectorfield(0);
						}
					}

					// Interaction Pairs
					if (myfile.Find("interaction_pairs_file")) myfile.iss >> interaction_pairs_file;
					if (interaction_pairs_file.length() > 0)
					{
						// The file name should be valid so we try to read it
						Pairs_from_File(interaction_pairs_file, geometry, n_pairs,
							Exchange_pairs, Exchange_magnitude,
							DMI_pairs, DMI_magnitude, DMI_normal);
					}
					//else
					//{
					//	Log(Log_Level::Warning, Log_Sender::IO, "Hamiltonian_Heisenberg: Default Interaction pairs have not been implemented yet.");
					//	throw Exception::System_not_Initialized;
					//	// Not implemented!
					//}
					
					//		Dipole-Dipole Pairs
					// Dipole Dipole radius
					myfile.Read_Single(dd_radius, "dd_radius");
					// if (dd_radius >0 ) Log(Log_Level::Error, Log_Sender::IO, "Hamiltonian_Heisenberg: Dipole-Dipole energy is not correctly implemented, but you chose a radius > 0! -- r=" + std::to_string(dd_radius));
					// Dipole Dipole neighbours of each spin neigh_dd[nos][max_n]
					// std::vector<std::vector<int>> dd_neigh;
					// // Dipole Dipole neighbour positions of each spin neigh_dd[dim][nos][max_n]
					// std::vector<std::vector<std::vector<scalar>>> dd_neigh_pos;
					// // Dipole Dipole normal vectors [dim][nos][max_n]
					// std::vector<std::vector<std::vector<scalar>>> dd_normal;
					// // Dipole Dipole distance [nos][max_n]
					// std::vector<std::vector<scalar>> dd_distance;
					// // Create the DD neighbours
					// Engine::Neighbours::Create_Dipole_Neighbours(geometry, std::vector<bool>{ true, true, true }, dd_radius, dd_neigh, dd_neigh_pos, dd_normal, dd_distance);
					// // Get the DD pairs from the neighbours
					// Engine::Neighbours::Create_DD_Pairs_from_Neighbours(geometry, dd_neigh, dd_neigh_pos, dd_distance, dd_normal, DD_pairs, DD_magnitude, DD_normal);
					
					
					// Engine::Neighbours::Create_Dipole_Pairs(*geometry, dd_radius, DD_pairs, DD_magnitude, DD_normal);


					// Interaction Quadruplets
					if (myfile.Find("interaction_quadruplets_file")) myfile.iss >> quadruplets_file;
					if (quadruplets_file.length() > 0)
					{
						// The file name should be valid so we try to read it
						Quadruplets_from_File(quadruplets_file, geometry, n_quadruplets,
							quadruplets, quadruplet_magnitude);
					}

				}// end try
				catch (Exception ex) {
					if (ex == Exception::File_not_Found)
					{
						Log(Log_Level::Error, Log_Sender::IO, "Hamiltonian_Heisenberg: Unable to open Config File " + configFile + " Leaving values at default.");
					}
					else throw ex;
				}// end catch
			}
			else Log(Log_Level::Warning, Log_Sender::IO, "Hamiltonian_Heisenberg: Using default configuration!");
			
			// Return
			Log(Log_Level::Parameter, Log_Sender::IO, "Hamiltonian_Heisenberg:");
			Log(Log_Level::Parameter, Log_Sender::IO, "        boundary conditions = " + std::to_string(boundary_conditions[0]) + " " + std::to_string(boundary_conditions[1]) + " " + std::to_string(boundary_conditions[2]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        B[0]                = " + std::to_string(B));
			Log(Log_Level::Parameter, Log_Sender::IO, "        B_normal[0]         = " + std::to_string(B_normal[0]) + " " + std::to_string(B_normal[1]) + " " + std::to_string(B_normal[2]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        mu_s[0]             = " + std::to_string(mu_s[0]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        K[0]                = " + std::to_string(K));
			Log(Log_Level::Parameter, Log_Sender::IO, "        K_normal[0]         = " + std::to_string(K_normal[0]) + " " + std::to_string(K_normal[1]) + " " + std::to_string(K_normal[2]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        dd_radius           = " + std::to_string(dd_radius));
			auto hamiltonian = std::unique_ptr<Engine::Hamiltonian_Heisenberg>(new Engine::Hamiltonian_Heisenberg(
				mu_s,
				external_field_index, external_field_magnitude, external_field_normal,
				anisotropy_index, anisotropy_magnitude, anisotropy_normal,
				Exchange_pairs, Exchange_magnitude,
				DMI_pairs, DMI_magnitude, DMI_normal,
				DD_pairs, DD_magnitude, DD_normal,
				quadruplets, quadruplet_magnitude,
				geometry,
				boundary_conditions
			));
			Log(Log_Level::Info, Log_Sender::IO, "Hamiltonian_Heisenberg: built");
			return hamiltonian;
		}// end Hamiltonian_Heisenberg_From_Config
		
		
		std::unique_ptr<Engine::Hamiltonian_Gaussian> Hamiltonian_Gaussian_from_Config(const std::string configFile, const std::shared_ptr<Data::Geometry> geometry)
		{
			//-------------- Insert default values here -----------------------------
			// Number of Gaussians
			int n_gaussians = 1;
			// Amplitudes
			std::vector<scalar> amplitude = { 1 };
			// Widths
			std::vector<scalar> width = { 1 };
			// Centers
			std::vector<Vector3> center = { Vector3{ 0, 0, 1 } };

			//------------------------------- Parser --------------------------------
			Log(Log_Level::Info, Log_Sender::IO, "Hamiltonian_Gaussian: building");
			
			if (configFile != "")
			{
				try {
					IO::Filter_File_Handle myfile(configFile);

					// N
					myfile.Read_Single(n_gaussians, "n_gaussians");

					// Allocate arrays
					amplitude = std::vector<scalar>(n_gaussians, 1.0);
					width = std::vector<scalar>(n_gaussians, 1.0);
					center = std::vector<Vector3>(n_gaussians, Vector3{0, 0, 1});
					// Read arrays
					if (myfile.Find("gaussians"))
					{
						for (int i = 0; i < n_gaussians; ++i)
						{
							myfile.GetLine();
							myfile.iss >> amplitude[i];
							myfile.iss >> width[i];
							for (int j = 0; j < 3; ++j)
							{
								myfile.iss >> center[i][j];
							}
							center[i].normalize();
						}
					}
					else Log(Log_Level::Error, Log_Sender::IO, "Hamiltonian_Gaussian: Keyword 'gaussians' not found. Using Default: {0, 0, 1}");
				}// end try
				catch (Exception ex) {
					if (ex == Exception::File_not_Found)
					{
						Log(Log_Level::Error, Log_Sender::IO, "Hamiltonian_Gaussian: Unable to open Config File " + configFile + " Leaving values at default.");
					}
					else throw ex;
				}// end catch
			}
			else Log(Log_Level::Warning, Log_Sender::IO, "Hamiltonian_Gaussian: Using default configuration!");


			// Return
			Log(Log_Level::Parameter, Log_Sender::IO, "Hamiltonian_Gaussian:");
			Log(Log_Level::Parameter, Log_Sender::IO, "        n_gaussians  = " + std::to_string(n_gaussians));
			Log(Log_Level::Parameter, Log_Sender::IO, "        amplitude[0] = " + std::to_string(amplitude[0]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        width[0]     = " + std::to_string(width[0]));
			Log(Log_Level::Parameter, Log_Sender::IO, "        center[0]    = " + std::to_string(center[0][0]) + " " + std::to_string(center[0][1]) + " " + std::to_string(center[0][2]));
			auto hamiltonian = std::unique_ptr<Engine::Hamiltonian_Gaussian>(new Engine::Hamiltonian_Gaussian(
				amplitude, width, center
			));
			Log(Log_Level::Info, Log_Sender::IO, "Hamiltonian_Gaussian: built");
			return hamiltonian;
		}
	}// end namespace IO
}// end namespace Utility