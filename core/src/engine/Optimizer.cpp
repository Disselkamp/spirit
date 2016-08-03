#include "Optimizer.h"
#include "Timing.h"
#include"Logging.h"

using namespace Utility;


namespace Engine
{
    Optimizer::Optimizer(std::vector<std::shared_ptr<Data::Spin_System>> systems, std::shared_ptr<Engine::Method> method)
    {
        this->systems = systems;

        this->noi = systems.size();
        this->nos = systems[0]->nos;

        this->method = method;

        // this->configurations = std::vector<std::vector<double>>(this->noi, std::vector<double>(3 * this->nos));
        // for (int i = 0; i < this->noi; ++i)
        // {
        //     this->configurations[i] = systems[i]->spins;
        // }

        // this->force = std::vector<std::vector<double>>(this->noi, std::vector<double>(3 * this->nos, 0));	// [noi][3*nos]

        // this->force_call = force_call;
        // // Calculate forces once, so that the Solver does not think it's converged
        // this->force_call->Calculate(this->configurations, this->force);

        // Setup Timings
        for (int i=0; i<7; ++i) this->t_iterations.push_back(system_clock::now());
        this->ips = 0;
        this->starttime = Timing::CurrentDateTime();
    }

    
    void Optimizer::Iterate()
    {
        // TODO: get this from Method
        auto sender = Log_Sender::GNEB;

		//------------------------ Init local vars ---------------------------------
		int n_iterations    = this->method->parameters->n_iterations;
		int log_steps       = this->method->parameters->log_steps;
		this->starttime     = Timing::CurrentDateTime();
        //----
		int i, step = 0, n_log = n_iterations/log_steps;
        // TODO: How to get the right suffix?
		std::string suffix = ""; // GNEB
		// std::string suffix = "_archive"; // LLG
		//------------------------ End Init ----------------------------------------

        // Log message
		Log.Send(Log_Level::ALL, sender, "-------------- Started " + this->method->Name() + " Simulation --------------");
		Log.Send(Log_Level::ALL, sender, "Going to iterate " + std::to_string(n_log) + " steps");
        Log.Send(Log_Level::ALL, sender, "            with " + std::to_string(log_steps) + " iterations per step");
		Log.Send(Log_Level::ALL, sender, "Optimizer: " + this->FullName());
		Log.Send(Log_Level::ALL, sender, "-----------------------------------------------------");

        // Start Timings
		auto t_start = system_clock::now();
		auto t_current = system_clock::now();
		auto t_last = system_clock::now();

        // Iteration loop
		for (i = 0; i < n_iterations && this->ContinueIterating(); ++i)
		{
            // Pre-Iteration hook
            this->method->Hook_Pre_Step();
			// Do one single Iteration
			this->Iteration();
            // Post-Iteration hook
            this->method->Hook_Pre_Step();

			// Recalculate FPS
			this->t_iterations.pop_front();
			this->t_iterations.push_back(system_clock::now());

			// Log Output every log_steps steps
			if (0 == fmod(i, log_steps))
			{
				++step;

				t_last = t_current;
				t_current = system_clock::now();

				Log.Send(Log_Level::ALL, sender, this->Name() + " Iteration step          " + std::to_string(step) + " / " + std::to_string(n_log));
				Log.Send(Log_Level::ALL, sender, "                           = " + std::to_string(i) + " / " + std::to_string(n_iterations));
				Log.Send(Log_Level::ALL, sender, "    Time since last step:    " + std::to_string(Timing::SecondsPassed(t_last, t_current)) + " seconds.");
				Log.Send(Log_Level::ALL, sender, "    Iterations / sec:        " + std::to_string(log_steps / Timing::SecondsPassed(t_last, t_current)));
				Log.Send(Log_Level::ALL, sender, "    Maximum force component: " + std::to_string(this->method->force_maxAbsComponent));

				this->method->Save_Step(0, i, suffix);

				//output_strings[step - 1] = IO::Spins_to_String(c->images[0].get());
			}// endif log_steps
		}// endif i

        // End timing
		auto t_end = system_clock::now();

		Log.Send(Log_Level::ALL, sender, "-------------- Finished " + this->method->Name() + " Simulation --------------");
		Log.Send(Log_Level::ALL, sender, "Terminated at                   " + std::to_string(i) + " / " + std::to_string(n_iterations) + " iterations.");
		if (this->method->Force_Converged())
			Log.Send(Log_Level::ALL, sender, "    The transition has converged to a maximum force component of " + std::to_string(this->method->force_maxAbsComponent));
		else
			Log.Send(Log_Level::ALL, sender, "    Maximum force component:    " + std::to_string(this->method->force_maxAbsComponent));
        Log.Send(Log_Level::ALL, sender, "    Force convergence parameter: " + std::to_string(this->method->parameters->force_convergence));
		if (this->StopFilePresent())
			Log.Send(Log_Level::ALL, sender, "    A STOP file has been found.");
		Log.Send(Log_Level::ALL, sender, "    " + this->method->Name() + " Simulation ran for     " + std::to_string(Timing::MinutesPassed(t_start, t_end)) + " minutes.");
		Log.Send(Log_Level::ALL, sender, "Optimizer: " + this->FullName());
		Log.Send(Log_Level::ALL, sender, "------------------------------------------------------");

        // TODO: How to get the right suffix?
		suffix = "_final"; // GNEB
        // suffix = "_" + IO::int_to_formatted_string(i, (int)log10(n)) + "_final"; // LLG
		this->method->Save_Step(0, i, suffix);
		//IO::Dump_to_File(output_strings, "spin_archieve.dat", c->images[0]->debug_parameters->output_notification, step);
    }

    
    void Optimizer::Iteration()
    {
        // Not Implemented!
        Log.Send(Log_Level::L_ERROR, Log_Sender::ALL, std::string("Tried to use Optimizer::Step() of the Optimizer base class!"));
    }


    bool Optimizer::ContinueIterating()
    {
        return !this->method->Force_Converged() && !this->StopFilePresent(); // && c->iteration_allowed;
    }


    double Optimizer::getIterationsPerSecond()
    {
        double l_ips = 0.0;
        for (unsigned int i = 0; i < t_iterations.size() - 1; ++i)
        {
            l_ips += Timing::SecondsPassed(t_iterations[i], t_iterations[i+1]);
        }
        this->ips = 1.0 / (l_ips / (t_iterations.size() - 1));
        return this->ips;
    }

    bool Optimizer::StopFilePresent()
    {
        std::ifstream f("STOP");
        return f.good();
    }

    // Optimizer name as string
    std::string Optimizer::Name()
    {
        // Not Implemented!
        Log.Send(Log_Level::L_ERROR, Log_Sender::ALL, std::string("Tried to use Optimizer::Name() of the Optimizer base class!"));
        return "--";
    }

    std::string Optimizer::FullName()
    {
        // Not Implemented!
        Log.Send(Log_Level::L_ERROR, Log_Sender::ALL, std::string("Tried to use Optimizer::Fullname() of the Optimizer base class!"));
        return "--";
    }
}