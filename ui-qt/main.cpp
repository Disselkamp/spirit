#include "MainWindow.hpp"

#include "Spirit/State.h"
#include "Spirit/Chain.h"
#include "Spirit/Configurations.h"
#include "Spirit/Transitions.h"
#include "Spirit/Log.h"

// Main
int main(int argc, char ** argv)
{
	//--- Register SigInt
	// signal(SIGINT, Signal::Handle_SigInt);
	
	//---------------------- file names ---------------------------------------------
	//--- Config Files
	 const char * cfgfile = "input/input.cfg";
	// const char * cfgfile = "input/anisotropic/markus.cfg";
	// const char * cfgfile = "input/anisotropic/markus-paper.cfg";
	// const char * cfgfile = "input/anisotropic/kagome-spin-ice.cfg";
	// const char * cfgfile = "input/anisotropic/gideon-master-thesis-anisotropic.cfg";
	// const char * cfgfile = "input/isotropic/gideon-master-thesis-isotropic.cfg";
	// const char * cfgfile = "input/isotropic/daniel-master-thesis-isotropic.cfg";
	// const char * cfgfile = "input/gaussian/example-1.cfg";
	// const char * cfgfile = "input/gaussian/gideon-paper.cfg";
	//--- Data Files
	// std::string spinsfile = "input/anisotropic/achiral.txt";
	// std::string chainfile = "input/chain.txt";
	//-------------------------------------------------------------------------------
	
	//--- Initialise State
	std::shared_ptr<State> state = std::shared_ptr<State>(State_Setup(cfgfile), State_Delete);

	//---------------------- initialize spin_systems --------------------------------
	// Copy the system a few times
	/*Chain_Image_to_Clipboard(state.get());
	for (int i=1; i<7; ++i)
	{
		Chain_Insert_Image_After(state.get());
	}*/
	//-------------------------------------------------------------------------------
	
	//----------------------- spin_system_chain -------------------------------------
	// Read Image from file
	//Configuration_from_File(state.get(), spinsfile, 0);
	// Read Chain from file
	//Chain_from_File(state.get(), chainfile);

	// First image is homogeneous with a Skyrmion at pos
	Configuration_PlusZ(state.get());
	Configuration_Skyrmion(state.get(), 6.0, 1.0, -90.0, false, false, false);
	// Last image is homogeneous
	Chain_Jump_To_Image(state.get(), Chain_Get_NOI(state.get())-1);
	Configuration_PlusZ(state.get());
	Chain_Jump_To_Image(state.get(), 0);

	// Create transition of images between first and last
	Transition_Homogeneous(state.get(), 0, Chain_Get_NOI(state.get())-1);

	// Update the Chain's Data'
	Chain_Update_Data(state.get());
	//-------------------------------------------------------------------------------

	//------------------------ User Interface ---------------------------------------
	// Initialise Application and MainWindow
	QApplication app(argc, argv);
	//app.setOrganizationName("--");
	//app.setApplicationName("Spirit - Atomistic Spin Code - OpenGL with Qt");

	// Format for all GL Surfaces
	QSurfaceFormat format;
	format.setSamples(16);
	format.setVersion(3, 3);
	//format.setVersion(4, 2);
	//glFormat.setVersion( 3, 3 );
	//glFormat.setProfile( QGLFormat::CoreProfile ); // Requires >=Qt-4.8.0
	//glFormat.setSampleBuffers( true );
	format.setProfile(QSurfaceFormat::CoreProfile);
	format.setDepthBufferSize(24);
	format.setStencilBufferSize(8);
	QSurfaceFormat::setDefaultFormat(format);
	Log_Send(state.get(), Log_Level_Info, Log_Sender_UI, ("QSurfaceFormat version: " + std::to_string(format.majorVersion()) + "." + std::to_string(format.minorVersion())).c_str());

	MainWindow window(state);
	window.setWindowTitle(app.applicationName());
	window.show();
	// Open the Application
	int exec = app.exec();
	// If Application is closed normally
	if (exec != 0) throw exec;
	// Finish
	return exec;
	//-------------------------------------------------------------------------------


	return 0;
}
