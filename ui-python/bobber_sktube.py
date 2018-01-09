import os
import sys

### Make sure to find the Spirit modules
### This is only needed if you did not install the package
# spirit_py_dir = os.path.dirname(os.path.realpath(__file__)) + "core/python/Spirit"
spirit_py_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), "../core/python"))
sys.path.insert(0, spirit_py_dir)

### Import numpy
import numpy as np
np.set_printoptions(threshold=np.nan)

### others
import datetime
import time

### Import Spirit modules
from spirit import state, system, geometry, chain, configuration, transition, simulation, quantities, io, log, hamiltonian, parameters


if len(sys.argv) < 4: sys.exit("execute with: 1) beta (eg. 0.00); 2) cfgfile (eg. \"12x12x7_b0.00\" has to be in /input); 3)N (number of Spins NxNxL); 4)borderspherical (eg. 3 (-1 = SkyrmionTube))  5) stt_magnitudes (eg. 0.025 0.050 0.100 0.150)")  # abort if not enough arguments are given
beta = float(sys.argv[1])
cfgfile = sys.argv[2]
N = int(sys.argv[3])
borderspherical = float(sys.argv[4])

with state.State("input/"+cfgfile+".cfg") as p_state:

    for stt_magnitude in [ float(i) for i in sys.argv[5:] ]:  # read stt magnitude from terminal input

        # Plus Z
        configuration.PlusZ(p_state, pos=[0.0,0.0,0.0], border_rectangular=[-1.0,-1.0,-1.0], border_cylindrical=-1.0, border_spherical=-1.0, inverted=False, idx_image=-1, idx_chain=-1)
        print("PlusZ")

        # Add Noise
        configuration.Add_Noise_Temperature(p_state, 5, pos=[0,0,0], border_rectangular=[-1,-1,-1], 
                          border_cylindrical=-1, border_spherical=-1, inverted=False, 
                          idx_image=-1, idx_chain=-1)
        print("Noise added")

        # deploy skyrmiontube/bobber
        radius = 5  # radius of skyrmiontube/bobber
        configuration.Skyrmion(p_state, radius, order=1, phase=0, upDown=False, achiral=False, rightleft=True, pos=[0,0,4], border_rectangular=[-1,-1,-1], border_cylindrical=-1, border_spherical=borderspherical, inverted=False, idx_image=-1, idx_chain=-1)  # position has to change according to system size
        print("Bobber/SkyrmionT initiated")

        np.savetxt("bobber_configuration.txt",system.Get_Spin_Directions(p_state, idx_chain=-1))

        # relax the system with SIB (no skyrmion/bobber with VP)
        simulation.PlayPause(p_state, "LLG", "VP", n_iterations=16000)
        print("Bobber/SkyrmionT relaxed")

        np.savetxt("relaxed_configuration.txt",system.Get_Spin_Directions(p_state, idx_chain=-1))

        # set current in x-direction (True = gradient method)
        parameters.llg.setSTT(p_state, True, stt_magnitude, [1.0,0.0,0.0], idx_image=-1, idx_chain=-1)
        print("STT_magnitude: %.3f"%round(stt_magnitude, 3))

        # create directory for output files
        directory = "output/"+cfgfile+"/beta_%.2f/stt_%.3f"%(beta, stt_magnitude)
        if not os.path.exists(directory): os.makedirs(directory)

        # file for position data
        filename1 = directory+"/"+'positions_STTmagn%.3f'%stt_magnitude
        file_positions = open(filename1, 'w')
        file_positions.write('Step    Position')


        steps = 500
        for i in range(0, steps):
            
            spins = system.Get_Spin_Directions(p_state, idx_chain=-1)  # spin directions
            m = min([spin[2] for spin in spins[-N*N:]])  # find spin with smallest z-component (of spins lying in the +z-surface 12*12=144)

            for j in range(0, len(spins[-N*N:])):  # search and save position of the skyrmiontube/bobber
                if spins[j-N*N][2] == m:
                    print(j)
                    file_positions.write('\n'+str(i)+'    '+str(j))
            
            if i%10 == 0 or i == steps-1:  # save spin directions every 10 steps and the last
                filename2 = directory+"/"+'Spins'+str(i)+'.txt'
                np.savetxt(filename2, spins)

            simulation.PlayPause(p_state, "LLG", "SIB", n_iterations=1000)  # calculate 1000 iterations
            print("PlayPause: ", i)
