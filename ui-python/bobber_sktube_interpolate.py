import os
import sys

import scipy
from scipy import interpolate, ndimage

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


if len(sys.argv) < 3: sys.exit("execute with: 1) cfgfile (eg. \"12x12x7_b0.00\" has to be in /input); 2) N (number of Spins NxNxL); 3) borderspherical (eg. 3 (-1 = SkyrmionTube))  4) stt_magnitudes (eg. 0.025 0.050 0.100 0.150)")  # abort if not enough arguments are given
# beta = float(sys.argv[1])
cfgfile = sys.argv[1]
N = int(sys.argv[2])
borderspherical = float(sys.argv[3])

# get beta from cfg file
for line in open("input/"+cfgfile+".cfg"):
    if 'llg_beta' in line:
        beta = float(line.split('\t')[-1])

# save used cfg file to folder
directory0 = "output/"+cfgfile
if not os.path.exists(directory0): os.makedirs(directory0)
open(directory0+"/"+cfgfile+".cfg", 'w').write(open("input/"+cfgfile+".cfg").read())

with state.State("input/"+cfgfile+".cfg") as p_state:

    for stt_magnitude in [ float(i) for i in sys.argv[4:] ]:  # read stt magnitude from terminal input

        # create directory for output files
        directory = "output/"+cfgfile+"/beta_%.2f/stt_%.4f"%(beta, stt_magnitude)
        if not os.path.exists(directory): os.makedirs(directory)

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

        np.savetxt(directory+"/bobber_configuration.txt",system.Get_Spin_Directions(p_state, idx_chain=-1))

        # relax the system with SIB (no skyrmion/bobber with VP)
        simulation.PlayPause(p_state, "LLG", "VP", n_iterations=16000)
        print("Bobber/SkyrmionT relaxed")

        np.savetxt(directory+"/relaxed_configuration.txt",system.Get_Spin_Directions(p_state, idx_chain=-1))

         # set current in x-direction (True = gradient method)
        parameters.llg.setSTT(p_state, True, stt_magnitude, [1.0,0.0,0.0], idx_image=-1, idx_chain=-1)
        print("STT_magnitude: %.4f"%round(stt_magnitude, 4))

        # file for position data
        filename1 = directory+"/"+'positions_STTmagn%.4f'%stt_magnitude
        file_positions = open(filename1, 'w')
        file_positions.write('x-position'+'\t'+'y-position')
        file_positions.close()

        npointsx=600
        npointsy=600
        bmin,bmax = geometry.Get_Bounds(p_state)
        grid_x, grid_y = np.mgrid[bmin[0]:bmax[0]:npointsx*1j, bmin[1]:bmax[1]:npointsy*1j]
        positions = geometry.Get_Spin_Positions(p_state)
        points = positions[-N*N:,0:2]

        steps = 20000
        for i in range(0, steps):

            spins = system.Get_Spin_Directions(p_state, idx_image=-1)[-N*N:]
            spins_z = spins[:,2]
            grid_z = interpolate.griddata(points, spins_z, (grid_x, grid_y), method='cubic')
            minima = (grid_z == ndimage.minimum_filter(grid_z, 8))
            valx,valy = np.nonzero(minima)

	    zt = 0
            for ki in range(valx.shape[0]-1,-1,-1):
                if (grid_z[valx[ki],valy[ki]] > -0.8):
                        valx = np.delete(valx, ki)
                        valy = np.delete(valy, ki)
		elif (grid_z[valx[ki], valy[ki]] < zt):
			zt = grid_z[valx[ki], valy[ki]]
			xvalue = grid_x[valx[ki],valy[ki]]
			yvalue = grid_y[valx[ki],valy[ki]]

            # plt.imshow(grid_z.T, extent=(bmin[0],bmax[0],bmin[1],bmax[1]), origin='lower', cmap="RdBu")
            # plt.plot(grid_x[valx,valy],grid_y[valx,valy],'o', color="white")
            # plt.show()

            file_positions = open(filename1, 'a')
            file_positions.write('\n'+str(xvalue)+'\t'+str(yvalue))
            file_positions.close()

            if i%2000 == 0 or i == steps-1:  # save spin directions every 10 steps and the last
                filename2 = directory+"/"+'Spins'+str(i)+'.txt'
                np.savetxt(filename2, spins)

            simulation.PlayPause(p_state, "LLG", "SIB", n_iterations=1000)  # calculate 1000 iterations
            print("PlayPause: ", i)
            print(str(grid_x[valx,valy])+'    '+str(grid_y[valx,valy]))
