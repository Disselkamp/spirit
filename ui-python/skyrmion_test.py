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
import datetime
import time

### Import Spirit modules
from spirit import state
from spirit import system
from spirit import geometry
from spirit import chain
from spirit import configuration
from spirit import transition
from spirit import simulation
from spirit import quantities
from spirit import io
from spirit import log
from spirit import hamiltonian

cfgfile = "input/50x50.cfg" # beta??

with state.State(cfgfile) as p_state:
    
    # WELCHE WERTE FUER STTMAGNITUDE?
    k = 0
    for stt_magnitude in [0.036]: #, 0.038, 0.040, 0.042, 0.044, 0.046, 0.050, 0.055]:
        
        configuration.PlusZ(p_state, pos=[0.0,0.0,0.0], border_rectangular=[-1.0,-1.0,-1.0], border_cylindrical=-1.0, border_spherical=-1.0, inverted=False, idx_image=-1, idx_chain=-1)
        print "PlusZ"

        radius = 5
        configuration.Skyrmion(p_state, radius, order=1, phase=0, upDown=True, achiral=False, rightleft=False, pos=[0,0,0], border_rectangular=[-1,-1,-1], border_cylindrical=-1, border_spherical=-1, inverted=False, idx_image=-1, idx_chain=-1)
        print "Skyrmion initiated"

        simulation.PlayPause(p_state, "LLG", "SIB", n_iterations=1000)
        print "Skyrmion relaxed"

        print "STT_magnitude: %.3f"%round(stt_magnitude, 3)
        hamiltonian.Set_STT(p_state, stt_magnitude, [1.0,0.0,0.0], idx_image=-1, idx_chain=-1)

        directory = "output/%.3f"%stt_magnitude
        if not os.path.exists(directory):
            os.makedirs(directory)

        filename1 = directory+"/"+'positions_STTmagn%.3f'%stt_magnitude
        f = open(filename1, 'w')
        f.write('Step    Position')

        steps = 500

        for i in range(0, steps):
            simulation.PlayPause(p_state, "LLG", "SIB", n_iterations=1000)
            spins = system.Get_Spin_Directions( p_state, idx_chain=-1)

            if i%10 == 0:
                filename2 = directory+"/"+'Spins'+str(i)+'.txt'
                np.savetxt(filename2, spins)

            print "PlayPause: ", i

            m = min([i[2] for i in spins])
            for j in range(0, len(spins)):
                if spins[j][2] == m:
                    print j
                    f.write('\n'+str(i)+'    '+str(j))
        k += 1
