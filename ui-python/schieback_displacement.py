import os
import sys

### Make sure to find the Spirit modules
### This is only needed if you did not install the package
# spirit_py_dir = os.path.dirname(os.path.realpath(__file__)) + "core/python/Spirit"
spirit_py_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), "../core/python"))
sys.path.insert(0, spirit_py_dir)

### Import numpy
import numpy as np
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

cfgfile = "input/input.cfg"

with state.State(cfgfile) as p_state:
    
    filename = "output/spins_schieback_relaxed_VP.txt"
    # WELCHE WERTE FUER STTMAGNITUDE?
    k = 0
    for stt_magnitude in [0.01, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008, 0.009, 0.01, 0.015, 0.02, 0.025, 0.03, 0.032, 0.034, 0.036, 0.038, 0.04]:
        
        # head to head domain with domain wall - relaxed with VP
        print "read relaxed System" # verschoben um -1 von Mitte
        io.Image_Read(p_state, filename, fileformat=0, idx_image=-1, idx_chain=-1)

        print "STT_magnitude: %.3f"%round(stt_magnitude, 3)
        hamiltonian.Set_STT(p_state, stt_magnitude, [1.0,0.0,0.0], idx_image=-1, idx_chain=-1)

        ts = time.time()
        filename = 'output/'+datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d_%H-%M-%S')+'_Wallpositions_STTmagn%.3f'%stt_magnitude
        f = open(filename, 'w')
        f.write('Step    Position')
        measurement = [range(0,10), range(0,1000), range(0,1000), range(0,1000), range(0,1000), range(0,1000), range(0,600), range(0,600), range(0,600), range(0,600), range(0,400), range(0,400), range(0,400), range(0,400), range(0,400), range(0,400), range(0,400), range(0,400), range(0,400)]
        for i in measurement[k]:
            simulation.PlayPause(p_state, "LLG", "SIB", n_iterations=1000)
            spins = system.Get_Spin_Directions( p_state, idx_chain=-1)

            print "PlayPause: ", i
            
            for j in range(0+2, len(spins)-2):
                if spins[(j-2)][0] > 0.0 and spins[(j-1)][0] > 0.0 and spins[j][0] < 0.0:
                    print j
                    f.write('\n'+str(i)+'    '+str(j))
        k += 1
