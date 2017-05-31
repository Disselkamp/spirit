import os
import sys

### Make sure to find the Spirit modules
### This is only needed if you did not install the package
# spirit_py_dir = os.path.dirname(os.path.realpath(__file__)) + "core/python/Spirit"
spirit_py_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), "../core/python"))
sys.path.insert(0, spirit_py_dir)


### Import numpy
import numpy as np

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

cfgfile = "input/schieback_displacement_single.cfg"

with state.State(cfgfile) as p_state:

    # configuration.PlusZ(p_state, border_cylindrical=-1)
    # configuration.MinusZ(p_state, border_cylindrical=128.0, pos=[-129.0,0.0,0.0])

    # head to head domain
    configuration.Domain(p_state, [-1.0,0.0,0.0], pos=[0,0,0], border_cylindrical=-1)
    configuration.Domain(p_state, [1.0,0.0,0.0], pos=[-129,0,0], border_cylindrical=128)

    # domain wall
    configuration.Domain(p_state, [-1.0,1.0,0.0], pos=[1,0,0], border_cylindrical=1)
    configuration.Domain(p_state, [-1.0,1.0,0.0], pos=[0,0,0], border_cylindrical=1)
    configuration.Domain(p_state, [0.0,1.0,0.0], pos=[-1,0,0], border_cylindrical=1)
    configuration.Domain(p_state, [1.0,1.0,0.0], pos=[-2,0,0], border_cylindrical=1)
    configuration.Domain(p_state, [-1.0,1.0,0.0], pos=[-3,0,0], border_cylindrical=1)
    

    spins = system.Get_Spin_Directions( p_state, idx_chain=-1 )
    print spins
    print spins[240:270]

    # relaxation
    # i = 0
    # E = [0,1]
    # while E[i%2] < E[(i+1)%2]:
    #         i = i + 1
    #         simulation.PlayPause(p_state, "LLG", "SIB", n_iterations=1000)
    #         E[i%2] = system.Get_Energy(p_state, idx_image=-1, idx_chain=-1)
    simulation.PlayPause(p_state, "LLG", "SIB", n_iterations=10000)
    print "System relaxed"
    
    # STT on 
    hamiltonian.Set_STT(p_state, 0.01, [1.0,0.0,0.0], idx_image=-1, idx_chain=-1)
   
    a = []
    b = []

    for i in range(0,1000):
        simulation.PlayPause(p_state, "LLG", "SIB", n_iterations=1000)
        spins = system.Get_Spin_Directions( p_state, idx_chain=-1 )
        
        print "PlayPause: ", i
        at = []
        
        for j in range(0+2, len(spins)-2):
            if spins[(j-2)][0] > 0.0 and spins[(j-1)][0] > 0.0 and spins[j][0] < 0.0:
                at = at + [j]
                print j
        a = a + [at]

    print b
    print a
