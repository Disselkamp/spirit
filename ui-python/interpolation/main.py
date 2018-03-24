import os
import sys

import scipy
from scipy import interpolate, ndimage
import matplotlib.pyplot as plt

### Make sure to find the Spirit modules
### This is only needed if you did not install the package
# spirit_py_dir = os.path.dirname(os.path.realpath(__file__)) + "core/python/Spirit"
spirit_py_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), "../../core/python"))
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
from spirit import parameters

npointsx=3000
npointsy=3000

cfgfile = "bobber-paper-L14_0.2.cfg"

p_state=state.setup(cfgfile)

#io.Chain_Read(p_state, 'Spins0.txt')
#chain.Update_Data(p_state)
io.Image_Read(p_state, 'Spins0.txt')

noi = chain.Get_NOI(p_state)

bmin,bmax = geometry.Get_Bounds(p_state)

grid_x, grid_y = np.mgrid[bmin[0]:bmax[0]:npointsx*1j, bmin[1]:bmax[1]:npointsy*1j]
# grid_x, grid_y = np.mgrid[0:30:npointsx*1j, 0:30:npointsy*1j]

positions = geometry.Get_Spin_Positions(p_state)
points = positions[-30*30:,0:2]

allmin = []

# energy = chain.Get_Energy(p_state)
# rx = chain.Get_Rx(p_state)

for i in range(0, noi):
    val = system.Get_Spin_Directions(p_state, idx_image=i)[-30*30:]
    values = val[:,2]
    grid_z = interpolate.griddata(points, values, (grid_x, grid_y), method='cubic')
    minima = (grid_z == ndimage.minimum_filter(grid_z, 8))
    valx,valy = np.nonzero(minima)
    for ki in range(valx.shape[0]-1,-1,-1):
	if (grid_z[valx[ki],valy[ki]] > -0.8):
            valx = np.delete(valx, ki)
            valy = np.delete(valy, ki)
    plt.imshow(grid_z.T, extent=(bmin[0],bmax[0],bmin[1],bmax[1]), origin='lower', cmap="RdBu")
    # plt.imshow(minima.T, extent=(bmin[0],bmax[0],bmin[1],bmax[1]), origin='lower')
    plt.plot(grid_x[valx,valy],grid_y[valx,valy],'o', color="white")
    print grid_x[valx,valy], grid_y[valx,valy]
    plt.show()
    allmin.append([])
    for ki in range(0,valx.shape[0]):
        allmin[i].append([grid_x[valx[ki],valy[ki]],grid_y[valx[ki],valy[ki]]])

allmin = np.array(allmin)

# for i in range(1, noi):
#     if (np.linalg.norm(allmin[i][0]-allmin[i-1][0]) > np.linalg.norm(allmin[i][0]-allmin[i-1][1])):
#         allmin[i] = np.array([allmin[i][1],allmin[i][0]])

# phi = np.zeros((noi))
# for i in range(0,noi):
#     dist=allmin[i][1]-allmin[i][0]
#     phi[i] = np.arctan2(dist[1],dist[0])

# plt.plot([1,2,3])
# plt.subplot(211)
# plt.plot(phi,energy)
# plt.subplot(212)
# plt.plot(rx,energy)
# plt.show()
