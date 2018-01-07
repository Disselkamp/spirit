#!/bin/sh -l
module unload gcc
module load gcc/6
python ui-python/bobber_sktube.py 0.02 "bobber-paper-L12" 30 3 0.030 0.035 0.040 0.045 0.050 0.055 0.060 0.065 0.070 0.075 0.080 0.085 0.090 0.100
