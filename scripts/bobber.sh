#!/bin/sh -l
python2.7 ui-python/bobber_sktube.py "bobber-paper-L14_0.2" 30 -1 0.001 0.002 0.003 > ./bobber_1_20180319.log 2>&1 &
python2.7 ui-python/bobber_sktube.py "bobber-paper-L14_0.2" 30 -1 0.004 0.005 0.006 > ./bobber_2_20180319.log 2>&1 &
python2.7 ui-python/bobber_sktube.py "bobber-paper-L14_0.2" 30 -1 0.007 0.008 > ./bobber_3_20180319.log 2>&1 &
python2.7 ui-python/bobber_sktube.py "bobber-paper-L14_0.2" 30 -1 0.009 0.010 > ./bobber_4_20180319.log 2>&1 &