#!/bin/sh -l
python2.7 ui-python/bobber_sktube.py "bobber-paper-L14_0.2" 30 -1 0.0001 0.0003 0.0005 > ./sktube_1_20180319.log 2>&1 &
python2.7 ui-python/bobber_sktube.py "bobber-paper-L14_0.2" 30 -1 0.0006 0.0008 0.0010 > ./sktube_2_20180319.log 2>&1 &
python2.7 ui-python/bobber_sktube.py "bobber-paper-L14_0.2" 30 -1 0.003 0.005 0.007 > ./sktube_3_20180319.log 2>&1 &
python2.7 ui-python/bobber_sktube.py "bobber-paper-L14_0.2" 30 -1 0.009 0.008 0.009 > ./sktube_4_20180319.log 2>&1 &