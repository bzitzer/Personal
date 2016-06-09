#!/bin/bash

time ./bayes.py --filelist $1                     --Mlisti $2 --dNdE_option $3 2>&1 > TMP/$3_E_$2.dat
