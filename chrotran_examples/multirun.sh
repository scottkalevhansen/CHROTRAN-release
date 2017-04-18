#!/bin/bash

./runScript.sh -t=test1 -f=test1 -np=1
./runScript.sh -t=test1 -f=test1_numerical -np=1
./runScript.sh -t=test2-column -f=test2smp -np=1

./runScript.sh -t=test8 -f=test8_fakeParams -np=16
./runScript.sh -t=test8 -f=test8_fakeParams_retard -np=16
./runScript.sh -t=test8 -f=test8_fakeParams_retard_direct -np=16

./runScript.sh -t=test9 -f=test9_noretard_nodirect -np=16
./runScript.sh -t=test9 -f=test9_retard_nodirect -np=16
./runScript.sh -t=test9 -f=test9_noretard_direct -np=16
./runScript.sh -t=test9 -f=test9_retard_direct -np=16

./runScript.sh -t=test10 -f=test10_biocide -np=16
./runScript.sh -t=test10 -f=test10_alcohol -np=16
./runScript.sh -t=test10 -f=test10_alcohol_biocide -np=16
