#!/bin/sh

cd /move_to_the_current_directory/

inputpk=./Planck2015_massiveneu_matterpower.dat 
Om=0.3156 
h0=0.6727 
ns=0.96
hod=cmass_param_1404_3742
z=0.55

./calc_tinker07 test.dat $inputpk $Om $h0 $ns $hod $z
