#!/bin/bash
#
#   md analysis script
#



#== this allows you to input varaibles into the script
echo "input the simulation temperature [K] and a nametag separated by a space:"
read temp tag
echo "# your input paramters:" 
echo "temperature=$temp"
echo "tag=$tag"
echo
echo "press enter to continue (ctrl-c to abort)"
read asdf




#== run some g tools
echo "0" | gmx trjconv -f traj_comp.xtc -s topol.tpr -dt 10 -conect -o traj-$tag.pdb
echo "5 0" | gmx energy -f ener.edr -o energy-potential.xvg
echo "2 0" | gmx energy -f ener.edr -o energy-angle.xvg
 
 
#== radius of gyration (and free energy)


#== end-to-end distance (and free energy)


#== zip together rg and dist into a 2D plot
#== cleanup
rm \#*

exit 0
