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
gmx gyrate -f traj_comp.xtc -s topol.tpr
xmgrace gyrate.xvg
../tools/xvg_trim.sh gyrate.xvg > rg.xvg
gmx analyze -f rg.xvg -dist hist-rg.xvg -bw 0.025
xmgrace hist-rg.xvg
../tools/pdf2pmf.py -f hist-rg.xvg -T 300 --units kJ/mol > pmf-rg.xvg
xmgrace pmf-rg.xvg

#== end-to-end distance (and free energy)
gmx distance -f traj_comp.xtc -s topol.tpr -n index-pairs.ndx -oall distance.xvg
../tools/xvg_trim.sh distance.xvg > dist.xvg
gmx analyze -f dist.xvg -dist hist-dist.xvg -bw 0.1
xmgrace hist-dist.xvg
../tools/pdf2pmf.py -f hist-dist.xvg -T 300 --units kJ/mol > pmf-dist.xvg
xmgrace pmf-dist.xvg
../tools/xvg_zip.sh rg.xvg dist.xvg > xy-rg-dist.xvg
xmgrace xy-rg-dist.xvg

../tools/xvg_zip.sh rg.xvg energy-potential.xvg > xy-rg-potential.xvg
xmgrace xy-rg-potential.xvg
#== zip together rg and dist into a 2D plot
#== cleanup
rm \#*

exit 0
