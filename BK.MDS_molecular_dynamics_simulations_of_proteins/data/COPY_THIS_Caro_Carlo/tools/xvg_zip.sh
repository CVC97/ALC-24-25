#/bin/bash



egrep  -v "@|#" $1 | awk '{print $2}' > aaaa.dat
egrep  -v "@|#" $2 | awk '{print $2}' > bbbb.dat
paste aaaa.dat bbbb.dat && rm aaaa.dat bbbb.dat

