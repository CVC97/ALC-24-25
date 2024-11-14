#!/usr/bin/env python3
#
#   reads: an xvg file with columns x, p(x)
#   writes: an xvg file with cols x, G(x)
#
#

import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-f', type=str, required=True, help='xvg file for a probability distribution [x p(x)]')
parser.add_argument('-T', type=float, default=None, help='temperature [K]')
parser.add_argument('--units', type=str, default='kBT', choices=['kBT', 'kJ/mol'], help='units: if kJ/mol, then -T is required')
parser.add_argument('--keepinf', action='store_true', help='keep inf values of free energy')
args = parser.parse_args()


if args.units == 'kJ/mol':
  kBT = 0.008315 * args.T
  unitstr = 'kJ/mol'
else: 
  kBT = 1.0
  unitstr = 'k_BT'
  
ylbl = 'Free Energy [%s]' % unitstr 


xp = []
#== printing header, and parsing numeric data
with open(args.f) as ff:
  for line in ff.readlines():
    if line[0] in ['#','@']:
      if line.find('yaxis')>0 and ylbl is not None :
        print ('@   yaxis  label  "%s"' %( ylbl ))
      else:
        print (line,)
    else:
      xp.append(line.split())

#== [x p(x)]
xp = np.asarray(xp, dtype=float)

#== drop p(x)=0,  to avoid log(0) warnings?
if not args.keepinf:
  keep = [ xp[:,1] != 0. ][0]
  xp = xp[keep,:]

#== free energy
G = -np.log(xp[:,1])*kBT
G -= min(G)
#== new array
xg = xp*1.0
xg[:,1] = G

#== print numerical data
for row in xg:
    print ('%5g   %5g'% tuple(row))
  


