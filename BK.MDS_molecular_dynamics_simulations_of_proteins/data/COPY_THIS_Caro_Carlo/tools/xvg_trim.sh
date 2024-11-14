#!/bin/bash



egrep     "@|#" $1 
egrep  -v "@|#" $1 | awk '{print $1,$2}'



