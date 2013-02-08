#!/bin/bash 

# A simple script test to launch multiple instances of matlab
THISPATH=`pwd`/..
BAPHYPATH=/home/ivar/matlab/baphy

MATLABPATH=$THISPATH:$BAPHYPATH
export MATLABPATH

xterm -fg white -bg black -geometry 80x24+0+0     -e matlab -nodesktop -r "analysis_batch_241(1,6)" &
xterm -fg white -bg black -geometry 80x24+600+0   -e matlab -nodesktop -r "analysis_batch_241(2,6)" &
xterm -fg white -bg black -geometry 80x24+0+400   -e matlab -nodesktop -r "analysis_batch_241(3,6)" &
xterm -fg white -bg black -geometry 80x24+600+400 -e matlab -nodesktop -r "analysis_batch_241(4,6)" &
xterm -fg white -bg black -geometry 80x24+0+800   -e matlab -nodesktop -r "analysis_batch_241(5,6)" &
xterm -fg white -bg black -geometry 80x24+600+800 -e matlab -nodesktop -r "analysis_batch_241(6,6)" &