#!/bin/bash 

# A script to launch multiple instances of matlab and start them analyzing
THISPATH=`pwd`/..
BAPHYPATH=/home/ivar/matlab/baphy

MATLABPATH=$THISPATH:$BAPHYPATH
export MATLABPATH

xterm -fg white -bg black -geometry 80x24+0+0     -e matlab -nodesktop -r "analyze_batches([240, 241, 242],1,6)" &
xterm -fg white -bg black -geometry 80x24+600+0   -e matlab -nodesktop -r "analyze_batches([240, 241, 242],2,6)" &
xterm -fg white -bg black -geometry 80x24+0+400   -e matlab -nodesktop -r "analyze_batches([240, 241, 242],3,6)" &
xterm -fg white -bg black -geometry 80x24+600+400 -e matlab -nodesktop -r "analyze_batches([240, 241, 242],4,6)" &
xterm -fg white -bg black -geometry 80x24+0+800   -e matlab -nodesktop -r "analyze_batches([240, 241, 242],5,6)" &
xterm -fg white -bg black -geometry 80x24+600+800 -e matlab -nodesktop -r "analyze_batches([240, 241, 242],6,6)" &