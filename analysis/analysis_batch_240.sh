#!/bin/bash 

# A simple script test to launch multiple instances of matlab
THISPATH=`pwd`/..
BAPHYPATH=/home/ivar/matlab/baphy

xterm -fg white -bg black -geometry 80x24+0+0     -e matlab $BAPHYPATH $THISPATH -nodesktop -r "analysis_batch_240(1,6)" &
xterm -fg white -bg black -geometry 80x24+600+0   -e matlab $BAPHYPATH $THISPATH -nodesktop -r "analysis_batch_240(2,6)" &
xterm -fg white -bg black -geometry 80x24+0+400   -e matlab $BAPHYPATH $THISPATH -nodesktop -r "analysis_batch_240(3,6)" &
xterm -fg white -bg black -geometry 80x24+600+400 -e matlab $BAPHYPATH $THISPATH -nodesktop -r "analysis_batch_240(4,6)" &
xterm -fg white -bg black -geometry 80x24+0+800   -e matlab $BAPHYPATH $THISPATH -nodesktop -r "analysis_batch_240(5,6)" &
xterm -fg white -bg black -geometry 80x24+600+800 -e matlab $BAPHYPATH $THISPATH -nodesktop -r "analysis_batch_240(6,6)" &