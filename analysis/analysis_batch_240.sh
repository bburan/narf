#!/bin/bash 

# A simple script test to launch four instances of matlab

xterm -fg white -bg black -geometry 80x24+0+0     -e matlab -nodesktop -r "analysis_batch_240(1,4)" &
xterm -fg white -bg black -geometry 80x24+600+0   -e matlab -nodesktop -r "analysis_batch_240(2,4)" &
xterm -fg white -bg black -geometry 80x24+0+400   -e matlab -nodesktop -r "analysis_batch_240(3,4)" &
xterm -fg white -bg black -geometry 80x24+600+400 -e matlab -nodesktop -r "analysis_batch_240(4,4)" &
