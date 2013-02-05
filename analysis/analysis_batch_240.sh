#!/bin/bash 

# A simple script test to launch four instances of matlab

konsole --workdir ./ -e matlab -nodesktop -r "analysis_batch_240(1,4)" &
konsole --workdir ./ -e matlab -nodesktop -r "analysis_batch_240(2,4)" &
konsole --workdir ./ -e matlab -nodesktop -r "analysis_batch_240(3,4)" &
konsole --workdir ./ -e matlab -nodesktop -r "analysis_batch_240(4,4)" &
