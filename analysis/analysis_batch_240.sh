#!/bin/bash 

# A simple script test to launch four instances of matlab

matlab -r "analysis_batch_240(1,4)" &
matlab -r "analysis_batch_240(2,4)" &
matlab -r "analysis_batch_240(3,4)" &
matlab -r "analysis_batch_240(4,4)" &