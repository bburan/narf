#!/bin/bash 

# A simple script test to launch three instances of matlab

matlab -r "analysis_batch_240(1)" &
matlab -r "analysis_batch_240(2)" &
matlab -r "analysis_batch_240(3)" &