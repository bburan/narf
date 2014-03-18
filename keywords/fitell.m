function fitell()
% Fitting snippet for fitting bandpass filterbanks

wc01(); % Weight channels provides variable gain
nmse();
for ii = 1:6  
    fit_scaat('StopAtAbsScoreDelta', 10^-ii, ...          
              'InitStepSize', 10.0, ...
              'StopAtStepsize', 10^-4);
end
pop_module(); % Remove NMSE
pop_module(); % Remove wc01