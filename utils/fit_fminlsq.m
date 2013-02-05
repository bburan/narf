function termcond = fit_fminlsq()
% TODO: What about non-default field names for the optimization?
% This function needs arguments!
    fit_fminsearch();
    termcond = fit_lsq();
end
