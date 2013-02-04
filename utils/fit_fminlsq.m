function termcond = fit_fminlsq()
% TODO: What about non-default field names for the optimization?
% This function needs arguments!
    fit_objective();
    termcond = fit_lsq();
end
