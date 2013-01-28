function termcond = fit_fminlsq()
    fit_objective();
    termcond = fit_with_lsqcurvefit();
end