function termcond = fit_slsqtwo()    
    termcond = twostep(@fit_slsq, @fit_lsq);
end
