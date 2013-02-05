function termcond = fit_twostep()    
    termcond = twostep(@fit_fminlsq, @fit_fminlsq);
end
