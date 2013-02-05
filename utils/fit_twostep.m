function termcond = fit_twostep()    
    termcond = twostep(@fminlsq, @fminlsq);
end
