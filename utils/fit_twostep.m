function termcond = fit_twostep()   
    termcond = twostep(@() fit_fminlsq('score', 'fittertempstim', 'respavg'),...
                       @fit_fminlsq);
end
