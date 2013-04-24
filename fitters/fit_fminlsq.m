function termcond = fit_fminlsq(objective_score, field1, field2)
if nargin < 3
    objective_score = 'score';
    field1 = 'stim';
    field2 = 'respavg';
end
    fit_fminsearch(objective_score);
    termcond = fit_lsq(field1, field2);
end
