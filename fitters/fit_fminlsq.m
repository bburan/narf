function [termcond, n_iters] = fit_fminlsq()
[~, n_iters1] = fit_fminsearch();
[termcond, n_iters2] = fit_lsq();
n_iters = n_iters1 + n_iters2;
