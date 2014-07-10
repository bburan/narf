function inc01()
% An incrementally fit model and fitter (intertwined)
% Starts with a simple model, and gradually increases its complexity.
% logn_wcg01_p1z0
% logn_wcgd01_p1z0
% logn_wcgd01_dep1pc_p1z0
% lognbz_wcgd01_dep1pc_p1z0
% lognbz_wcgd01_dep1pc_p1z0_dexp

global MODULES STACK FLATSTACK XXX FLATXXX;

% Build the simple model first
logn();
wcg01();
p3z2();
fit09d(); pop_module(); % Score: 9.643168e-0

% Change the wcg01 to wcgd01
[m_wcg, idx_wcg] = find_modules(STACK, 'weight_channels', 1);
STACK{idx_wcg}{1}.phi   = [STACK{idx_wcg}{1}.phi 0 1.6];
STACK{idx_wcg}{1}.phifn = @wc_gaussian_differences;
XXX = XXX(1:2); FLATXXX = FLATXXX(1:2); calc_xxx(2);  % HACKY RECALC
% fit09d(); pop_module(); % Score: 9.642545e-01

% Add depression before the PZ filter
[pzmod, idx_pz] = find_modules(STACK, 'pole_zeros', 1);
pop_module();
dep1pcw();
append_module(pzmod{1});

% Re-fit with wcgd and depression
fit09d(); pop_module(); %  9.620246e-01,

% Make logn into lognbz, which is pretty nonlinear
[~, idx_log] = find_modules(STACK, 'nonlinearity', 1);
STACK{idx_log}{1}.phi   = [STACK{idx_log}{1}.phi 0 0];
XXX = XXX(1:2); FLATXXX = FLATXXX(1:2); calc_xxx(2);  % HACKY RECALC

% Add the output nonlinearity, which is also pretty nonlinear
dexp();

% Fit everything 
fit09d();