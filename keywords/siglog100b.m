function siglog100b()

global MODULES XXX;

m_respavg = [];
for sf = fieldnames(XXX{end}.dat)', sf=sf{1};
    m_respavg = [m_respavg; mean(XXX{end}.dat.(sf).respavg,2)];
end
p1 = prctile(m_respavg,5);
p2 = prctile(m_respavg,95);
p3 = median(m_respavg); % would be better if computed with stim
p4 = 100;
p5 = 100;

append_module(MODULES.nonlinearity.mdl(...
    struct('fit_fields', {{'phi'}}, ...
            'fit_constraints', {{struct( ...
                                     'var', 'phi', ...
                                     'lower', [0     0   0   1   1], ...
                                     'upper', [20 180 200 500 500], ...
                                     'fitter', 'custom')}}, ...
           'phi', [p1 p2 p3 p4 p5], ...
           'nlfn', @nl_sig_logistic100b)));

fitSubstack();

