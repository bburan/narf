function fit09k()

global STACK;

semse();

phi = pack_fittables(STACK)'; 

opts = saoptimset('MaxIter', length(phi) * 100, ...
                  'MaxFunEvals', length(phi) * 100, ...
                  'InitialTemperature', 1, ...
                  'Display', 'diagnose', ...
                  'TolFun', 0.01);
fit_anneal(opts);

end
