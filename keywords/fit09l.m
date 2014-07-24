function fit09l()

global STACK;

semse();

phi = pack_fittables(STACK)'; 

opts = saoptimset('MaxIter', length(phi) * 100, ...
                  'MaxFunEvals', length(phi) * 100, ...
                  'InitialTemperature', 5, ...
                  'Display', 'diagnose', ...
                  'TolFun', 0.01);
fit_anneal(opts);

end
