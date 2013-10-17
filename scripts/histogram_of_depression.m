function histogram_of_depression(batch, cellids, modelnames)

if length(modelnames) > 1
	error('It does not make sense to compare depressions from different model structures. Try again!');
end

    function ret = extract_tau(stack)
        mod = find_modules(stack, 'depression_filter_bank', true);
        if isempty(mod)
            ret = nan;
        else
            ret = mod{1}.tau; 
            if ret > 50 || ret < 0
                ret = NaN;
            end
        end
    end

[taus, ~, ~] = load_model_batch(batch, cellids, modelnames, @extract_tau);

taus = cell2mat(taus);
taus = taus(~isnan(taus));

figure('Name', 'Histogram of Depression', 'NumberTitle', 'off', 'Position', [20 50 900 900]);
hist(taus, 20); 
% set(h, 'Interpreter', 'none');
xlabel('Synaptic Depression Time Constant [ms]');
ylabel('Occurances');
title('Synaptic Depression Time Constant Histogram');

end