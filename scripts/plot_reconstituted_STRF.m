function plot_reconstituted_STRF(batch, cellids, modelnames)
global STACK META XXX;

function x = extract_x(x)        
end

[stacks, metas, x0s] = load_model_batch(batch, cellids, modelnames, ...
                       @extract_x, @extract_x, @extract_x);

for ii = 1:length(stacks)    
    load_model(metas{ii}.modelpath);
    calc_xxx(1);
    
    [m, idx] = find_modules(STACK, 'pole_zeros', 1);
    
    if isempty(m)
        [mods, mods_idx] = find_modules(STACK, 'fir_filter');
        if isempty(mods)
            continue;
        end
        m = mods{end};
        idx = mods_idx{end};
    end    
    
    m = m{1};
    
    figure('Name', META.modelname, ...
                'NumberTitle', 'off', 'Position', [20 50 900 900]);
    
    m.auto_plot({}, STACK(1:idx), XXX(1:idx));
    
end
end