function mixfit3nlperfile()

global STACK XXX

% Relative boost algorithm across all Passive files to initialize
XXX_save=XXX;
[xxxs,unique_codes]=split_by_filecodeP(XXX);
pidx=find(strcmp(unique_codes,'P'));
XXX=xxxs{pidx};
calc_xxx(1);

fitSubstack([],10^-2);

nmse(); % We need a performance metric, so add one 
mixfit3nomse();

% now restore XXX to include active and passive data
XXX=XXX_save;
calc_xxx(2);

[~, mod_idxs] = find_modules(STACK, 'fir_filter', false);
for ii=1:mod_idxs{end}(end),
    if isfield(STACK{ii}{1},'fit_fields'),
        STACK{ii}{1}.fit_fields={};
    end
end

% Then boost on each file individually
fit_split_simply(@mixfit3nomse, @split_by_respfile, @unify_respfiles); 


