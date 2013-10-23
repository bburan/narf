function fit05nlperfile()

global STACK XXX

% Relative boost algorithm across all Passive files to initialize
% $$$ disp('FIRST FIT ONLY ON PASSIVE FILES');
% $$$ XXX_save=XXX;
% $$$ [xxxs,unique_codes]=split_by_filecodeP(XXX);
% $$$ pidx=find(strcmp(unique_codes,'P'));
% $$$ XXX=xxxs{pidx};
% $$$ calc_xxx(1);
disp('FIRST FIT ACROSS ALL FILES');

fit05();

% now restore XXX to include active and passive data
% $$$ XXX=XXX_save;
% $$$ calc_xxx(2);

disp('NOW FITTING POST-FIR STAGES PER FILE');
[~, mod_idxs] = find_modules(STACK, 'fir_filter', false);
save_fitfields=cell(length(STACK),1);
for ii=1:(mod_idxs{end}(1)), % 1:mod_idxs{end}(end),
    if isfield(STACK{ii}{1},'fit_fields'),
        save_fitfields{ii}=STACK{ii}{1}.fit_fields;
        STACK{ii}{1}.fit_fields={};
    end
end

% Then boost on each file individually
fit_split_simply(@fit05, @split_by_filecode, @unify_respfiles); 
%fit_split_simply(@fit05a, @split_by_filecode, @unify_respfiles); 

% restore fit fields settings
for ii=1:length(save_fitfields),
    if ~isempty(save_fitfields{ii}),
        for jj=1:length(STACK{ii}),
            STACK{ii}{jj}.fit_fields=save_fitfields{ii};
        end
    end
end
