function fit05firperfile()

global STACK XXX

% Relative boost algorithm across all Passive files to initialize
% $$$ disp('FIRST FIT ONLY ON PASSIVE FILES');
% $$$ XXX_save=XXX;
% $$$ [xxxs,unique_codes]=split_by_filecodeP(XXX);
% $$$ pidx=find(strcmp(unique_codes,'P'));
% $$$ XXX=xxxs{pidx};
% $$$ calc_xxx(1);

disp('FIRST SINGLE FIT FOR ALL DATA (NO SPLITTING)');
fit05();

% now restore XXX to include active and passive data
% $$$ XXX=XXX_save;
% $$$ calc_xxx(2);

disp('NOW FITTING FIR AND LATER STAGES PER FILE');
[~, mod_idxs] = find_modules(STACK, 'fir_filter', false);
for ii=1:(mod_idxs{1}(1)-1), % 1:mod_idxs{end}(end),
    if isfield(STACK{ii}{1},'fit_fields'),
        STACK{ii}{1}.fit_fields={};
    end
end

% Then boost on each file individually
fit_split_simply(@fit05, @split_by_respfile, @unify_respfiles); 
