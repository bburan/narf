function fit05firperfile()

global STACK XXX

disp('FIRST SINGLE FIT FOR ALL DATA (NO SPLITTING)');
fit05();

disp('NOW FITTING FIR AND LATER STAGES PER FILE');
[~, mod_idxs] = find_modules(STACK, 'fir_filter', false);
save_fitfields=cell(length(STACK),1);
for ii=1:(mod_idxs{1}(1)-1), % 1:mod_idxs{end}(end),
    if isfield(STACK{ii}{1},'fit_fields'),
        save_fitfields{ii}=STACK{ii}{1}.fit_fields;
        STACK{ii}{1}.fit_fields={};
    end
end

% Then boost on each file individually
fit_split_simply(@fit05, @split_by_filecode, @unify_respfiles); 
%fit_split_simply(@fit00, @split_by_filecode, @unify_respfiles); 

% restore fit fields settings
for ii=1:length(save_fitfields),
    if ~isempty(save_fitfields{ii}),
        for jj=1:length(STACK{ii}),
            STACK{ii}{jj}.fit_fields=save_fitfields{ii};
        end
    end
end

