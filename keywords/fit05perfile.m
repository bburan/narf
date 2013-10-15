function fit05perfile()

global STACK XXX

if 1,
    disp('FIRST SINGLE FIT FOR ALL DATA (NO SPLITTING)');
    fit05;
else
    disp('FIRST FIT FIR PER FILE');
    savefitfields=cell(length(STACK),1);
    for ii=1:length(STACK),
        if isfield(STACK{ii}{1},'fit_fields'),
            savefitfields{ii}=STACK{ii}{1}.fit_fields;
        end
    end
    fit05firperfile();
    keyboard
    for ii=1:length(savefitfields),
        if isfield(STACK{ii}{1},'fit_fields'),
            STACK{ii}{1}.fit_fields=savefitfields{ii};
        end
    end
end

% Then boost on each file individually
disp('THEN FIT EACH FILE SEPARATELY');
fit_split_simply(@fit05, @split_by_respfile, @unify_respfiles); 
