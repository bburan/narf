function m = file_masker(args)
% Module which REMOVES ALL files from the training set, test set, and XXX{}.dat
% structure, unless their filecodes match exactly what is in only_filecodes.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @file_masker;
m.name = 'file_masker';
m.fn = @do_file_masker;
m.pretty_name = 'Mask Files';
m.editable_fields = {'only_filecode'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.only_filecode = ''; 

% Optional fields
m.plot_fns = {};

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

function x = do_file_masker(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Remove an entry unless it's filecode match only_filecode   
    for ii = 1:length(x.filecodes),    
        if ~strcmp(mdl.only_filecode, x.filecodes{ii})
            x.dat = rmfield(x.dat, x.training_set{ii});      
            x.dat = rmfield(x.dat, x.test_set{ii});       
        end
    end
    
    msk = strcmp(x.filecodes, mdl.only_filecode);
    if any(msk)
        x.training_set = x.training_set(msk);
        x.test_set = x.test_set(msk);
        x.filecodes = x.filecodes(msk);
    end
    
end

end