function m = file_masker(args)
% Module which REMOVES ALL files from the training set, test set, and XXX{}.dat
% structure, unless their indexes, as positioned in the training set, 
% are explicitly included in only_indexes.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @file_masker;
m.name = 'file_masker';
m.fn = @do_file_masker;
m.pretty_name = 'Mask Files';
m.editable_fields = {'inputs', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.inputs = {'stim1', 'stim2'}; 
m.time  = 'stim_time';
m.only_indexes = []; % Default masks out everything
m.output = 'stim';

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @(stack, xxx) do_plot_signal(stack, xxx, m.time, m.output);
m.plot_fns{1}.pretty_name = 'Output vs Time';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

function x = do_file_masker(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Remove entries indexed in only_indexes
    fns = fieldnames(x.dat);
    for ii = 1:length(fns), 
        sf = fns{ii};
        if ~isempty(find(ii == mdl.only_indexes, 1))
            x.dat = rmfield(x.dat, sf); 
        end
    end
    
    x.training_set = x.training_set(mdl.only_indexes);
    x.test_set = x.test_set(mdl.only_indexes);
    
end

end