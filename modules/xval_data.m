function m = split_signal(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @split_signal;
m.name = 'split_signal';
m.fn = @do_split_signal;
m.pretty_name = 'Split Signal';
m.editable_fields = {'splitter', 'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.output};          % These signals are modified
% ------------------------------------------------------------------------
% INSTANCE METHODS

function x = do_split_signal(mdl, x)   
    fns = x.training_set;
    save_fraction=mdl.save_fraction;
    if save_fraction>1,save_fraction=1;end
    if save_fraction<0.01,save_fraction=0.01;end
    
    for ii = 1:length(fns)
         sf=fns{ii};
         fields=fieldnames(x.dat.(sf));
         for f=1:length(fields),
             tmp=x.dat.(sf).(fields{f});
             [T, S, C] = size(tmp);
             %savetimebins=round(T.*save_fraction);
             if S>1,
                savestimcount=round(S.*save_fraction);
             
                x.dat.(sf).(fields{f})=tmp(:,1:savestimcount,:);
             end
         end
    end
end

% ------------------------------------------------------------------------
% Plot methods

% use defaults
unique_codes = unique(filecodes);
xxxs = cell(1, length(unique_codes));
estfiles = xxx{end}.training_set;
valfiles = xxx{end}.test_set;

% Otherwise, group by filecode
for ii = 1:length(unique_codes);
    fc = unique_codes{ii};
    
    matches = strcmp(fc, filecodes);   
    
    efs = estfiles(matches);
    vfs = valfiles(matches(matches<=length(valfiles))); 
    
    xxxs{ii} = xxx;
    xxxs{ii}{end}.dat = [];   
    xxxs{ii}{end}.training_set={};
    xxxs{ii}{end}.test_set={};
    xxxs{ii}{end}.filecodes={};
    
    for jj = 1:length(efs)
        f = efs{jj};
        xxxs{ii}{end}.dat.(f) = xxx{end}.dat.(f);
        xxxs{ii}{end}.training_set{jj} = f;
        xxxs{ii}{end}.filecodes{jj} = fc;
    end
    for jj = 1:length(vfs)
        f = vfs{jj};
        xxxs{ii}{end}.dat.(f) = xxx{end}.dat.(f);
        xxxs{ii}{end}.test_set{jj} = f;
    end
    
end




end
