function m = split_signal(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @split_signal;
m.name = 'split_signal';
m.fn = @do_split_signal;
m.pretty_name = 'Split Signal';
m.editable_fields = {'splitter', 'input', 'time', 'outputs'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.mapfn = @split_by_filecode; % Also good: unify
m.splitter='';
m.input =  'stim';
m.time =   'stim_time';
m.outputs = {'stim01', 'stim02'};
m.unique_codes = [1 2];

% Optional fields
m.auto_init = @auto_init_split_signal;

m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_channels_as_heatmap;
m.plot_fns{1}.pretty_name = 'Normalized Channels (Heatmap)';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.outputs{:}};      % These signals are modified
% ------------------------------------------------------------------------
% INSTANCE METHODS

function mdl = auto_init_split_signal(stack, xxx)
    % NOTE: Unlike most plot functions, auto_init functions get a 
    % STACK and XXX which do not yet have this module or module's data
    % added to them.
    
    mdl = m;
    
    % figure out how many unique states should be fit from trial_codes
    dat=xxx{end}.dat;
    files=fields(dat);
    input=mdl.input;
    trial_code_all=[];
    for sf=1:length(files),
       if isfield(dat.(files{sf}),'trial_code'),
          trial_code_all=cat(1,trial_code_all,dat.(files{sf}).trial_code);
       end
    end
    mdl.unique_codes=unique(trial_code_all)';
    
    % generate one output per unique trial_code value
    if length(mdl.unique_codes)<=1,
       % only one file code, effectively revery to pass-thru
       mdl.outputs={input};
    else
       mdl.outputs=cell(1,length(mdl.unique_codes));
       for ii=1:length(mdl.unique_codes),
          mdl.outputs{ii}=sprintf('%s%02d',input,mdl.unique_codes(ii));
       end
    end
    mdl.modifies=mdl.outputs;
end

function x = do_split_signal(mdl, x)            
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
         sf=fns{ii};
         for kk=1:length(mdl.unique_codes),
             ff=find(x.dat.(sf).trial_code==mdl.unique_codes(kk));
             x.dat.(sf).(mdl.outputs{kk})=...
                 x.dat.(sf).(mdl.input)(:,ff,:);
         end
    end
end

% 
end
