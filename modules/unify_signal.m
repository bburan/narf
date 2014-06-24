function m = unify_signal(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @unify_signal;
m.name = 'unify_signal';
m.fn = @do_unify_signal;
m.pretty_name = 'Unify Signal';
m.editable_fields = {'unifier','inputs', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.unifier='';
m.time =   'stim_time';
m.inputs = {'stim01', 'stim02'};
m.output = 'stim';
m.unique_codes = [1 2];

% Optional fields
m.auto_init = @auto_init_unify_signal;

m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_channels_as_heatmap;
m.plot_fns{1}.pretty_name = 'Normalized Channels (Heatmap)';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.inputs{:}, m.time};   % Signal dependencies
m.modifies = {m.output};      % These signals are modified
% ------------------------------------------------------------------------
% INSTANCE METHODS

function mdl = auto_init_unify_signal(stack, xxx)
% NOTE: Unlike most plot functions, auto_init functions get a 
% STACK and XXX which do not yet have this module or module's data
% added to them.
    
    mdl = m;
    
    % figure out how many unique states should be fit from trial_codes
    dat=xxx{end}.dat;
    files=fields(dat);
    trial_code_all=[];
    for sf=1:length(files),
        if isfield(dat.(files{sf}),'trial_code'),
            trial_code_all=cat(1,trial_code_all,dat.(files{sf}).trial_code);
        end
    end
    mdl.unique_codes=unique(trial_code_all)';
    
    % generate one output per unique trial_code value
    output=mdl.output;
    if length(mdl.unique_codes)<=1,
        % only one file code, effectively revery to pass-thru
        mdl.inputs={output};
    else
        mdl.inputs=cell(1,length(mdl.unique_codes));
        for ii=1:length(mdl.unique_codes),
            mdl.inputs{ii}=sprintf('%s%02d',output,mdl.unique_codes(ii));
        end
    end
    mdl.required={mdl.inputs{:}, mdl.time};
    mdl.modifies={mdl.output};
end

function x = do_unify_signal(mdl, x)            
    fns = fieldnames(x.dat);
    
    for ii = 1:length(fns)
        sf=fns{ii};
                
        % Find the size of the input signals
        for jj = 1:length(mdl.inputs)
           insize=size(x.dat.(sf).(mdl.inputs{jj}));    
           if insize(2) ~= 0 
                break;
           end        
        end 
        
        if size(x.dat.(sf).(mdl.inputs{jj}), 3) > 1;
            fprintf('Insize weird.\n');            
            keyboard;
        end
        x.dat.(sf).(mdl.output)=zeros([insize(1),...
                                      length(x.dat.(sf).trial_code),...
                                      insize(3:end)]);
               
        for kk=1:length(mdl.unique_codes),
            ff = find(x.dat.(sf).trial_code == mdl.unique_codes(kk));
            if ~isempty(ff),
               x.dat.(sf).(mdl.output)(:,ff,:) = x.dat.(sf).(mdl.inputs{kk});
            end
        end
    end
end

end
