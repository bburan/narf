function m = truncate_data(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @truncate_data;
m.name = 'truncate_data';
m.fn = @do_truncate_data;
m.pretty_name = 'Weight Channels';
m.editable_fields = {'save_fraction', 'truncate_reps', ...
                    'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.weights = [1]; % Each column weights several channels to produce
m.y_offset = [0]; % A column of y-offsets to be added to each output chan. 
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields
m.is_splittable = true;
m.auto_plot = @do_plot_channels_as_heatmap;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_channels_as_heatmap;
m.plot_fns{1}.pretty_name = 'All Channels (Heatmap)';
m.plot_fns{1}.fn = @do_plot_all_default_outputs;
m.plot_fns{1}.pretty_name = 'Output Channels (All)';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% ------------------------------------------------------------------------
% INSTANCE METHODS

function x = do_truncate_data(mdl, x, stack, xxx)   
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
             savetimebins=round(T.*save_fraction);
             
             x.dat.(sf).(fields{f})=tmp(1:savetimebins,:,:);
         end
    end
end

% ------------------------------------------------------------------------
% Plot methods

% use defaults


end