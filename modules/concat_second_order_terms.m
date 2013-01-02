function m = concat_second_order_terms(args)
% Concatenates second order terms only in the channel dimension.
% 
% If there are N channels, then nchoosek(N, 2) channels will be appended to
% the existing channels. This can be used to model second order spatial
% (but not temporal) multiplicitive interactions between channels.
%

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @concat_second_order_terms;
m.name = 'concat_second_order_terms';
m.fn = @do_concat_second_order_terms;
m.pretty_name = 'Concat Second Order Terms';
m.editable_fields = {'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input = 'stim'; 
m.time  = 'stim_time';
m.output = 'stim';

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @(stack, xxx) do_plot_channel_vs_time(stack, xxx, m.time, m.output);
m.plot_fns{1}.pretty_name = 'Channel vs Time';
m.plot_fns{2}.fn = @(stack, xxx) do_plot_channels_as_heatmap(stack, xxx, m.output);
m.plot_fns{2}.pretty_name = 'Channel vs Time (Heatmap)';
m.plot_gui_create_fn = @create_chan_selector_gui;

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

function x = do_concat_second_order_terms(stack, xxx)
    mdl = stack{end};
    x = xxx{end};

    for sf = fieldnames(x.dat)', sf=sf{1};           
        [T, S, C] = size(x.dat.(sf).(mdl.input));
        
        pairs = nchoosek([1:C], 2);
        n_pairs = nchoosek(C, 2);
        
        % Build the second order terms
        for ii = 1:n_pairs;
            a = x.dat.(sf).(mdl.input)(:, :, pairs(1));
            b = x.dat.(sf).(mdl.input)(:, :, pairs(2));
            x.dat.(sf).(mdl.output)(:,:, ii+C) = a .* b;
        end
    end
end

end