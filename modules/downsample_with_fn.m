function m = downsample_with_fn(args)
% A function to create a (possibly nonlinear) downsampling module.
%
% Very simple algorithm to downsample a signal:
%    Apply PRECONV_FN to every value of the input signal.
%    Take N of those values at a time, where N = ceil(raw_freq / downsampled_freq)
%    Apply those N values to function CONV_FN. (Do the convolution, baby!)
%    Apply POSTCONV_FN to each of those values. 
%
% The only PRECONV_FN function I can think is useful is @abs, but maybe
% somebody will want to use another function.
%
% Good CONV_FN functions include @max, @mean, and @sum, but you may prefer to
% define your own.
%
% Good POSTCONV_FN functions are @sqrt, and maybe @log if your data never
% goes to zero.
%
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information. TODO.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @downsample_with_fn;
m.name = 'downsample_with_fn';
m.fn = @do_downsampling;
m.pretty_name = 'Downsample with Arbitrary Function';
m.editable_fields = {'downsampled_freq', 'preconv_fn', 'conv_fn', 'postconv_fn'};
m.isready_pred = @downsampler_isready;

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_downsampled_stimulus;
m.plot_fns{1}.pretty_name = 'Downsampled Stimulus vs Time';

% Module fields that are specific to THIS MODULE
m.downsampled_freq = 200;
m.preconv_fn = @abs;
m.conv_fn = @mean;
m.postconv_fn = @sqrt;

% Overwrite the default module fields with arguments 
if nargin == 1
    m = merge_structs(m, args);
end

function x = do_downsampling(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    
    ds_dim = 2; % Dimension to downsample on. 
    scale = ceil(baphy_mod.raw_stim_fs / mdl.downsampled_freq);

    for sf = fieldnames(x.dat)', sf=sf{1};
        x.dat.(sf).ds_stim = ...
            mdl.postconv_fn(conv_fn(mdl.preconv_fn(x.dat.(sf).pp_stim), ...
                                    ds_dim, mdl.conv_fn, scale, 0));
        x.dat.(sf).ds_stim_time = ...
            linspace(1/mdl.downsampled_freq, ...
                     x.dat.(sf).raw_stim_time(end), ...
                     length(x.dat.(sf).ds_stim));
    end
end

% Plot the filter responses
function do_plot_downsampled_stimulus(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    % Find the GUI controls
    baphy_mod = find_module(stack, 'load_stim_resps_from_baphy');
    filt_pop = find_module_gui_control(stack, 'selected_filter_popup');
    
    c = cellstr(get(baphy_mod.plot_gui.selected_stimfile_popup, 'String'));
    sf = c{get(baphy_mod.plot_gui.selected_stimfile_popup, 'Value')};
    stim_idx = get(baphy_mod.plot_gui.selected_stim_idx_popup, 'Value');
    filt_idx = get(filt_pop, 'Value');
       
    dat = x.dat.(sf);
    
    plot(dat.ds_stim_time, ...
         squeeze(dat.ds_stim(stim_idx,:,filt_idx)), ...
         pickcolor(filt_idx));
    axis tight;
    drawnow;
end

end