function m = weight_channels_with_gaussians(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @weight_channels_with_gaussians;
m.name = 'weight_channels_with_gaussians';
m.fn = @do_weight_channels_with_gaussians;
m.pretty_name = 'Weight Channels (Gaussians)';
m.editable_fields = {'mu', 'sigma', 'y_offset', 'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.mu      = [1]; % Gaussian centers in units of CHANNELS
m.sigma   = [0]; % Std. Dev in units of CHANNELS
m.y_offset = [0];
m.input =  'stim';
m.time =   'stim_time';
m.output = 'stim';

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_all_default_outputs;
m.plot_fns{1}.pretty_name = 'Output Channels (All)';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.output};          % These signals are modified

% ------------------------------------------------------------------------
% INSTANCE METHODS

function x = do_weight_channels_with_gaussians(mdl, x)   
    fns = fieldnames(x.dat);
       
    weights = [];
    [T, S, C] = size(x.dat.(fns{1}).(mdl.input));
    for ii = 1:length(mdl.mu)
       weights(:, ii) = gauss1([mdl.mu(ii) mdl.sigma(ii)], [1:C]);
    end
                 
    for ii = 1:length(fns)
         sf=fns{ii};
         
         % Check dimensions are OK
         [T, S, C] = size(x.dat.(sf).(mdl.input));
         if ~isequal(C, size(weights, 1))
            error('Dimensions of (mdl.input) don''t match weights.');
         end

         % Sum the weighted values of each stream to create new channels
         tmp = zeros(T, S, size(weights, 2));
         for s = 1:S
             in = squeeze(x.dat.(sf).(mdl.input)(:, s, :));
             tmp(:,s,:) = in * weights;
         end         
         x.dat.(sf).(mdl.output) = tmp + mdl.y_offset;
    end
end

% ------------------------------------------------------------------------
% Plot methods
% N/A

end
