function m = weight_channels(args)

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @weight_channels;
m.name = 'weight_channels';
m.fn = @do_weight_channels;
m.pretty_name = 'Weight Channels';
m.editable_fields = {'weights', 'y_offset', 'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.weights = [1]; % Each column weights several channels to produce
m.y_offset = [0]; % A column of y-offsets to be added to each output chan. 
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

% ------------------------------------------------------------------------
% INSTANCE METHODS

function x = do_weight_channels(mdl, x, stack, xxx)   
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
         sf=fns{ii};
         
         % Check dimensions are OK
         [T, S, C] = size(x.dat.(sf).(mdl.input));
         if ~isequal(C, size(mdl.weights, 1))
            error('Dimensions of (mdl.input) don''t match weights.');
         end

         % Sum the weighted values of each stream to create new channels
         tmp = zeros(T, S, size(mdl.weights, 2));
         for s = 1:S
             in = squeeze(x.dat.(sf).(mdl.input)(:, s, :));
             tmp(:,s,:) = in * mdl.weights;
         end         
         x.dat.(sf).(mdl.output) = bsxfun(@plus, tmp', mdl.y_offset)';
    end
end

% ------------------------------------------------------------------------
% Plot methods
% N/A

end