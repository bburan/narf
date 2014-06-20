function m = add_nth_order_terms(args)
% Concatenates nth order terms only in the channel dimension.
% 
% If there are N channels, then nchoosek(N, order) channels will be appended to
% the existing channels when order = 2. Used to model second order spatial
% (but not temporal) multiplicitive interactions between channels. Higher
% order terms can be added by setting order to a higher value.
%

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @add_nth_order_terms;
m.name = 'add_nth_order_terms';
m.fn = @do_add_nth_order_terms;
m.pretty_name = 'Add Nth Order Terms';
m.editable_fields = {'input', 'time', 'output'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input = 'stim'; 
m.time  = 'stim_time';
m.output = 'stim';
m.order = 2;
m.selfterms = false; % Include terms like AA and BB instead of just AB

% Optional fields
m.plot_gui_create_fn = @create_chan_selector_gui;
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_all_default_outputs;
m.plot_fns{1}.pretty_name = 'Channels (All)';
m.plot_fns{2}.fn = @do_plot_single_default_output;
m.plot_fns{2}.pretty_name = 'Channel (Single)';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optimize this module for tree traversal  
m.required = {m.input, m.time};   % Signal dependencies
m.modifies = {m.output};          % These signals are modified

function x = do_add_nth_order_terms(mdl, x)    

    for sf = fieldnames(x.dat)', sf=sf{1};           
        [T, S, C] = size(x.dat.(sf).(mdl.input));
        
        terms = {};
        for o = 2:min(mdl.order, C)
            pairings = nchoosek([1:C], o);
            
            % Add pairings AB, AC, BC, etc
            for ii = 1:size(pairings, 1)
                tmp = num2cell(pairings(ii,:));
                terms{end+1} = tmp;
            end
            
            if isfield(mdl, 'selfterms') && mdl.selfterms
                % Add AA, BB, CC, etc
                for ss = 1:C
                    terms{end+1} = num2cell(ss*ones(o,1));
                end
            end
        end        
        
        % Build and append the nth order terms to the output
        for ii = 1:length(terms),
            p = ones(size(x.dat.(sf).(mdl.input)(:, :, terms{ii}{1})));
            for jj = 1:length(terms{ii})
                p = p .* x.dat.(sf).(mdl.input)(:, :, terms{ii}{jj});
            end
            x.dat.(sf).(mdl.output)(:,:, ii+C) = p;
        end               
        
    end
end

end
