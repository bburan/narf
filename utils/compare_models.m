function fig = compare_models(filenames) 
% Given a cell array of file names, display a dense summary of information
% about each one in a scrollable window.  Not very general yet; right now,
% it simply plots:
%    1. The default plot of the module in STACK{2}
%    2. The first FIR filter coefs that it finds.
%    3. The default plot of the module immediately after the FIR filter
% Remember that filenames must be the FULL PATH TO THE FILES.
% Returns the figure handle to the window displayed.

global STACK XXX NARF_SAVED_MODELS_PATH;

ax = cell(size(filenames));

% Size parameters of GUI
w = 1200;  
h = 120;
y0 = 10;  % For scrolling
vspace = 0.2; % Relative
hspace = 0.05; 
h_win = 1200;

fig = figure('Menubar','figure', 'Resize','off', ...
             'doublebuffer', 'on', 'Units','pixels', ...
             'Position', [50 50 w+20 h_win]);
         
function ret = text_pos (index)
    ret = [w*0.00 y0+h*(index-1) w*0.20 (1-vspace)*h];
end

function ret = compression_pos (index)
    ret = [w*0.25 y0+h*(index-1) w*0.20 (1-vspace)*h];
end

function ret = fircoefs_pos (index)
    ret = [w*0.50 y0+h*(index-1) w*0.20 (1-vspace)*h];
end

function ret = nonlinearity_pos (index)
    ret = [w*0.75 y0+h*(index-1) w*0.20 (1-vspace)*h];
end

score = [];

% For each model, 
for ii = 1:length(filenames)
    
    % load the model STACK and XXX{1} and extract needed info
    f = filenames{ii};
    load_model(f);
    f = regexprep(f, [NARF_SAVED_MODELS_PATH filesep], ''); % Remove redundant prefix just for printing
    
    recalc_xxx(1);
    [firmod, firmod_idx] = find_module(STACK, 'fir_filter');
    
    % Buld some plot axes for compression, coefficients, and nonlinearity
    ax{ii}.compr = axes('Parent', fig, 'Units','pixels', 'Position', compression_pos(ii));
    ax{ii}.coefs = axes('Parent', fig, 'Units','pixels', 'Position', fircoefs_pos(ii));
    ax{ii}.nlin  = axes('Parent', fig, 'Units','pixels', 'Position', nonlinearity_pos(ii));
    ax{ii}.text  = uicontrol('Parent', fig, 'Style', 'text', ...
            'Units', 'pixels', 'FontName', 'FixedWidth', ...
            'HorizontalAlignment', 'left', ...
            'Position', text_pos(ii));
    
    % Use the STACK's functions to plot prettily
    axes(ax{ii}.compr);
    if strcmp(STACK{2}.name, 'nonlinearity')
        STACK{2}.plot_fns{3}.fn(STACK(1:2), XXX(1:3));
    end
    axes(ax{ii}.coefs);
    firmod.plot_fns{1}.fn(STACK(1:firmod_idx), XXX(1:firmod_idx+1));
    axes(ax{ii}.nlin);
    if strcmp(STACK{firmod_idx+1}.name, 'nonlinearity') | ...
       strcmp(STACK{firmod_idx+1}.name, 'sparse_empirical_nonlinearity') | ...
       strcmp(STACK{firmod_idx+1}.name, 'gmm_nonlinearity') | ...    
       strcmp(STACK{firmod_idx+1}.name, 'nonparm_nonlinearity')     
        STACK{firmod_idx+1}.plot_fns{1}.fn(STACK(1:firmod_idx+1), XXX(1:firmod_idx+2));
    end
    
    % Compute the moment of the distribution
    firmod = find_module(STACK, 'fir_filter');
    sparsity = sparsity_metric(firmod.coefs);
    
    % Write to the text box some interesting values
    set(ax{ii}.text, 'String', ...
        sprintf('%s\nTrain r^2: %.5f\nTest r^2:  %.5f\nSparsity: %.5f', ...
                 f, ...
                 XXX{end}.score_train_corr, ...
                 XXX{end}.score_test_corr, ...
                 sparsity));
    
    score(ii,1) = XXX{end}.score_test_corr;
    score(ii,2) = ii;
    score(ii,3) = XXX{end}.score_train_corr;             
    
    % At the end of the day, the only thing that knows how to print a model
    % properly is _itself_. However, a model is not a distinct object, it
    % is a collection of modules. So unless these modules are encapsulated
    % by another structure, there is no default 'summarize model' method
    % that could be defined and we have to assume many things about the
    % the model as we have done in the above.
    
end

% Sort all the plots according to their score_test_corr (if it exists)
score = sortrows(score);
newax = {};
sorted_filenames = {};
for ii = 1:size(score, 1)
    newax{ii} = ax{score(ii,2)};
    sorted_filenames{ii} = filenames{score(ii,2)};
end
ax = newax;

% Set up the scrollbar and its callback. 
function scroll_callback(h, evts, hds)
    
    y0 = - get(h,'Value');
    
    % Update all the positions of the graphs
    for ii = 1:length(filenames)
        set(ax{ii}.compr, 'Position', compression_pos(ii));
        set(ax{ii}.coefs, 'Position', fircoefs_pos(ii));
        set(ax{ii}.nlin,  'Position', nonlinearity_pos(ii));
        set(ax{ii}.text,  'Position', text_pos(ii));
    end
end

y0max = h*length(filenames) - h_win;
scr = uicontrol('Parent', fig, 'Style','slider', 'Units', 'pixels', ...
          'Position', [w 0 20 h_win], 'Min', 0-eps, 'Max', y0max, 'Value', 0, ...
          'Callback', @scroll_callback);
% 'SliderStep', [0.02 (h_win / y0max)], ...  

scroll_callback(scr, [], []);

% % ----------------------
% % % Also, show the test/train curves
% fh_tt = figure;
% plot(1:length(filenames), score(:,1), ...
%      1:length(filenames), score(:,3));
% xticks(1:length(filenames));
% 
% % label the x axis with filenames at a 90 degree angle
% lab = char(sorted_filenames);
% hx = get(gca,'XLabel');
% set(hx, 'Units', 'data');
% pos = get(hx, 'Position');
% 
% for ii = 1:size(lab,1)
%     t(ii) = text(ii, pos(2), lab(ii,:), 'Interpreter', 'none');
% end
% 
% P=get(gca,'Position');
% set(t, 'Rotation', 90, 'HorizontalAlignment', 'right') 
% set(gca, 'Position', [0.1300    0.4100    0.7750    0.5150])
% 
% legend('Test Score', 'Training Score', 'Location', 'NorthWest');
% 
% title(filenames{1}, 'Interpreter', 'none');

% fh1 = fig;
% fh2 = fh_tt;
end
