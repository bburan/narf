function plotpath = plot_model_summary()
% Plots the model currently loaded into global memory. 
% Saves the plot to the appropriate path.
% Returns the path of the saved plot.
% Does NOT do anything else, so make sure your model is completely ready.
% Usually that means load_model(), recalc_xxx(1), and then this function.

global META XXX STACK NARF_SAVED_IMAGES_PATH;

lb = 30;  % Left border
bb = 30;  % Bottom border
th = 180;  % text height
ph = 200; % Plot height

w = 800; % Pixels
h = 3*ph + th;

vspace = 0.2; % Relative sizes
hspace = 0.05; 

% Create a new, invisible figure
fig = figure('Menubar', 'figure', 'Resize','off', 'Visible', 'off', ...
             'Units','pixels', 'Position', [50 50 w+20 h]);

% FIR FILTER PLOT
ax_fir = axes('Parent', fig, 'Units','pixels', 'Position', [lb 2*ph+bb w-lb*2 ph-bb]);
[firmod, firmod_idx] = find_module(STACK, 'fir_filter');
sparsity = NaN;
if ~isempty(firmod) 
    sparsity = sparsity_metric(firmod.coefs);
    firmod.plot_fns{1}.fn(STACK(1:firmod_idx), XXX(1:firmod_idx+1));
end

% NONLINEARITY PLOT
ax_nl = axes('Parent', fig, 'Units','pixels', 'Position', [lb 1*ph+bb w-lb*2 ph-bb]);
nlidx = firmod_idx + 2; % TODO: This assumes firn or depn, which is not always true!
if strcmp(STACK{nlidx}.name, 'nonlinearity') | ...
        strcmp(STACK{nlidx}.name, 'sparse_empirical_nonlinearity') | ...
        strcmp(STACK{nlidx}.name, 'gmm_nonlinearity') | ...
        strcmp(STACK{nlidx}.name, 'nonparm_filter_nonlinearity') | ...
        strcmp(STACK{nlidx}.name, 'nonparm_nonlinearity')
    STACK{nlidx}.plot_fns{1}.fn(STACK(1:nlidx), XXX(1:nlidx+1));
end

% PREDICTION/REALITY
ax_pred = axes('Parent', fig, 'Units','pixels', 'Position', [lb bb w-lb*2 ph-bb]);
[corrmod, corrmod_idx] = find_module(STACK, 'correlation');
if ~isempty(corrmod) 
    corrmod.plot_fns{1}.fn(STACK(1:corrmod_idx), XXX(1:corrmod_idx+1));
end

% TEXT AT TOP
if ~isfield(META, 'batch')
    META.batch = 240;
end

% This is the right way to do it, but MATLAB is buggy with PNG exports when
% the figure is hidden. I guess a callback isn't getting called.
% ax_text  = uicontrol('Parent', fig, 'Style', 'text', 'HandleVisibility', 'on', ...
%     'Units', 'pixels', 'FontName', 'FixedWidth', ...
%     'HorizontalAlignment', 'left', ...
%     'Position', [lb (h-th) w th], ...
%     'String', sprintf('Batch: %d\nCellid: %s\nModel: %s\nTrain Set: %s\nTest Set: %s\nTrain r: %.5f\nTest r: %.5f\nSparsity: %.5f', ...
%                        META.batch, XXX{1}.cellid, META.modelname, ...
%                        [XXX{1}.training_set{:}], [XXX{1}.test_set{:}], ...
%                        XXX{end}.score_train_corr, ...
%                        XXX{end}.score_test_corr, sparsity));

% This is the wrong way to do it, but works kind of
axtt = axes('Parent', fig , 'Units','pixels', 'Position', [lb 3*ph+bb w-lb*2 ph-bb]);
set(gca,'xtick',[]); set(gca,'xticklabel',[]);
set(gca,'ytick',[]); set(gca,'yticklabel',[]);
ax_text  = text('Interpreter', 'none', 'Position', [0.05, 0.45], ...
   'String', sprintf('Batch:     %d\nCellid:    %s\nModel:     %s\nTrain Set: %s\nTest Set:  %s\nTrain r:   %.5f\nTest r:    %.5f\nSparsity:  %.5f', ...
                       META.batch, XXX{1}.cellid, META.modelname, ...
                       [XXX{1}.training_set{:}], [XXX{1}.test_set{:}], ...
                       XXX{end}.score_train_corr, ...
                       XXX{end}.score_test_corr, sparsity));
                   
% Save the file
celldir = [NARF_SAVED_IMAGES_PATH filesep XXX{1}.cellid];
if ~exist(celldir)
    mkdir(celldir);
end
pngfile = [celldir filesep META.modelfile '.png'];
set(gcf,'PaperPositionMode','auto');
set(gcf,'InvertHardcopy','off');
set(0,'defaultTextFontName', 'Arial');
print(fig, pngfile, '-dpng');
close(fig);

plotpath = pngfile;