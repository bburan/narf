function plotpath = plot_model_summary()
% plotpath = plot_model_summary()
% 
% Plots the model currently loaded in NARF and saves a PNG image at the
% appropriate path. Does NOT do anything else, so make sure your model is
% ready to be plotted.
%
% ARGUMENTS: none
%  
% RETURNS:
%    plotpath    The path to the PNG file that was just generated.
%
% Example:
%   load_model('/path/to/mymodel.mat');
%   recalc_xxx(1);
%   plot_model_summary();

global META XXX STACK NARF_SAVED_IMAGES_PATH;

lb = 30;  % Left border
bb = 30;  % Bottom border
th = 180; % text height
ph = 200; % Plot height

w = 800; % Pixels

vspace = 0.2; % Relative sizes
hspace = 0.05; 

% Scan through the STACK looking for things to .auto_plot
ap = [];
for ii = 1:length(STACK)
    if iscell(STACK{ii})
        m = STACK{ii}{1};
    else
        m = STACK{ii};
    end
    
    if isfield(Sm, 'auto_plot')
        ap(end+1) = ii;
    end
end       
nplots = length(ap);

% Create a new, invisible figure
h = (nplots+1)*ph + th;
fig = figure('Menubar', 'figure', 'Resize','off', 'Visible', 'on', ...
             'Units','pixels', 'Position', [50 50 w+20 h]);

% Call the auto-plot functions
for ii = 1:nplots
    idx = ap(ii);  
    mm = STACK{idx};
    if iscell(mm)    
        lm = length(mm);
        for jj = 1:lm,
            m = mm{jj};
        
            plotfn = m.auto_plot;
    
            aw = (w-lb*2) / lm;
            ax = axes('Parent', fig, 'Units', 'pixels', ...
                 'Position', [lb+(aw*(jj-1)) (nplots-ii)*ph+bb aw-lb ph-bb]);
    
            plotfn(STACK(1:idx), XXX(1:idx+1));          
        end
    else
        m = mm;
        
        plotfn = m.auto_plot;
    
        ax = axes('Parent', fig, 'Units', 'pixels', ...
            'Position', [lb (nplots-ii)*ph+bb w-lb*2 ph-bb]);
    
        plotfn(STACK(1:idx), XXX(1:idx+1));
    end
end

% TEXT AT TOP
if ~isfield(META, 'batch')
    META.batch = 0;
end

% Print the text at the top
axtt = axes('Parent', fig , 'Units','pixels', ...
    'Position', [lb nplots*ph+bb w-lb*2 ph-bb]);
set(gca,'xtick',[]); set(gca,'xticklabel',[]);
set(gca,'ytick',[]); set(gca,'yticklabel',[]);
ax_text  = text('Interpreter', 'none', 'Position', [0.05, 0.45], ...
   'String', sprintf('Batch:     %d\nCellid:    %s\nModel:     %s\nTrain Set: %s\nTest Set:  %s\nTrain r:   %.5f\nTest r:    %.5f', ...
                       META.batch, XXX{1}.cellid, META.modelname, ...
                       [XXX{1}.training_set{:}], [XXX{1}.test_set{:}], ...
                       XXX{end}.score_train_corr, ...
                       XXX{end}.score_test_corr));
                   
% Save the file
celldir = [NARF_SAVED_IMAGES_PATH filesep num2str(META.batch) filesep XXX{1}.cellid];
if ~exist(celldir)
    mkdir(celldir);
end
pngfile = [celldir filesep META.modelfile '.png'];
set(gcf,'PaperPositionMode','auto');
set(gcf,'InvertHardcopy','off');
set(0,'defaultTextFontName', 'Arial');
print(fig, pngfile, '-dpng');

% SVD- Don't close figure so that it gets saved by queuerun.m as well
%close(fig);

plotpath = pngfile;