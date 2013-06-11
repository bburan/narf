function handles = narf_browser(parent_handle)

global STACK XXX META MODULES;

narf_set_path
dbopen;
MODULES = scan_directory_for_modules();   

h = 1000;  % Total window height
w = 1900; % Total window width
rh = 550; % Results height
rw = 800; % Results width
mw = 190; % Menu width
iw = 450; % Image width
ih = 650; % Image height
ts = 25;  % Text spacing
tw = 70;  % Text label width
ta = -3;  % Text alignment vertical spacing

if ~exist('parent_handle', 'var')
    parent_handle = figure('Menubar','none', 'Resize','off', ...
       'Units','pixels', 'Position', [20 50 w h],...
       'Name', 'NARF Model Browser', 'NumberTitle', 'off');
end

db_results = [];
sortby = 'r_test';
sortdir = 'DESC';
selected = [];

handles.image_panel = uipanel('Parent',parent_handle, ...
    'Units','pixels', 'Position', [mw rh-40 rw-mw ih]);

handles.viewer_panel = uipanel('Parent',parent_handle, ...
    'Units','pixels', 'Position', [rw 0 w-(rw) h]);

handles.db_results = uitable('Parent', parent_handle, ...
        'Enable', 'on',  'Units', 'pixels', 'RowName', [],...
        'ColumnWidth', {50, 40, 100, 280, 50, 50, 70, 150}, ...
        'ColumnName', {'ID', 'Batch', 'CellID', 'ModelName', 'R_test', 'R_fit', 'Sparsity', 'Last Mod.'}, ...
        'Position', [0 0 rw rh]);
    
hJScroll = findjobj(handles.db_results); 
hJTable = hJScroll.getViewport.getView; 
hJTable.setNonContiguousCellSelection(false);
hJTable.setColumnSelectionAllowed(false);
hJTable.setRowSelectionAllowed(true);
hJTablecb = handle(hJTable, 'CallbackProperties');
set(hJTablecb, 'MousePressedCallback', {@results_row_select, gcf});
set(hJTablecb, 'KeyPressedCallback', {@results_row_select, gcf});
panax = axes('Units','normal', 'Position', [0 0 1 1],...
             'Parent', handles.image_panel);

    function results_row_select(a,~,~)
        r = a.getSelectedRows();
        res = db_results(r+1);
        selected = res; 
        [I,map] = imread(char(res.figurefile),'png');  
        I = imresize(I, (ih/size(I,2)));      
        axes(panax);
        imshow(I, map);
    end

global STACK XXX META;
delete_all_module_guis();
STACK = {};
XXX = {};
META = {};
modelpane_updater = narf_modelpane(handles.viewer_panel); 

CustomFunctionHandle=uicontrol('Parent', parent_handle, 'Style', 'edit',...
    'String', 'dummy',...
    'Units','pixels', 'Position', [round(mw./2)+1 rh round(mw./2) 25], ...
    'Callback', @custom_model_analysis);

uicontrol('Parent', parent_handle, 'Style', 'pushbutton',...
    'String', 'Custom:', ...
    'Units','pixels', 'Position', [0 rh round(mw./2) 25], ...
    'Callback', @custom_model_analysis);

function custom_model_analysis(~,~,~)
%global STACK XXX META;
    cust_function_name=get(CustomFunctionHandle,'String');
    if ~isempty(selected)
        delete_all_module_guis();
        load_model(char(selected.modelpath));  
        
        feval(cust_function_name,STACK,XXX,META);
    end
end

uicontrol('Parent', parent_handle, 'Style', 'pushbutton',...
    'String', 'Inspect Selected Model', ...
    'Units','pixels', 'Position', [0 rh+25 mw 25], ...
    'Callback', @view_selected_in_modelpane);

    function view_selected_in_modelpane(~,~,~)
        if ~isempty(selected)
            delete_all_module_guis();
            load_model(char(selected.modelpath));  
            modelpane_updater();
        end
    end

uicontrol('Parent', parent_handle, 'Style', 'pushbutton',...
    'String', 'Preview in new Window', ...
    'Units','pixels', 'Position', [0 rh+50 mw 25], ...
    'Callback', @view_image_in_new_window);

    function view_image_in_new_window(~,~,~)
        if ~isempty(selected)
            figure('Menubar', 'none');
            [I,map] = imread(char(selected.figurefile),'png');  
            I = imresize(I, 0.666); % ...the number of the beast! Rowr!
            imshow(I, map);
        end
    end

uicontrol('Parent', parent_handle, 'Style', 'pushbutton',...
    'String', 'View Selected Cell STRF', ...
    'Units','pixels', 'Position', [0 rh+75 mw 25], ...
    'Callback', @view_strf_button_callback);

    function view_strf_button_callback(~, ~, ~)
        if isempty(selected)
            return
        end
        [cfd, ~, ~] = dbgetscellfile('cellid', selected.cellid);
        for i = 1:length(cfd);
            % TODO: Replace magic number  1 with better description of TOR files
            if (cfd(i).runclassid == 1)                
                figure;
                strf_offline2([cfd(i).stimpath cfd(i).stimfile], ...
                    [cfd(i).path cfd(i).respfile], ...
                    cfd(i).channum, cfd(i).unit);
            end
        end
    end


uicontrol('Parent', parent_handle, 'Style', 'pushbutton',...
    'String', 'View Selected in CellDB', ...
    'Units','pixels', 'Position', [0 rh+100 mw 25], ...
    'Callback', @launch_celldb_in_browser);
    
    function launch_celldb_in_browser(~,~,~)
        if ~isempty(selected.cellid)
            % If you want to try launching firefox directly, you will need
            % to set the environment variables and create symlinks to
            % libstdc++, and it's a pain.
            %setenv('LD_LIBRARY_PATH', [getenv('LD_LIBRARY_PATH'), ':/usr/lib/firefox'])
            %setenv('LD_LIBRARY_PATH', [getenv('LD_LIBRARY_PATH'), ':/usr/lib/x86_64-linux-gnu'])
            % etc...
            
            web(['http://hyrax.ohsu.edu/celldb/peninfo.php?penname=' ...
                 selected.cellid(1:6) '#' selected.cellid(1:7)]);
        end
    end

uicontrol('Parent', parent_handle, 'Style', 'pushbutton',...
    'String', 'Plot Cumulative Curves', ...
    'Units','pixels', 'Position', [0 rh+125 mw 25], ...
    'Callback', @plot_cumulative_curves);
    
    function plot_cumulative_curves(~,~,~)
        
        sel_batch  = popup2str(condition_handles(1));
        if strcmp(sel_batch, '*'),
            return
        end
        
        sel_cellid = popup2str(condition_handles(2));
        if strcmp(sel_cellid, '*'),
            sel_cellid = '';
        end
        
        inner_sql = sql_query_builder(); 
        ret = mysql(['SELECT DISTINCT modelname FROM (' inner_sql ') AS sq']);
        ret = cellstr(char(ret(:).modelname));
        ret2 = sprintf('%s_', ret{:}); % Big long string
        toks = tokenize_string(ret2); % Break it up into tokens
        toks = unique([toks{:}]); % Only unique tokens
            
        % Display only tokens not common across all modelfiles
        disptoks = {};
        for tt = 1:length(toks)
            t = toks{tt};
            
            tmp={};
            for rr = 1:length(ret)
                garr = tokenize_string(ret{rr}); 
                tmp{rr} = [garr{:}];                
            end
            
            if ~all(cellfun(@(x) any(strcmp(x, t)), tmp))
                disptoks{end+1} = t;
            end
        end
                
        holdtoks = {popup2str(condition_handles(3)), ...
                    popup2str(condition_handles(4)), ...
                    popup2str(condition_handles(5)), ...
                    popup2str(condition_handles(6))};
                
        star_cells = cellfun(@(x)strcmp('*', x) , holdtoks);
        holdtoks(star_cells) = []; % Remove *'d cells
        
        plot_cumulative_performance(sel_batch, ...
            sel_cellid, ...
            holdtoks, ...
            disptoks, ...
            popup2str(condition_handles(8)));
        
    end
    
% SVD: scatter plot button
uicontrol('Parent', parent_handle, 'Style', 'pushbutton',...
    'String', 'Scatter Plot', ...
    'Units','pixels', 'Position', [0 rh+150 mw./2 25], ...
    'Callback', @plot_scatter);
    
    function plot_scatter(~,~,~)
        
        sel_batch  = popup2str(condition_handles(1));
        if strcmp(sel_batch, '*'),
            return
        end
        
        sel_cellid = popup2str(condition_handles(2));
        if strcmp(sel_cellid, '*'),
            sel_cellid = '';
        end
        
        inner_sql = sql_query_builder(); 
        ret = mysql(['SELECT DISTINCT modelname FROM (' inner_sql ') AS sq']);
        ret = cellstr(char(ret(:).modelname));
        ret2 = sprintf('%s_', ret{:}); % Big long string
        toks = tokenize_string(ret2); % Break it up into tokens
        toks = unique([toks{:}]); % Only unique tokens
        
        % Display only tokens not common across all modelfiles
        disptoks = {};
        for tt = 1:length(toks)
            t = toks{tt};
            
            if any(cellfun(@(x) isempty(strfind(x, t)), ret))
                % If any are empty, add it
                disptoks{end+1} = t;
            end
        end
        
        holdtoks = {popup2str(condition_handles(3)), ...
                    popup2str(condition_handles(4)), ...
                    popup2str(condition_handles(5)), ...
                    popup2str(condition_handles(6))};
                
        star_cells = cellfun(@(x)strcmp('*', x) , holdtoks);
        holdtoks(star_cells) = []; % Remove *'d cells
        plot_performance_scatter(sel_batch, ...
            sel_cellid, ...
            holdtoks,...
            char(selected.modelname), ...
            popup2str(condition_handles(8)));
        
    end
    
% SVD: bar plot summary button
uicontrol('Parent', parent_handle, 'Style', 'pushbutton',...
    'String', 'Bar Plot', ...
    'Units','pixels', 'Position', [mw./2 rh+150 mw./2 25], ...
    'Callback', @plot_bar);
    
    function plot_bar(~,~,~)
        
        sel_batch  = popup2str(condition_handles(1));
        if strcmp(sel_batch, '*'),
            return
        end
        
        sel_cellid = popup2str(condition_handles(2));
        if strcmp(sel_cellid, '*'),
            sel_cellid = '';
        end
        
        inner_sql = sql_query_builder(); 
        ret = mysql(['SELECT DISTINCT modelname FROM (' inner_sql ') AS sq']);
        ret = cellstr(char(ret(:).modelname));
        ret2 = sprintf('%s_', ret{:}); % Big long string
        toks = tokenize_string(ret2); % Break it up into tokens
        toks = unique([toks{:}]); % Only unique tokens
        
        holdtoks = {popup2str(condition_handles(3)), ...
                    popup2str(condition_handles(4)), ...
                    popup2str(condition_handles(5)), ...
                    popup2str(condition_handles(6))};
                
        star_cells = cellfun(@(x)strcmp('*', x) , holdtoks);
        holdtoks(star_cells) = []; % Remove *'d cells
        
        plot_performance_bar(sel_batch, ...
            sel_cellid, ...
            holdtoks, ...
            toks, ...
            popup2str(condition_handles(8)));
        
    end
    


% TODO: Heat map button

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

condition_handles = [];

% Topmost label
uicontrol('Parent', parent_handle, 'Style', 'text', 'Units', 'pixels',...
          'String', 'DB Query Conditions', ...
          'Position', [0 h-ts+ta mw ts]);

    function make_condition_widget(n, label)
        
        uicontrol('Parent', parent_handle, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', label, ...
          'Position', [0 h-(ts*(n+1))+ta tw ts]);
      
        condition_handles(n) = uicontrol('Parent', parent_handle, ...
        'Style', 'popupmenu', 'Enable', 'on', ...
        'String', '*', ...
        'Units', 'pixels', 'Position', [tw h-(ts*(n+1)) (mw-tw) ts], ...
        'Callback', @any_condition_changed_callback);
    
    end

make_condition_widget(1, 'Batch #:');
make_condition_widget(2, 'CellID:'); 
make_condition_widget(3, 'Token 1:');
make_condition_widget(4, 'Token 2:');
make_condition_widget(5, 'Token 3:');
make_condition_widget(6, 'Token 4:');

make_condition_widget(8, 'Sort By:');
make_condition_widget(9, 'Direction:');

    function sql = sql_query_builder()       
        sql = 'SELECT * FROM NarfResults';
        
        sel_batch  = popup2str(condition_handles(1));
        sel_cellid = popup2str(condition_handles(2));
        sel_token1 = popup2str(condition_handles(3));
        sel_token2 = popup2str(condition_handles(4));
        sel_token3 = popup2str(condition_handles(5));
        sel_token4 = popup2str(condition_handles(6));
        
        isfirst = true;
        function s = whereand()
            if isfirst
                isfirst = false;
                s = ' WHERE ';
            else 
                s = ' AND ';
            end
        end
        
        if ~isequal(sel_batch, '*')
            sql = [sql whereand() 'batch=' sel_batch ''];
        end
        
        if ~isequal(sel_cellid, '*')
            sql = [sql whereand() 'cellid="', sel_cellid '"'];
        end
        
        if ~isequal(sel_token1, '*')
            sql = [sql whereand() 'modelname REGEXP "(^|[^[:alnum:]])' sel_token1 '([^[:alnum:]]|$)"'];
        end    
                
        if ~isequal(sel_token2, '*')
            sql = [sql whereand() 'modelname REGEXP "(^|[^[:alnum:]])' sel_token2 '([^[:alnum:]]|$)"'];
        end    
        
        if ~isequal(sel_token3, '*')
            sql = [sql whereand() 'modelname REGEXP "(^|[^[:alnum:]])' sel_token3 '([^[:alnum:]]|$)"'];
        end    
        
        if ~isequal(sel_token4, '*')
            sql = [sql whereand() 'modelname REGEXP "(^|[^[:alnum:]])' sel_token4 '([^[:alnum:]]|$)"'];
        end
    end

    function rebuild_condition_options()
        % Query DB with current settings
        inner_sql = sql_query_builder();        

        % Record current positions of selected things
        sel_batch  = popup2str(condition_handles(1));
        sel_cellid = popup2str(condition_handles(2));
        sel_token1 = popup2str(condition_handles(3));
        sel_token2 = popup2str(condition_handles(4));
        sel_token3 = popup2str(condition_handles(5));
        sel_token4 = popup2str(condition_handles(6));
        
        ret = mysql(['SELECT DISTINCT batch FROM (' inner_sql ') AS sq ORDER BY batch']);
        set(condition_handles(1), 'String', {'*', ret(:).batch});

        ret = mysql(['SELECT DISTINCT cellid FROM (' inner_sql ') AS sq ORDER BY cellid']);
        set(condition_handles(2), 'String', {'*', ret(:).cellid});
                
        ret = mysql(['SELECT DISTINCT modelname FROM (' inner_sql ') AS sq']);
        ret = cellstr(char(ret(:).modelname));
        ret = sprintf('%s_', ret{:}); % Big long string
        toks = tokenize_string(ret); % Break it up into chunks
        toks = cat(2, {'*'}, unique([toks{:}])); % Create unique options
        
        set(condition_handles(3), 'String', toks);
        set(condition_handles(4), 'String', toks);
        set(condition_handles(5), 'String', toks);
        set(condition_handles(6), 'String', toks);
        
        set(condition_handles(8), 'String', {'batch', 'cellid', 'id', 'lastmod', 'modelname' 'r_fit', 'r_test', 'sparsity'});
        sortby = popup2str(condition_handles(8));
        set(condition_handles(9), 'String', {'DESC', 'ASC'});
        sortdir = popup2str(condition_handles(9));
               
        % Update the popups so that selected values don't change
        set(condition_handles(1), 'Value', ...
            find(ismember(cellstr(get(condition_handles(1), 'String')), sel_batch)));
        set(condition_handles(2), 'Value', ...
            find(ismember(cellstr(get(condition_handles(2), 'String')), sel_cellid)));
        set(condition_handles(3), 'Value', ...
            find(ismember(cellstr(get(condition_handles(3), 'String')), sel_token1)));
        set(condition_handles(4), 'Value', ...
            find(ismember(cellstr(get(condition_handles(4), 'String')), sel_token2)));
        set(condition_handles(5), 'Value', ...
            find(ismember(cellstr(get(condition_handles(5), 'String')), sel_token3)));
        set(condition_handles(6), 'Value', ...
            find(ismember(cellstr(get(condition_handles(6), 'String')), sel_token4)));
    end

    function update_query_results_table()
        l = length(db_results);
        c = cell(l,8);
        for i = 1:l
            c{i,1} = db_results(i).id;
            c{i,2} = db_results(i).batch;
            c{i,3} = db_results(i).cellid;
            c{i,4} = char(db_results(i).modelname);
            c{i,5} = db_results(i).r_test;
            c{i,6} = db_results(i).r_fit;
            c{i,7} = db_results(i).sparsity;
            c{i,8} = db_results(i).lastmod;
        end
        set(handles.db_results, 'Data', c);
        drawnow;
    end

    function any_condition_changed_callback(~,~,~)
        rebuild_condition_options();
        dbopen; 
        db_results = mysql([sql_query_builder() ' ORDER BY ' sortby ' ' sortdir ' LIMIT 0, 1000']);
        update_query_results_table();        
        selected = [];
        axes(panax); cla;
    end

% Call the callback once to get things running
any_condition_changed_callback();

end