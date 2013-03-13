function handles = narf_browser(parent_handle)

global STACK XXX META NARF_PATH NARF_SAVED_MODELS_PATH;

if ~exist('parent_handle', 'var')
    parent_handle = figure('Menubar','figure', 'Resize','off', ...
       'Units','pixels', 'Position', [20 50 1200 900]);
end

pos = get(parent_handle, 'Position');
w = pos(3); % Width of the parent panel
h = pos(4); % Height of the parent panel
mw = 250;   % Menu width
ih = 600;   % Image section height
lb = 10;    % Left border spacing
ts = 25;    % Text spacing
tw = 80;    % Text label width
ta = -3;    % Text alignment vertical spacing

% Struct array of DB hits
db_results = [];
sortby = 'R_test';
sortdir = 'DESC';

% Create a panel inside it which will slide
handles.container_panel = uipanel('Parent',parent_handle, ...
    'Units','pixels', 'Position', [mw 0 w-20 ih]);

% Create the scroll bar
handles.container_slider = uicontrol('Parent',parent_handle, ...
    'Style','slider', 'Enable','off', ...
    'Units','pixels', 'Position', [w-20 0 20 ih], ...
    'Min', 0-eps, 'Max', 0, 'Value', 0, ...
    'Callback', @(h,evts,hds) update_panel_positions());

% Create the data table
handles.db_results = uitable('Parent', parent_handle, ...
        'Enable', 'on',  'Units', 'pixels', 'RowName', [],...
        'ColumnWidth', {50, 40, 100, 300, 50, 50, 70, 150}, ...
        'ColumnName', {'ID', 'Batch', 'CellID', 'ModelName', 'R_test', 'R_fit', 'Sparsity', 'Last Mod.'}, ...
        'Position', [mw ih w-mw h-ih]);
    
hJScroll = findjobj(handles.db_results); 
hJTable = hJScroll.getViewport.getView; 
hJTable.setNonContiguousCellSelection(false);
hJTable.setColumnSelectionAllowed(false);
hJTable.setRowSelectionAllowed(true);
hJTablecb = handle(hJTable, 'CallbackProperties');
set(hJTablecb, 'MousePressedCallback', {@results_row_select, gcf});
set(hJTablecb, 'KeyPressedCallback', {@results_row_select, gcf});
panax = axes('Units','normal', 'Position', [0 0 1 1],...
    'Parent', handles.container_panel);

    function results_row_select(a,~,~)
        r = a.getSelectedRows();
        res = db_results(r+1);
        [I,map] = imread(char(res.figurefile),'png');        
        axes(panax);
        imshow(I, map, 'InitialMagnification', 25);
    end

    % Callback to allow scrolling of image summary
    function update_panel_positions()
        hSld = handles.container_slider;
        hPan = handles.container_panel;
        offset = get(hSld,'Value');
        p = get(hPan, 'Position');
        set(hPan, 'Position', [p(1) -offset p(3) p(4)]);
        top = 1000;
        for yi = 1:N
            p = get(STACK{yi}.gh.fn_panel, 'Position');
            set(STACK{yi}.gh.fn_panel, 'Position', [0 (top-ph*yi-offset) p(3) p(4)]);
        end
        drawnow;
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

condition_handles = [];

% Topmost label
uicontrol('Parent', parent_handle, 'Style', 'text', 'Units', 'pixels',...
          'String', 'DB Query Conditions', ...
          'Position', [lb h-ts+ta mw-lb ts]);

    function make_condition_widget(n, label)
        
        uicontrol('Parent', parent_handle, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', label, ...
          'Position', [lb h-(ts*(n+1))+ta tw-lb ts]);
      
        condition_handles(n) = uicontrol('Parent', parent_handle, ...
        'Style', 'popupmenu', 'Enable', 'on', ...
        'String', '*', ...
        'Units', 'pixels', 'Position', [tw h-(ts*(n+1)) mw-tw ts], ...
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
            sql = [sql whereand() 'modelname LIKE "%' sel_token1 '%"'];
        end    
                
        if ~isequal(sel_token2, '*')
            sql = [sql whereand() 'modelname LIKE "%' sel_token2 '%"'];
        end    
        
        if ~isequal(sel_token3, '*')
            sql = [sql whereand() 'modelname LIKE "%' sel_token3 '%"'];
        end    
        
        if ~isequal(sel_token4, '*')
            sql = [sql whereand() 'modelname LIKE "%' sel_token4 '%"'];
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
        
        % IFF CELLID & BATCH ARE DEFINED
        ret = mysql(['SELECT DISTINCT modelname FROM (' inner_sql ') AS sq']);
        ret = cellstr(char(ret(:).modelname));
        ret = sprintf('%s_', ret{:}); % Big long string
        toks = tokenize_modelname(ret); % Break it up into chunks
        toks = cat(2, {'*'}, unique([toks{:}])); % Create options
        
        set(condition_handles(3), 'String', toks);
        set(condition_handles(4), 'String', toks);
        set(condition_handles(5), 'String', toks);
        set(condition_handles(6), 'String', toks);
        
        set(condition_handles(8), 'String', {'batch', 'cellid', 'id', 'lastmod', 'modelname' 'r_fit', 'r_test', 'sparsity'});
        sortby = popup2str(condition_handles(8));
        set(condition_handles(9), 'String', {'DESC', 'ASC'});
        sortdir = popup2str(condition_handles(9));
        
        
        % Update the popup selected values
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
        %fprintf('Updating Query Results Table\n');
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
        db_results = mysql([sql_query_builder() ' ORDER BY ' sortby ' ' sortdir ' LIMIT 0, 1000']);
        % TODO: Perhaps adding a 'sort' field would be useful?
        update_query_results_table();        
    end

% Call the callback once to get things running
any_condition_changed_callback();

end