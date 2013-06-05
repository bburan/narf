function handle = narf_analysis(parent_handle)
% A GUI to edit the contents of the MySQL NarfAnalysis table

narf_set_path;
dbopen;

h = 800;  % Total window height
w = 1200; % Total window width
dh = 400; % Vertical division height
mh = 200; % Middle division height
lw = 300; % Left/Right division point
cw = 400; % Right panel width
rw = w-lw-cw;
ts = 25;  % Text spacing
tw = 70;  % Text label width
pad = 3;

if ~exist('parent_handle', 'var')
    parent_handle = figure('Menubar','none', 'Resize','off', ...
       'Units','pixels', 'Position', [20 50 w h]);
end

db_results = []; % Shared amongst local functions that need common SQL query results

% Three things are important enough to keep the handle around
left_panel = uipanel('Parent', parent_handle, ...
    'Units','pixels', 'Position', [0 h-dh lw dh]);

center_panel = uipanel('Parent', parent_handle, ...
    'Units','pixels', 'Position', [lw h-dh cw dh]);

right_panel = uipanel('Parent', parent_handle, ...
    'Units','pixels', 'Position', [lw+cw h-dh rw dh]);

db_results_table = uitable('Parent', parent_handle, ...
        'Enable', 'on',  'Units', 'pixels', 'RowName', [],...
        'ColumnWidth', {50, 40, 60, 50, 70, 200, 100, 150, 70, 70, 70, 70}, ...
        'ColumnName', {'ID', 'Batch', 'CellID', 'EstSet', 'ValSet', 'Model Name', 'Val R', 'Val NLL', 'Est R', 'Est NLL', 'Last Mod.', 'Notes'}, ...
        'Position', [0 0 w h-dh-ts]);

% Left Panel
uicontrol('Parent', left_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'center', 'String', 'ANALYSIS SELECTOR', ...
          'Position', [pad dh-(ts*1)-pad lw-pad ts]);
      
uicontrol('Parent', left_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Status Filter:', ...
          'Position', [pad dh-(ts*2)-pad 100 ts]);
      
handle.status_filter = uicontrol('Parent', left_panel, 'Style', 'popupmenu', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', '*', ...
          'Position', [100 dh-(ts*2) lw-100-pad ts]);    

uicontrol('Parent', left_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Tag Filter:', ...
          'Position', [pad dh-(ts*3)-pad 100 ts]);
      
handle.tag_filter = uicontrol('Parent', left_panel, 'Style', 'popupmenu', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', '*', ...
          'Position', [100 dh-(ts*3) lw-100-pad ts], ...
          'Callback', @any_condition_changed_callback);    

handle.analyses_table = uitable('Parent', left_panel, 'Enable', 'on', 'Units', 'pixels',...
        'RowName', [], ...
        'ColumnWidth', {30, 60, lw-pad*2-120}, ...
        'ColumnName', {'ID', 'Status', 'Analysis Name'}, ...
        'Position', [pad pad+ts lw-pad*2 dh-4*ts-pad]);
    
handle.delete_analysis = uicontrol('Parent', left_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Delete', ...
          'Position', [pad pad (1/3)*(lw-pad*2) ts], ...
          'Callback', @any_condition_changed_callback);

handle.refresh_analysis = uicontrol('Parent', left_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Refresh', ...
          'Position', [pad+(2/3)*(lw-pad*2) pad (1/3)*(lw-pad*2) ts], ...
          'Callback', @any_condition_changed_callback); 
      
% Center panel      
uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Name:', ...
          'Position', [pad dh-(ts*1)-pad 70 ts-pad]);
      
handle.name = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Sample Analysis', ...
          'Position', [70+pad dh-(ts*1) 200 ts-pad], ...
          'Callback', @any_condition_changed_callback);
      
uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Status:', ...
          'Position', [70+200+pad dh-(ts*1)-pad 70 ts-pad]);
      
handle.status = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'active', ...
          'Position', [70+200+70+pad dh-(ts*1) cw-(70+200+70+2*pad) ts-pad], ...
          'Callback', @any_condition_changed_callback);

uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Q:', ...
          'Position', [pad dh-(ts*2)-pad 20 ts-pad]);
      
handle.short_desc = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'Max', 3, 'Min', 1, ...
          'HorizontalAlignment', 'left', 'String', 'Motivating question or experimental summary description goes here.', ...
          'Position', [20+pad dh-(ts*4) cw-20-pad ts*3-pad]);    

uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...         
          'HorizontalAlignment', 'left', 'String', 'A:', ...
          'Position', [pad dh-(ts*6) 20 ts*2]);
      
handle.long_desc = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'Max', 10, 'Min', 1, ...
          'HorizontalAlignment', 'left', 'String', 'Answer or long description goes here.', ...
          'Position', [20+pad ts cw-20-pad dh-ts*5]);   
    
uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...          
          'HorizontalAlignment', 'left', 'String', 'Tags:', ...
          'Position', [pad 0 50 ts-pad]);

handle.tags = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'unsorted, weird, goodresult, etc', ...
          'Position', [50+pad pad cw-50-pad*2 ts-pad]);   

% Right Panel is for the model structures

uicontrol('Parent', right_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'ModelTree:', ...
          'Position', [pad dh-(ts*1)-pad 70 ts-pad]);
      
handle.modeltree = uicontrol('Parent', right_panel, 'Style', 'edit', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', '{''env100'', ''log2b'', ''firn'', {''npnl'', ''npfnl''}, ''boost''}', ...
          'Position', [pad dh-(ts*2) rw-pad*2 ts-pad], ...
          'Callback', @any_condition_changed_callback);

handle.regene_modeltree = uicontrol('Parent', right_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Regenerate List', ...
          'Position', [rw-150-pad dh-(ts)-pad 150 ts], ...
          'Callback', @any_condition_changed_callback); 
      
handle.modellist = uicontrol('Parent', right_panel, 'Style', 'listbox', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', '', ...
          'Position', [pad ts rw-pad*2 dh-3*ts-pad], ...
          'Callback', @any_condition_changed_callback);

% Dividing button row
uicontrol('Parent', parent_handle, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Add Above Selected to Batch', ...
          'Position', [400 h-dh-ts 200 ts-pad], ...
          'Callback', @any_condition_changed_callback);

uicontrol('Parent', parent_handle, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Delete Below Selected from Batch', ...
          'Position', [700 h-dh-ts 250 ts-pad], ...
          'Callback', @any_condition_changed_callback);   

uicontrol('Parent', parent_handle, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Save and Update Batch', ...
          'Position', [w-200 h-dh-ts 200 ts-pad], ...
          'Callback', @any_condition_changed_callback);
      
% Set up the DB Results table widget behavior
hJScroll = findjobj(db_results_table); 
hJTable = hJScroll.getViewport.getView; 
hJTable.setNonContiguousCellSelection(false);
hJTable.setColumnSelectionAllowed(false);
hJTable.setRowSelectionAllowed(true);
% hJTablecb = handle(hJTable, 'CallbackProperties');
% set(hJTablecb, 'MousePressedCallback', {@get_selected_row, gcf});
% set(hJTablecb, 'KeyPressedCallback', {@get_selected_row, gcf});
% 
% % TODO: When the delete button is pressed, delete selected rows
% 
%     function get_selected_row(a,~,~)
%         r = a.getSelectedRows();
%         res = db_results(r+1);
%     end

    function sql = sql_query_builder()       
        sql = 'SELECT * FROM NarfAnalysis';
        
        sel_batch  = popup2str(condition_handle(1));
        sel_cellid = popup2str(condition_handle(2));
        sel_token1 = popup2str(condition_handle(3));
        sel_token2 = popup2str(condition_handle(4));
        sel_token3 = popup2str(condition_handle(5));
        sel_token4 = popup2str(condition_handle(6));
        
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
    end
% 
%     function rebuild_animal_options()
%         sql = ['SELECT DISTINCT animal FROM (SELECT * FROM gAnimal WHERE lab LIKE "lbhb") AS sq'];
%         ret = mysql(sql);
%         animals = cellstr(char(ret(:).animal));
%         set(handle.animals, 'String', cat(2, {'*'}, animals));
%     end
% 
%     function any_condition_changed_callback()        
%         idx = get(handle.animals, 'Value');
%         list = get(hObject,'String');
%         animal = list{idx};
%         
%         % Convert that selection into a short prefix        
%         sql = ['SELECT DISTINCT cellprefix FROM (SELECT * FROM gAnimal WHERE animal LIKE "' animal '") AS sq'];
%         ret = mysql(sql);
%         if (length(ret) > 1)
%             error('How did I get more than one prefix for an animal!?');
%         end
%         prefix = char(ret(1).cellprefix);
%         
%         sql =  ['SELECT * from sRunData WHERE cellid REGEXP "' prefix '*"'];
%         ret = mysql(sql);
%         cellids = cellstr(char(ret(:).cellid));
%         set(handle.cellids, 'String', cellids);
%     end

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
        set(db_results_table, 'Data', c);
        drawnow;
    end

    function any_condition_changed_callback(~,~,~)
%         rebuild_condition_options();
%         dbopen; 
%         db_results = mysql([sql_query_builder() ' ORDER BY ' sortby ' ' sortdir ' LIMIT 0, 1000']);
%         update_query_results_table();        
%         selected = [];
%         axes(panax); cla;
        fprintf('Callback called\n');
    end

%rebuild_animal_options();
% Call the callback once to get things running
% any_condition_changed_callback();

end



