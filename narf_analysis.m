function parent_handle = narf_analysis(parent_handle)
% A GUI to edit the contents of the MySQL NarfAnalysis table

narf_set_path;
dbopen;

h = 800;  % Total window height
w = 1200; % Total window width
dh = 400; % Vertical division height
bh = h-dh; % Bottom division height
lw = 300; % Left/Right division point
cw = 400; % Right panel width
rw = w-lw-cw;
ts = 25;  % Text spacing
tw = 70;  % Text label width
pad = 3;

if ~exist('parent_handle', 'var')
    parent_handle = figure('Menubar','none', 'Resize','off', ...
       'Units','pixels', 'Position', [20 50 w h],...
       'Name', 'NARF Analysis Browser', 'NumberTitle', 'off');
end

db_results = []; % Shared amongst local functions that need common SQL query results
analyses_found = [];
sel_analysis = [];

left_panel = uipanel('Parent', parent_handle, ...
    'Units','pixels', 'Position', [0 h-dh lw dh]);

center_panel = uipanel('Parent', parent_handle, ...
    'Units','pixels', 'Position', [lw h-dh cw dh]);

right_panel = uipanel('Parent', parent_handle, ...
    'Units','pixels', 'Position', [lw+cw h-dh rw dh]);

bottom_panel = uipanel('Parent', parent_handle, ...
    'Units','pixels', 'Position', [0 0 w bh]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Left Panel
uicontrol('Parent', left_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'center', 'String', 'ANALYSIS SELECTOR', ...
          'Position', [pad dh-(ts*1)-pad lw-pad ts]);
      
uicontrol('Parent', left_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Status Filter:', ...
          'Position', [pad dh-(ts*2)-pad 100 ts]);
      
handles.status_filter = uicontrol('Parent', left_panel, 'Style', 'popupmenu', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', '*', ...
          'Position', [100 dh-(ts*2) lw-100-pad ts], ...
          'Callback', @any_condition_changed_callback);    

uicontrol('Parent', left_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Tag Filter:', ...
          'Position', [pad dh-(ts*3)-pad 100 ts]);
      
handles.tag_filter = uicontrol('Parent', left_panel, 'Style', 'popupmenu', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', '*', ...
          'Position', [100 dh-(ts*3) lw-100-pad ts], ...
          'Callback', @any_condition_changed_callback);    

handles.analyses_table = uitable('Parent', left_panel, ...
        'Enable', 'on', 'Units', 'pixels',...
        'RowName', [], ...
        'ColumnWidth', {30, 80, lw-pad*2-135}, ...
        'ColumnName', {'ID', 'Status', 'Analysis Name'}, ...
        'Position', [pad pad+ts lw-pad*2 dh-4*ts-pad]);

handles.new_analysis = uicontrol('Parent', left_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'New', ...
          'Position', [pad pad (1/3)*(lw-pad*2) ts], ...
          'Callback', @add_new_analysis_callback);

     function add_new_analysis_callback(~,~,~)
         sqlinsert('NarfAnalysis', ...
          'name',      'Untitled Analysis', ...
          'status',    'uninitialized', ...
	      'question',  'Please write a brief description or question here.', ...
	      'answer',    'Please write a longer discussion or answer here. ', ...
          'summaryfig', '');
         any_condition_changed_callback();
     end
  
 
handles.delete_analysis = uicontrol('Parent', left_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Delete', ...
          'Position', [pad+(1/3)*(lw-pad*2) pad (1/3)*(lw-pad*2) ts], ...
          'Callback', @delete_analysis_callback);
      
     function delete_analysis_callback(~,~,~)
         sel = questdlg(['Are you sure you want to delete ID=' num2str(sel_analysis.id) ...
             '? This cannot be undone.'], 'PAY ATTENTION!', ...
             'Delete', 'Cancel', 'Cancel'); 
         
         if strcmp(sel, 'Delete')
            sql = ['DELETE FROM NarfAnalysis WHERE id=' num2str(sel_analysis.id)];
            mysql(sql);
            any_condition_changed_callback();
         end         
     end        

handles.refresh_analysis = uicontrol('Parent', left_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Refresh', ...
          'Position', [pad+(2/3)*(lw-pad*2) pad (1/3)*(lw-pad*2) ts], ...
          'Callback', @any_condition_changed_callback); 

% Configure the analyses table selection to update the center and right panels
hJS = findjobj(handles.analyses_table); 
hJT = hJS.getViewport.getView;
hJT.setNonContiguousCellSelection(false);
hJT.setColumnSelectionAllowed(false);
hJT.setRowSelectionAllowed(true);
hJT.setSelectionMode(0); % Allow only a single row to be selected at once
hJTcb = handle(hJT, 'CallbackProperties');
set(hJTcb, 'MousePressedCallback', {@analyses_table_row_selected, gcf});
set(hJTcb, 'KeyPressedCallback', {@analyses_table_row_selected, gcf});

    function analyses_table_row_selected(a,~,~)
        r = a.getSelectedRows();
        sel_analysis = analyses_found(r+1);
        
        % Update
        set(handles.name, 'String', sel_analysis.name);
        set(handles.status, 'String', sel_analysis.status);
        set(handles.short_desc, 'String', char(sel_analysis.question));
        set(handles.long_desc, 'String', char(sel_analysis.answer));
        if isempty(sel_analysis.batch)
            set(handles.batch, 'String', 'SELECT A BATCH');
        else
            set(handles.batch, 'String', sel_analysis.batch);           
        end
        set(handles.tags, 'String', char(sel_analysis.tags));
        if isempty(char(sel_analysis.modeltree))
            set(handles.modeltree, 'String', '');
            set(handles.modellist, 'String', 'DEFINE A MODELTREE ABOVE');
            set(handles.modellist, 'Value', 1);
        else
            set(handles.modeltree, 'String', char(sel_analysis.modeltree));
            set(handles.modellist, 'Value', 1);
            rebuild_modeltree();
        end
    end


     function rebuild_status_filter()
         sql = ['SELECT DISTINCT status FROM NarfAnalysis AS sq'];
         ret = mysql(sql);
         statuses = cellstr(char(ret(:).status));
         set(handles.status_filter, 'String', cat(1, {'*'}, statuses));
     end

     function rebuild_tag_filter()
         sql = ['SELECT DISTINCT tags FROM NarfAnalysis AS sq'];
         ret = mysql(sql);
         tags = cellstr(char(ret(:).tags));
         % Split the tags by spaces or commas
         stags = cellfun(@(s) strtrim(strsep(s, ',')), tags, 'UniformOutput', false);        
         stags = cat(2, {'*'}, stags{:});
         stags = unique(stags);
         set(handles.tag_filter, 'String', stags);
     end
 
     function rebuild_analyses_table()            
         
         sql = ['SELECT * FROM NarfAnalysis'];
         
         isfirst = true;
         function s = whereand()
             if isfirst
                 isfirst = false;
                 s = ' WHERE ';
             else
                 s = ' AND ';
             end
         end 
         
         sel_status = popup2str(handles.status_filter);
         sel_tag = popup2str(handles.tag_filter);
         
         if ~isequal(sel_status, '*') 
             sql = [sql whereand() 'status="' sel_status '"'];
         end
         
         if ~isequal(sel_tag, '*') 
             sql = [sql whereand() 'tags LIKE "%' sel_tag '%"'];
         end
         
         analyses_found = mysql(sql);
         l = length(analyses_found);
         c = cell(l,3);
         for i = 1:l
             c{i,1} = analyses_found(i).id;
             c{i,2} = analyses_found(i).status;
             c{i,3} = analyses_found(i).name;  
         end
         set(handles.analyses_table, 'Data', c);
         drawnow; 
     end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Center panel      
uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Name:', ...
          'Position', [pad dh-(ts*1)-pad 50 ts-pad]);
      
handles.name = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Sample Analysis', ...
          'Position', [50+pad dh-(ts*1)-pad 200 ts], ...
          'Callback', @analysis_changed_callback);
      
uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Status:', ...
          'Position', [50+200+pad dh-(ts*1)-pad 50 ts-pad]);
      
handles.status = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'active', ...
          'Position', [50+200+50+pad dh-(ts*1)-pad cw-(50+200+50+2*pad) ts], ...
          'Callback', @analysis_changed_callback);

uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Q:', ...
          'Position', [pad dh-(ts*2)-pad 20 ts-pad]);
      
handles.short_desc = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'Max', 3, 'Min', 1, ...
          'HorizontalAlignment', 'left', 'String', 'Motivating question or experimental summary description goes here.', ...
          'Position', [20+pad dh-(ts*4) cw-20-pad ts*3-pad], ...
          'Callback', @analysis_changed_callback);    

uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...         
          'HorizontalAlignment', 'left', 'String', 'A:', ...
          'Position', [pad dh-(ts*6) 20 ts*2]);
      
handles.long_desc = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'Max', 10, 'Min', 1, ...
          'HorizontalAlignment', 'left', 'String', 'Answer or long description goes here.', ...
          'Position', [20+pad ts cw-20-pad dh-ts*5], ...
          'Callback', @analysis_changed_callback);   
    
uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...          
          'HorizontalAlignment', 'left', 'String', 'Tags:', ...
          'Position', [pad 0 50 ts-pad]);

handles.tags = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'unsorted, weird, goodresult, etc', ...
          'Position', [50+pad pad cw-50-pad*2 ts-pad], ...
          'Callback', @analysis_changed_callback);   
      
	% Whenever the analysis changes, update the database. 
    function analysis_changed_callback(~,~,~)
        
        % Throw a popup message if somebody else has updated the analysis
        % between when you last read it.  
        sql = ['SELECT * FROM NarfAnalysis WHERE ID=' num2str(sel_analysis.id)]; 
        ret = mysql(sql);
        if length(ret) ~= 1
            error('Selected analysis ID number not found\n');
        end
        if ret(1).lastmod ~= sel_analysis.lastmod;
            sel = questdlg('ERROR: analysis timestamp mismatch! Is somebody else editing this analysis too?', ...
                    'Timestamp Error', 'Discard Your Changes', 'Debug', 'Debug');
            if strcmp(sel, 'Debug')
                keyboard;
            else
                any_condition_changed_callback();
                %analyses_table_row_selected(findjobj(handles.analyses_table)); 
                return;
            end
        end
                
        % Otherwise, commit the changes to the database
        sql = ['UPDATE NarfAnalysis SET ' ...
               'name="' sql_sanitize(get(handles.name, 'String')) '", ' ...
               'status="' sql_sanitize(get(handles.status, 'String')) '", ' ...
               'question="' sql_sanitize(get(handles.short_desc, 'String')) '", ' ...
               'answer="' sql_sanitize(get(handles.long_desc, 'String')) '",' ...
               'tags="' sql_sanitize(get(handles.tags, 'String')) '", ' ...
               'modeltree="' sql_sanitize(get(handles.modeltree, 'String')) '", ' ...
               'batch="' sql_sanitize(get(handles.batch, 'String')) '" ' ...
               'WHERE id=' num2str(sel_analysis.id)]; 
        
        [r, affected] = mysql(sql);
        if affected ~= 1
            error('Why was the SQL database not updated!!?');
        end
        any_condition_changed_callback();
        %analyses_table_row_selected(findjobj(handles.analyses_table)); 
    end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Right Panel is for the model structures

uicontrol('Parent', right_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'ModelTree:', ...
          'Position', [pad dh-(ts*1)-pad 70 ts-pad]);
      
handles.modeltree = uicontrol('Parent', right_panel, 'Style', 'edit', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', '{''env100'', ''log2b'', ''firn'', {''npnl'', ''npfnl''}, ''boost''}', ...
          'Position', [pad dh-(ts*2) rw-pad*2 ts-pad], ...
          'Callback', @update_and_recalc_modeltree);

    function update_and_recalc_modeltree(~,~,~)
        analysis_changed_callback();
        rebuild_modeltree();
    end
handles.regen_modeltree = uicontrol('Parent', right_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Regenerate Model List', ...
          'Position', [rw-150-pad dh-(ts)-pad 150 ts], ...
          'Callback', @rebuild_modeltree); 
      
handles.modellist = uicontrol('Parent', right_panel, 'Style', 'listbox', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', '', ...
          'Max', 3, 'Min', 1, ...
          'Position', [pad dh-ts*7+pad rw-pad*2 5*ts-pad], ...
          'Callback', @any_condition_changed_callback);

handles.select_all_models = uicontrol('Parent', right_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Select All', ...
          'Position', [pad dh-(ts*8) 100 ts], ...
          'Callback', @any_condition_changed_callback); 
      
handles.select_none_models = uicontrol('Parent', right_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Select None', ...
          'Position', [100+pad dh-(ts*8) 100 ts], ...
          'Callback', @any_condition_changed_callback); 
      
uicontrol('Parent', right_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Batch:', ...
          'Position', [pad dh-(ts*9) 70 ts-pad]);
      
handles.batch = uicontrol('Parent', right_panel, 'Style', 'popupmenu', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', '*', ...
          'Position', [70 dh-(ts*9) rw-200-pad ts]);
      
handles.refresh_batch = uicontrol('Parent', right_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Refresh', ...
          'Position', [rw-100-pad dh-(ts*9) 100 ts], ...
          'Callback', @any_condition_changed_callback); 

handles.batchs_table = uitable('Parent', right_panel, 'Enable', 'on', 'Units', 'pixels',...
        'RowName', [], ...
        'ColumnWidth', {40, rw-310, 100, 100, 60}, ...
        'ColumnName', {'Show?', 'CellID', 'Est Set', 'Val Set', 'Filecodes'}, ...
        'Position', [pad pad+ts rw-pad*2 dh-10*ts-pad]);

handles.select_all_cellids = uicontrol('Parent', right_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Select All', ...
          'Position', [pad pad 100 ts], ...
          'Callback', @any_condition_changed_callback); 

handles.select_none_cellids = uicontrol('Parent', right_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Select None', ...
          'Position', [100+pad pad 100 ts], ...
          'Callback', @any_condition_changed_callback); 
      
    function rebuild_modeltree(~,~,~)
        try     
            mm = eval(get(handles.modeltree, 'String'));
            modulekeys = keyword_combos(mm);
            models = {};
            for ii = 1:length(modulekeys)
                tmp = cellfun(@(n) sprintf('%s_', n), modulekeys{ii}, ...
                    'UniformOutput', false);
                modelname = strcat(tmp{:});
                modelname = modelname(1:end-1);
                models{end+1} = modelname;
            end
            set(handles.modellist, 'String', models);
        catch            
            set(handles.modellist, 'String', {'ERROR IN MODELTREE EXPRESSION'});
        end 
    end      
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bottom Panel

uicontrol('Parent', bottom_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'SortBy:', ...
          'Position', [pad bh-ts-pad 50 ts-pad]);
      
handles.sortby = uicontrol('Parent', bottom_panel, 'Style', 'popupmenu', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'cellid', ...
          'Position', [50+pad bh-ts 150 ts-pad], ...
          'Callback', @any_condition_changed_callback);

handles.sortdir = uicontrol('Parent', bottom_panel, 'Style', 'popupmenu', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'ASC', ...
          'Position', [50+150+pad bh-ts 60 ts-pad], ...
          'Callback', @any_condition_changed_callback);

uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Scatter Plot', ...
          'Position', [300 bh-ts 100 ts-pad], ...
          'Callback', @any_condition_changed_callback);
      
uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Preview Model', ...
          'Position', [300 bh-ts 100 ts-pad], ...
          'Callback', @any_condition_changed_callback);
      
uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'View STRF', ...
          'Position', [400 bh-ts 100 ts-pad], ...
          'Callback', @any_condition_changed_callback);

uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'View CellDB', ...
          'Position', [500 bh-ts 100 ts-pad], ...
          'Callback', @any_condition_changed_callback);
      
uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Inspect Model', ...
          'Position', [600 bh-ts 100 ts-pad], ...
          'Callback', @any_condition_changed_callback);
      
uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Completion Report', ...
          'Position', [w-400 bh-ts 200 ts-pad], ...
          'Callback', @any_condition_changed_callback);
      
uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Enqueue Incomplete Models', ...
          'Position', [w-200 bh-ts 200-pad ts-pad], ...
          'Callback', @any_condition_changed_callback);
      
db_results_table = uitable('Parent', bottom_panel, ...
        'Enable', 'on',  'Units', 'pixels', 'RowName', [],...
        'ColumnWidth', {50, 40, 60, 50, 70, 200, 100, 150, 70, 70, 70, 70}, ...
        'ColumnName', {'ID', 'Batch', 'CellID', 'EstSet', 'ValSet', 'Model Name', 'Val R', 'Val NLL', 'Est R', 'Est NLL', 'Last Mod.', 'Notes'}, ...
        'Position', [pad pad w bh-ts-pad*2]);

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
        rebuild_status_filter();
        rebuild_tag_filter();
        rebuild_analyses_table();
        fprintf('Callback called\n');
    end

% Call the callback once to get things running
any_condition_changed_callback();

end



