function parent_handle = narf_analysis(parent_handle)
% A GUI to edit the contents of the MySQL NarfAnalysis table

global MODULES;
narf_set_path;
dbopen;
MODULES = scan_directory_for_modules();   

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
cellids_found = {};
models_found = {};
sel_analysis = [];
sel_models = {};
sel_batch = [];
sel_cellids = {};
sel_results = [];
preview_fig = [];

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
          'summaryfig', '', ...
          'batch', 'SELECT A BATCH');
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
            dbopen; mysql(sql);
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
set(hJTcb, 'MouseReleasedCallback', {@analyses_table_row_selected, gcf});
set(hJTcb, 'KeyPressedCallback', {@analyses_table_row_selected, gcf});

    function analyses_table_row_selected(a,~,~)
        r = a.getSelectedRows();
        sel_analysis = analyses_found(r+1);
        sel_models = {};
        sel_cellids = {};
        sel_results = [];
        
        set(handles.name, 'String', sel_analysis.name);
        set(handles.status, 'String', sel_analysis.status);
        set(handles.short_desc, 'String', char(sel_analysis.question));
        set(handles.long_desc, 'String', char(sel_analysis.answer));
        
        available_batches = get(handles.batch, 'String');
        if isempty(sel_analysis.batch) || isempty(strcmp(available_batches, sel_analysis.batch))
            set(handles.batch, 'Value', 1);
        else
            vec = strcmp(available_batches, sel_analysis.batch);
            idx = 1:length(vec);
            idx = idx(vec);
            if isempty(idx)
                idx(1) = 1;
            end
            set(handles.batch, 'Value', idx(1));
            if strcmp(sel_analysis.batch, 'SELECT A BATCH')
                sel_batch = [];
            else
                sel_batch = str2num(sel_analysis.batch);
            end
            rebuild_batch_table();
        end
        
        set(handles.tags, 'String', char(sel_analysis.tags));
        if isempty(char(sel_analysis.modeltree))
            set(handles.modeltree, 'String', '');
            set(handles.modellist, 'String', 'DEFINE A MODELTREE ABOVE');
            set(handles.modellist, 'Value', []);
        else
            set(handles.modeltree, 'String', char(sel_analysis.modeltree));
            set(handles.modellist, 'Value', []);
            rebuild_modeltree();
        end
        update_query_results_table();
    end


     function rebuild_status_filter()
         sql = ['SELECT DISTINCT status FROM NarfAnalysis AS sq'];
         dbopen;
         ret = mysql(sql);
         statuses = cellstr(char(ret(:).status));
         set(handles.status_filter, 'String', cat(1, {'*'}, statuses));
     end

     function rebuild_tag_filter()
         sql = ['SELECT DISTINCT tags FROM NarfAnalysis AS sq'];
         dbopen;
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
         
         dbopen;
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

    function rebuild_batches_popup()
        sql = ['SELECT DISTINCT batch FROM NarfBatches AS sq ORDER BY batch'];
        dbopen;
        ret = mysql(sql);       
        batches = cellstr(char(ret(:).batch));
        set(handles.batch, 'String', cat(1, {'SELECT A BATCH'}, batches));
    end

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Center panel      
uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Name:', ...
          'Position', [pad dh-(ts*1)-pad 50 ts-pad]);
      
handles.name = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', '', ...  
          'Max', 1, 'Min', 1, ...
          'Position', [50+pad dh-(ts*1)-pad 200 ts], ...
          'Callback', @analysis_changed_callback);
      
uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Status:', ...
          'Position', [50+200+pad dh-(ts*1)-pad 50 ts-pad]);
      
handles.status = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', '', ...
          'Position', [50+200+50+pad dh-(ts*1)-pad cw-(50+200+50+2*pad) ts], ...
          'Callback', @status_changed_callback);

    function status_changed_callback(~,~,~)  
        set(handles.status_filter, 'Value', 1);
        analysis_changed_callback;
    end

uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Q:', ...
          'Position', [pad dh-(ts*2)-pad 20 ts-pad]);
      
handles.short_desc = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'Max', 3, 'Min', 1, ...
          'HorizontalAlignment', 'left', 'String', '', ...
          'Position', [20+pad dh-(ts*4) cw-20-pad ts*3-pad], ...
          'Callback', @analysis_changed_callback);    

uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...         
          'HorizontalAlignment', 'left', 'String', 'A:', ...
          'Position', [pad dh-(ts*6) 20 ts*2]);
      
handles.long_desc = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'Max', 10, 'Min', 1, ...
          'HorizontalAlignment', 'left', 'String', '', ...
          'Position', [20+pad ts cw-20-pad dh-ts*5], ...
          'Callback', @analysis_changed_callback);   
    
uicontrol('Parent', center_panel, 'Style', 'text', 'Units', 'pixels',...          
          'HorizontalAlignment', 'left', 'String', 'Tags:', ...
          'Position', [pad 0 50 ts-pad]);

handles.tags = uicontrol('Parent', center_panel, 'Style', 'edit', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', '', ...
          'Position', [50+pad pad cw-50-pad*2 ts-pad], ...
          'Callback', @tag_changed_callback);
      
    function tag_changed_callback(~,~,~)  
        set(handles.tag_filter, 'Value', 1);
        analysis_changed_callback;
    end

	% Whenever the analysis changes, update the database. 
    function analysis_changed_callback(~,~,~)
        %  Do nothing if no analysis is selected
        if isempty(sel_analysis)
            return;
        end
        
        % Throw a popup message if somebody else has updated the analysis
        % between when you last read it.  
        sql = ['SELECT * FROM NarfAnalysis WHERE ID=' num2str(sel_analysis.id)]; 
        dbopen;
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
               'batch="' sql_sanitize(popup2str(handles.batch)) '" ' ...
               'WHERE id=' num2str(sel_analysis.id)]; 
        
        dbopen;
        [r, affected] = mysql(sql);
        % Commented out because sometimes you make NO changes to an entry,
        % so then you get no results affected.
%        if affected ~= 1 
%            error('Why was the SQL database not updated!!?');
        %end
        any_condition_changed_callback();
        update_query_results_table();
        %analyses_table_row_selected(findjobj(handles.analyses_table)); 
    end
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Right Panel is for the model structures

uicontrol('Parent', right_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'ModelTree:', ...
          'Position', [pad dh-(ts*1)-pad 70 ts-pad]);
      
handles.modeltree = uicontrol('Parent', right_panel, 'Style', 'edit', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', ' ', ...
          'Position', [pad dh-(ts*2) rw-pad*2 ts-pad], ...
          'Callback', @update_and_recalc_modeltree);

    function update_and_recalc_modeltree(~,~,~)
        analysis_changed_callback();
        rebuild_modeltree();
    end

handles.modellist = uicontrol('Parent', right_panel, 'Style', 'listbox', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', '', ...
          'Max', 3, 'Min', 1, ...
          'Position', [pad dh-ts*7+pad rw-pad*2 5*ts-pad], ...
          'Callback', @modellist_callback);
      
    function modellist_callback(~,~,~)  
        indexes = get(handles.modellist, 'Value');
        models = get(handles.modellist, 'String');
        sel_models = models(indexes);
        update_query_results_table();
    end

handles.select_all_models = uicontrol('Parent', right_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Select All', ...
          'Position', [rw-100-pad dh-(ts*1)-pad 100 ts], ...
          'Callback', @sel_all_models_callback); 
      
    function sel_all_models_callback(~,~,~)
        len = length(get(handles.modellist, 'String'));
        set(handles.modellist, 'Value', 1:len);        
        modellist_callback();
    end      
      
uicontrol('Parent', right_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Batch:', ...
          'Position', [pad dh-(ts*8) 50 ts-pad]);
      
handles.batch = uicontrol('Parent', right_panel, 'Style', 'popupmenu', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'SELECT A BATCH', ...
          'Position', [50+pad dh-(ts*8) 200 ts], ...
          'Callback', @batch_callback);
      
    function batch_callback(~,~,~)
        bat = popup2str(handles.batch);        
        sel_batch = str2num(bat);
        analysis_changed_callback();
        rebuild_batch_table();
    end
    
handles.refresh_batches = uicontrol('Parent', right_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Refresh Batch', ...
          'Position', [250+pad dh-(ts*8) 100 ts], ...
          'Callback', @refresh_batch_callback); 
      
    function refresh_batch_callback(~,~,~)
        if isempty(sel_batch)
            return
        end
        enable_or_disable_children(parent_handle, 'off');
        cells = request_celldb_batch(sel_batch);
        dbopen;
        sql = ['DELETE from NarfBatches WHERE batch="' num2str(sel_batch) '"'];
        mysql(sql);
        for ii = 1:length(cells)
            sqlinsert('NarfBatches', ...
                      'batch',      num2str(sel_batch), ...
                      'cellid',     cells{ii}.cellid, ...
                      'est_set',    write_readably(cells{ii}.training_set), ...
                      'val_set',    write_readably(cells{ii}.test_set), ...
                      'filecodes',  write_readably(cells{ii}.filecode));
        end
        rebuild_batch_table();
        enable_or_disable_children(parent_handle, 'on');
    end

handles.batch_table = uitable('Parent', right_panel, 'Enable', 'on', 'Units', 'pixels',...
        'RowName', [], ...
        'ColumnWidth', {rw-400, 150, 150, 60}, ...
        'ColumnName', {'CellID', 'Est Set', 'Val Set', 'Filecodes'}, ...
        'Position', [pad pad+ts rw-pad*2 dh-9*ts-pad]);

drawnow;
batch_JS = findjobj(handles.batch_table); 
batch_JT = batch_JS.getViewport.getView;
batch_JT.setNonContiguousCellSelection(false);
batch_JT.setColumnSelectionAllowed(false);
batch_JT.setRowSelectionAllowed(true);
batch_JTcb = handle(batch_JT, 'CallbackProperties');
set(batch_JTcb, 'MouseReleasedCallback', {@batch_table_row_selected, gcf});
set(batch_JTcb, 'KeyPressedCallback', {@batch_table_row_selected, gcf});

    function s = batch_table_row_selected(a, ~, ~)
        r = a.getSelectedRows();
        sel_cellids = cellids_found(r+1);
        update_query_results_table();
    end

    function s = remove_crap(s)
        s = regexprep(s, '''', '');
        s = regexprep(s, '{', '');
        s = regexprep(s, '}', '');
        s = regexprep(s, ',$', '');
    end
    
    function rebuild_batch_table(~,~,~)        
        sql = ['SELECT * from NarfBatches WHERE batch="' num2str(sel_batch) '" ORDER BY cellid'];
        dbopen;
        ret = mysql(sql);
        dat = cell(length(ret), 4);
        for ii = 1:length(ret)
            dat{ii, 1} = remove_crap(ret(ii).cellid);
            dat{ii, 2} = remove_crap(char(ret(ii).est_set));
            dat{ii, 3} = remove_crap(char(ret(ii).val_set));
            dat{ii, 4} = remove_crap(char(ret(ii).filecodes));
        end
        set(handles.batch_table, 'Data', dat);
        cellids_found = cellstr(char(ret(:).cellid));
        
        update_query_results_table();        
    end
    
handles.select_all_cellids = uicontrol('Parent', right_panel, ...
          'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Select All', ...
          'Position', [rw-100-pad dh-(ts*8) 100 ts], ...
          'Callback', @select_all_cellids_callback); 

    function select_all_cellids_callback(~,~,~)
        batch_JT.selectAll();
        batch_table_row_selected(batch_JT);
    end
      
    function rebuild_modeltree(~,~,~)
        try     
            mm = eval(get(handles.modeltree, 'String'));
            modulekeys = keyword_combos(mm);
            models_found = {};
            for ii = 1:length(modulekeys)
                tmp = cellfun(@(n) sprintf('%s_', n), modulekeys{ii}, ...
                    'UniformOutput', false);
                modelname = strcat(tmp{:});
                modelname = modelname(1:end-1);
                models_found{end+1} = modelname;
            end
            set(handles.modellist, 'String', models_found);
            update_query_results_table();
        catch            
            set(handles.modellist, 'String', {'ERROR IN MODELTREE EXPRESSION'});
        end 
    end

uicontrol('Parent', right_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Status Report', ...
          'Position', [pad pad 100 ts], ...
          'Callback', @status_report_callback);
    
    function yn = model_exists_in_db(bat, cellid, modelname, est_set, val_set, fcs) 
        sql = ['SELECT * FROM NarfResults WHERE batch=' num2str(bat) ...
               ' AND cellid="' cellid '" AND modelname="' modelname '"'];
           
        % TODO: Eventually this query should also test est, val, and fcs!
        ret = mysql(sql);
        yn = (1 == length(ret));
    end

    function status_report_callback(~,~,~)
%         questdlg('Are you sure you want to generate an analysis status report?', ...
%                  'Status Report?', 'Yes', 'No', 'No');
        
        enable_or_disable_children(parent_handle, 'off');
        fprintf('\n-------------------------------------------------------------------------------\n');             
        fprintf('Model/CellID Completion Matrix; an X indicates the model exists already.');
        fprintf('\n-------------------------------------------------------------------------------\n');
        fprintf('%-20s|1|2|3|4|5|6|7|8|9|0|1|2|3|4|5|6|7|8|9|0|1|2|3|4|5|6|7|8|9|0\n', 'CELLID');
        
        dbopen;
        total = 0;
        complete = 0;
        for ii = 1:length(cellids_found)
            fprintf('%-20s|', cellids_found{ii});
            for jj = 1:length(models_found)       
                if model_exists_in_db(sel_batch, cellids_found{ii}, models_found{jj})
                    fprintf('X|');
                    complete = complete+1;
                else
                    fprintf(' |');
                end
                total = total+1;
            end
            fprintf('\n');
        end
        
        fprintf('-------------------------------------------------------------------------------\n');        
        fprintf('CellIDs: %d\n', length(cellids_found));
        fprintf('Models: %d\n', length(models_found));
        fprintf('Status: [%d/%d] (complete/total), %d models not yet processed.\n', complete, total, total - complete);       
        fprintf('-------------------------------------------------------------------------------\n');        

        enable_or_disable_children(parent_handle, 'on');
    end

uicontrol('Parent', right_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Fit Model Now', ...
          'Position', [100+pad pad 100 ts], ...
          'Callback', @fit_selected_now);
      
    function fit_selected_now(~,~,~)        
        if length(sel_models) ~= 1 || length(sel_cellids) ~= 1               
            warndlg('Only a single model and cellid must be selected.');
        end
        
        reply = questdlg('Are you sure you want to fit the selected model to the selected cellid, on this machine only, immediately? While it is running, please do not click on any of the GUI buttons.', ...
                 'Fit this model now?', 'Yes', 'No', 'No');
        if strcmp(reply, 'No')
            return;
        end
        enable_or_disable_children(parent_handle, 'off');
        mm = eval(get(handles.modeltree, 'String'));
        modulekeys = keyword_combos(mm);
        
        cells = request_celldb_batch(sel_batch, sel_cellids{1});       
        ii=1;
        
        indexes = get(handles.modellist, 'Value');        
        jj=indexes(1);
        
        try 
            fit_single_model(sel_batch, cells{ii}.cellid, modulekeys{jj}, ...
            cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode, false);     
        catch err
            enable_or_disable_children(parent_handle, 'on');    
            rethrow(err);
        end
        
        enable_or_disable_children(parent_handle, 'on');    
    end

handles.force = uicontrol('Parent', right_panel, 'Style', 'checkbox', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Force?', ...
          'Position', [rw-220+pad pad 70-pad ts]);

uicontrol('Parent', right_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Enqueue Models', ...
          'Position', [rw-150-pad pad 150-pad ts], ...
          'Callback', @enqueue_models_callback);
      
    function enqueue_models_callback(~,~,~)
        reply = questdlg('Have you pushed your code to GIT? Are you sure you want to enqueue POSSIBLY MANY models? Have you run a status report first?', ...
                 'STOP AND THINK A MOMENT!', 'Yes', 'No', 'No');
        if strcmp(reply, 'No')
            return;
        end
        enable_or_disable_children(parent_handle, 'off');
        thebatch = sel_batch;        
        force = get(handles.force, 'Value');
        
        mm = eval(get(handles.modeltree, 'String'));
        modulekeys = keyword_combos(mm);
        
        cells = request_celldb_batch(thebatch);
        
        for ii = 1:length(cells)
            for jj = 1:length(modulekeys)                     
                enqueue_single_model(thebatch, cells{ii}.cellid, modulekeys{jj}, ...
                    cells{ii}.training_set, cells{ii}.test_set, cells{ii}.filecode, force);
            end
        end
        enable_or_disable_children(parent_handle, 'on');
             
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bottom Panel

ButtonWidth=75;

uicontrol('Parent', bottom_panel, 'Style', 'text', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'SortBy:', ...
          'Position', [pad bh-ts-pad 50 ts-pad]);
              
handles.sortby = uicontrol('Parent', bottom_panel, 'Style', 'popupmenu', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', {'batch', 'cellid', 'id', 'lastmod', 'modelname' 'r_fit', 'r_test', 'sparsity'}, ...
          'Value', 2, 'Position', [50+pad bh-ts 150 ts-pad], ...
          'Callback', @sort_callback);

    function sort_callback(~,~,~)
        update_query_results_table();
        any_condition_changed_callback();
    end

handles.sortdir = uicontrol('Parent', bottom_panel, 'Style', 'popupmenu', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', {'ASC', 'DESC'}, ...
          'Value', 1, 'Position', [50+150+pad bh-ts 60 ts-pad], ...
          'Callback', @sort_callback);

uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Preview', ...
          'Position', [300 bh-ts ButtonWidth ts-pad], ...
          'Callback', @preview_model_callback);
    function clear_preview_fig_handle(f,~)
        preview_fig = [];
        close(f);
    end
      
    function preview_model_callback(~,~,~)
        if ~isempty(sel_results)
            for ii = 1:length(sel_results)
                preview_fig = figure('Menubar', 'none', ...
                                     'CloseRequestFcn', @clear_preview_fig_handle);
                [I,map] = imread(char(sel_results(ii).figurefile),'png');  
                imshow(I, map);
            end
        end
    end

uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'STRF', ...
          'Position', [300+ButtonWidth*1 bh-ts ButtonWidth ts-pad], ...
          'Callback', @view_strf_button_callback);

    function view_strf_button_callback(~, ~, ~)
        if isempty(sel_results)
            return
        end
        for ii = 1:length(sel_results)
            [cfd, ~, ~] = dbgetscellfile('cellid', sel_results(ii).cellid);
            for jj = 1:length(cfd);
                % TODO: Replace magic number  1 with better description of TOR files
                if (cfd(jj).runclassid == 1)
                    figure;
                    strf_offline2([cfd(jj).stimpath cfd(jj).stimfile], ...
                        [cfd(jj).path cfd(jj).respfile], ...
                        cfd(jj).channum, cfd(jj).unit);
                end
            end
        end
    end

 
uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'CellDB', ...
          'Position', [300+ButtonWidth*2 bh-ts ButtonWidth ts-pad], ...
          'Callback', @launch_celldb_in_browser);
    
    function launch_celldb_in_browser(~,~,~)
        if isempty(sel_results)
            return
        end
        for ii = 1:length(sel_results)
            web(['http://hyrax.ohsu.edu/celldb/peninfo.php?penname=' ...
                 sel_results(ii).cellid(1:6) '#' sel_results(ii).cellid(1:7)]);
        end
    end

uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Inspect', ...
          'Position', [300+ButtonWidth*3 bh-ts ButtonWidth ts-pad], ...
          'Callback', @inspect_model_callback);
      
    function inspect_model_callback(~,~,~)
        if isempty(sel_results) || length(sel_results) > 1
            return
        end
        load_model(char(sel_results(1).modelpath));
        narf_modelpane();
    end

uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Scatter Plot', ...
          'Position', [300+ButtonWidth*4 bh-ts ButtonWidth ts-pad], ...
          'Callback', @scatter_plot_callback);
      
    function data = compute_data_matrix(fieldstoget)
        if isempty(db_results) || isempty(sel_batch) || isempty(sel_cellids) || isempty(sel_models)
            data = [];
            return;
        end
        dbopen;                
        
        enable_or_disable_children(parent_handle, 'off');
        data = nan(length(sel_cellids), length(sel_models));        
        for ii = 1:length(sel_models)
            fprintf('Querying DB for model: %s\n', sel_models{ii});
            for jj = 1:length(sel_cellids)
                sql = ['SELECT * FROM NarfResults WHERE batch=' num2str(sel_batch) ''];
                sql = [sql ' AND cellid in (' interleave_commas(sel_cellids(jj)) ')'];
                sql = [sql ' AND modelname in (' interleave_commas(sel_models(ii)) ')'];
                ret = mysql(sql); 
                
                if isempty(ret) 
                    continue;
                end
                if length(ret) > 1
                    error('More than one matching model found!');
                end
                
                for kk = 1:length(fieldstoget)
                    data(jj, ii, kk) = ret(1).(fieldstoget{kk});
                end
            end
        end
        enable_or_disable_children(parent_handle, 'on');
    end
      
    function scatter_plot_callback(~,~,~)        
        data = compute_data_matrix({'r_test'});
        if isempty(data)
            return;
        end
        % Query the DB and record a single data point
        figure('Name', 'Validation Set Comparison', 'NumberTitle', 'off', ...
               'Position', [10 10 900 900]);
        plot_scatter(data, sel_models); 
    end

    uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Heat Map', ...
          'Position', [300+ButtonWidth*5 bh-ts ButtonWidth ts-pad], ...
          'Callback', @heatmap_plot_callback);
      
    function heatmap_plot_callback(~,~,~)
        data = compute_data_matrix({'r_test'});
        if isempty(data)
            return;
        end
        figure('Name', 'Heat Map Comparison', 'NumberTitle', 'off', ...
               'Position', [10 10 900 900]);
        heatmap(data', sel_cellids, sel_models,  '', 'TickAngle', 90, 'ShowAllTicks', true); 
    end

    uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Bar Plot', ...
          'Position', [300+ButtonWidth*6 bh-ts ButtonWidth ts-pad], ...
          'Callback', @bar_plot_callback);
      
    function enable_or_disable_children(stuff, value)
        for ii = 1:length(stuff)
            h = stuff(ii);
            if isprop(h, 'Enable')
                set(h, 'Enable', value);
                drawnow;
            end
            children = get(stuff, 'children');
            if ~isempty(children);
                for jj = 1:length(children)
                    enable_or_disable_children(children(jj), value);                    
                end
            end
        end
    end
      
    function bar_plot_callback(~,~,~)
        data = compute_data_matrix({'r_test'});
        if isempty(data)
            return;
        end
        % Sort the data
        D = [nanmean(data)' data'];
        [sD, idxs] = sortrows(D, -1);
        data = sD(:, 2:end)';
        figure('Name', 'Bar Plot Comparison', 'NumberTitle', 'off', ...
               'Position', [10 10 900 300]);
        hold on;
        len = length(sel_models);
        bar(1:len, nanmean(data), 'r'); 
        plot(1:len, data, 'k.');
        errorbar(1:len, nanmean(data), nanvar(data), max(data) - nanmean(data), 'xk');
        hold off;
        set(gca,'XTick', 1:len);
        set(gca,'XTickLabel', sel_models(idxs));
        set(gca,'CameraUpVector',[-1,0,0]);
    end

    uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Rank Plot', ...
          'Position', [300+ButtonWidth*7 bh-ts ButtonWidth ts-pad], ...
          'Callback', @elite_plot_callback);
      
    function elite_plot_callback(~,~,~)
        data = compute_data_matrix({'r_fit', 'r_test'});
        if isempty(data)
            return;
        end
        
        vals = data(:,:,2);
        num_neurons = size(data,1);
        
        len = length(sel_models);
        val_rank = zeros(size(vals));
        for ii = 1:size(vals, 1) 
            [~, idxs] = sort(vals(ii,:), 2, 'descend');
            val_rank(ii,idxs) = 1:len;
        end
        
        score = nanmean(val_rank);
        [~, order] = sort(score, 'ascend');
        
        Ys = zeros(len);
        for ii = 1:len
            Ys(ii,:) = hist(val_rank(:, order(ii)), 1:len);
        end
        
        figure('Name', 'Est/Val Comparison', 'NumberTitle', 'off', ...
               'Position', [10 10 1000 500]);
        bar(Ys,'stacked');
        set(gca,'XTick', 1:length(sel_models));
        set(gca,'XTickLabel', sel_models(order));
        set(gca,'CameraUpVector',[-1,0,0]);
        legend(cellfun(@num2str, num2cell(1:len), 'UniformOutput', false), ...
            'Location','EastOutside');
        title(sprintf('Val. Set Ranking of Models (Blue is good) [%d cellids]', num_neurons));
        axis tight;
        %xlabel('CellID');
        %ylabel('Rankings'); 
    end

    CustomFunctionHandle=uicontrol('Parent', parent_handle, 'Style', 'edit',...
                'String', 'dummy',...
                'Units','pixels', 'Position', [300+ButtonWidth*9  bh-ts ButtonWidth*2 ts-pad], ...
                'Callback', @custom_model_analysis);

    uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
          'HorizontalAlignment', 'left', 'String', 'Custom:', ...
          'Position', [300+ButtonWidth*8 bh-ts ButtonWidth ts-pad], ...
          'Callback', @custom_model_analysis);

function custom_model_analysis(~,~,~)
    global STACK XXX META;
    if isempty(sel_results) || length(sel_results) > 1
        return
    end
    cust_function_name=get(CustomFunctionHandle,'String');
    load_model(char(sel_results(1).modelpath));
    feval(cust_function_name,STACK,XXX,META,sel_results);
end

    
    
    
db_results_table = uitable('Parent', bottom_panel, ...
        'Enable', 'on',  'Units', 'pixels', 'RowName', [],...
        'ColumnWidth', {60, 40, 100, 300, ...
                        10, 10, ...
                        60, 60, 60, 60, ...
                        60, 60, 60, ...
                        150, 50}, ...
        'ColumnName', {'ID', 'Batch', 'CellID', 'Modelname', ...
                       'est_set', 'val_set', ...
                       'val_corr', 'val_l1err', 'val_l2err', 'val_nlogl',...
                       'est_corr', 'Sparse', 'Smooth', ...
                       'Last Mod.', 'Note'}, ...
        'Position', [pad pad w-pad*2 bh-ts-pad*2]);
    
% Set up the DB Results table widget behavior
hJScroll = findjobj(db_results_table); 
hJTable = hJScroll.getViewport.getView; 
hJTable.setNonContiguousCellSelection(false);
hJTable.setColumnSelectionAllowed(false);
hJTable.setRowSelectionAllowed(true);
hJTablecb = handle(hJTable, 'CallbackProperties');
set(hJTablecb, 'MouseReleasedCallback', {@get_selected_row, gcf});
set(hJTablecb, 'KeyPressedCallback', {@get_selected_row, gcf});

    function get_selected_row(a,~,~)
        r = a.getSelectedRows();
        sel_results = db_results(r+1);
        if ~isempty(sel_results) && any(ishandle(preview_fig))
            sfigure(preview_fig);
            [I,map] = imread(char(sel_results(1).figurefile),'png');  
            imshow(I, map);
        end
    end

    function s = interleave_commas(c)
        s = ['"' c{1} '"']; 
        for ii = 2:length(c)
            s = cat(2, s, [', "' c{ii} '"']);
        end
    end
    function sql = update_query_results_table()
        if isempty(sel_batch) || isempty(sel_cellids) || isempty(sel_models)
            db_results = [];
            set(db_results_table, 'Data', {});
            sel_results = [];
            drawnow;
            return
        end        
         
        sql = ['SELECT * FROM NarfResults WHERE batch=' num2str(sel_batch) ''];
        sql = [sql ' AND cellid in (' interleave_commas(sel_cellids) ')'];
        sql = [sql ' AND modelname in (' interleave_commas(sel_models) ')'];      

        sortby = popup2str(handles.sortby);
        sortdir = popup2str(handles.sortdir);
        sql = [sql ' ORDER BY ' sortby ' ' sortdir ' LIMIT 0, 500'];
        
        dbopen;        
        db_results = mysql(sql);       
        
        l = length(db_results);
        c = cell(l,12);
        for i = 1:l
            c{i,1} = db_results(i).id;
            c{i,2} = db_results(i).batch; 
            c{i,3} = char(db_results(i).cellid);
            c{i,4} = char(db_results(i).modelname);
            %c{i,5} = char(db_results(i).est_set);
            %c{i,6} = char(db_results(i).val_set);
            %c{i,7} = db_results(i).val_corr;
            c{i,7} = db_results(i).r_test;
            %c{i,8} = db_results(i).val_l1err;
            %c{i,9} = db_results(i).val_l2err;
            %c{i,10} = db_results(i).val_nlogl;
            %c{i,11} = db_results(i).est_corr;
            c{i,11} = db_results(i).r_fit;
            c{i,12} = db_results(i).sparsity;
            %c{i,13} = db_results(i).smoothness;
            c{i,14} = db_results(i).lastmod; 
            %c{i,15} = char(db_results(i).notes);           
        end
        set(db_results_table, 'Data', c);
        sel_results = [];
        drawnow;        
    end

    function any_condition_changed_callback(~,~,~)
        rebuild_status_filter();
        rebuild_tag_filter();
        rebuild_analyses_table();
    end

% Call the callbacks once to get things running
pause(0.1);
rebuild_batches_popup();
any_condition_changed_callback();

end
