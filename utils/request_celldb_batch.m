function ret = request_celldb_batch(batch, cellid)
% ret = request_celldb_batch(batch, cellid)
% 
% Returns a cell array of structs containing the cellids, training sets,
% and test set that correspond with a particular batch number. Each struct
% has the following fields:
%      .cellid
%      .training_set
%      .test_set
%      .filecode
%
% If the 'cellid' parameter is provided, then the cell array will contain
% only one struct corresponding to the matching cellid. 

global NARF_DEBUG

ret = {};

if ~exist('batch', 'var'),
   error('syntax error: request_celldb_batch(batch) parameters required');
end

if batch==243,
    celllist={...
        'out121012_002_002',...
        'out121112_004_003_E',...
        'out121112_004_003_I',...
        'out121112_004_005_E',...
        'out121112_004_005_I',...
        'out121412_007_002',...
        'out121812_002_003',...
        'out121812_002_005',...
        'out121812_004_003',...
        'out121812_004_005',...
        'out122012_002_001_E',...
        'out122012_002_001_I',...
        'out122012_003_001',...
        'out122012_007_001_E',...
        'out122012_007_001_I'};
    ret=cell(length(celllist),1);
    for ii=1:length(celllist),
        ret{ii}=struct();
        ret{ii}.cellid=celllist{ii};
        ret{ii}.training_set={[celllist{ii} '_est']};
        ret{ii}.test_set={[celllist{ii} '_val']};
        ret{ii}.filecode={};
    end
    return
end

dbopen;
if exist('cellid','var'),
    sql = ['SELECT * from sRunData WHERE batch=', num2str(batch),...
          ' AND cellid="',cellid,'"'];
else
    sql = ['SELECT * from sRunData WHERE batch=', num2str(batch)];
end
rundata = mysql(sql);

if length(rundata) == 0,
   fprintf('batch number %d not found!\n', batch);
   return
end


for nn = 1:length(rundata)
    cellid = rundata(nn).cellid;
    
    fprintf('CELLID=%s BATCH=%d\n', cellid, batch);

    train_set={};
    test_set={};
    file_code={};
    rawid=[];
    
    if ismember(rundata(nn).batch,[241 244 251 252])
        
        % find a file that matches specific behavioral selection
        % criteria for this batch
        cellfiledata=dbgettspfiles(cellid,rundata(nn).batch);
        activefile=strcmpi({cellfiledata.behavior},'active');
           
        acounter=1;
        pcounter=0;
        for ii=1:length(cellfiledata),
            train_set{ii}=[cellfiledata(ii).stimfile,'_est'];
            test_set{ii}=[cellfiledata(ii).stimfile,'_val'];
            rawid(ii)=cellfiledata(ii).rawid;
            
            if activefile(ii),
                file_code{ii}=['A',num2str(acounter)];
                acounter=acounter+1;
                pcounter=0;
            elseif ~pcounter,
                file_code{ii}=['P',num2str(acounter)];
                pcounter=1;
            else
                file_code{ii}=['P',num2str(acounter),char('a'-1+pcounter)];
                pcounter=pcounter+1;
            end                                             
        end
    else
        % figure out what files to use for what stage of the analysis
        %[cellfiledata, times, params] = cellfiletimes(cellid, rundata(nn).batch);
        % get all valid files
        cellfiledata=dbbatchcells(rundata(nn).batch,cellid);
        
        % sort so that file with max reps is at the end of the list
        % and will be used for validation/test
        ss=[cat(1,cellfiledata.stimspeedid),...
            cat(1,cellfiledata.repcount),...
            -cat(1,cellfiledata.resplen)];
        [~,sortcell]=sortrows(ss);
        cellfiledata=cellfiledata(sortcell);
        
        runclassids=cat(cellfiledata.runclassid);
        if sum(runclassids~=103)==0,
            cellfiledata=cellfiledata(1);
        end
        
        if length(cellfiledata)==1,
            train_set{1}=[cellfiledata.stimfile,'_est'];
            test_set{1}=[cellfiledata.stimfile,'_val'];
            rawid(1)=cellfiledata.rawid;
       else
            for ii=1:length(cellfiledata)-1,
                train_set{end+1}=cellfiledata(ii).stimfile;
                rawid(end+1)=cellfiledata(ii).rawid;
            end
            test_set{1}=cellfiledata(end).stimfile;
        end
    end
    
    ret{nn} = [];
    ret{nn}.cellid = cellid;
    ret{nn}.training_set = train_set;
    ret{nn}.test_set = test_set;
    ret{nn}.filecode = file_code;
    ret{nn}.training_rawid = rawid;
end

