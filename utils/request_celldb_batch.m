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

if batch==243 | batch==256,
    if batch==243,
        celllist={...
            'out121012_002_002',...
            'out121412_007_002',...
            'out121812_002_003',...
            'out121812_002_005',...
            'out121812_004_003',...
            'out121812_004_005',...
            'out122012_003_001',...
            'out012213_007_002',...
            'out031113_002_002',...
            'out040213_003_002',...
            'out040213_004_001',...
            'out051113_003_004',...
            'out051113_003_009',...
            'out051113_004_003',...
            'out061913_002_002',...
            'out062113_005_001',...
            'out121112_004_003_E',...
            'out121112_004_003_I',...
            'out121112_004_005_E',...
            'out121112_004_005_I',...
            'out122012_002_001_E',...
            'out122012_002_001_I',...
            'out122012_007_001_E',...
            'out122012_007_001_I',...
            'out022313_005_002_E',...
            'out022313_005_002_I',...
            'out022313_008_001_E',...
            'out022313_008_001_I',...
            'out052213_003_001_E',...
            'out052213_003_001_I',...
            'out052213_004_001_E',...
            'out052213_004_001_I',...
            'out071113_003_002_E',...
            'out071113_003_002_I',...
            'out071113_004_001_E',...
            'out071113_004_001_I',...
            'out082613_002_002_E',...
            'out082613_002_002_I',...
            'out082613_003_001_E',...
            'out082613_003_001_I',...
            'out082713_005_002_E',...
            'out082713_005_002_I',...
            'out082813_002_003_E',...
            'out082813_002_003_I',...
            'out121812_003_002_E',...
            'out121812_003_002_I',...
                 };
    elseif batch==256,
        celllist={...
            'out121112_004_003_E',...
            'out121112_004_005_E',...
            'out122012_002_001_E',...
            'out122012_007_001_E',...
            'out022313_005_002_E',...
            'out022313_008_001_E',...
            'out052213_003_001_E',...
            'out052213_004_001_E',...
            'out071113_003_002_E',...
            'out071113_004_001_E',...
            'out082613_002_002_E',...
            'out082613_003_001_E',...
            'out082713_005_002_E',...
            'out082813_002_003_E',...
            'out121812_003_002_E',...
                 };
    end
    celllist=sort(celllist);
    ret=cell(length(celllist),1);
    cc=0;
    for ii=1:length(celllist),
        if ~exist('cellid','var') ||...
                strcmp(cellid,celllist{ii}),
            cc=cc+1;
            ret{cc}=struct();
            ret{cc}.cellid=celllist{ii};
            if batch==256,
                ifile=strrep(celllist{ii},'_E','_I');
                ret{cc}.training_set={[celllist{ii} '_est'],[ifile '_est']};
                ret{cc}.test_set={[celllist{ii} '_val'],[ifile '_val']};
                ret{cc}.filecode={'E','I'};
            else
                ret{cc}.training_set={[celllist{ii} '_est']};
                ret{cc}.test_set={[celllist{ii} '_val']};
                ret{cc}.filecode={};
            end
        end
    end
    return
elseif batch==256
    
end

dbopen;
if exist('cellid','var'),
    sql = ['SELECT * from sRunData WHERE batch=', num2str(batch),...
          ' AND cellid="',cellid,'"'];
else
    sql = ['SELECT * from sRunData WHERE batch=', num2str(batch),...
           ' ORDER BY cellid'];
end
rundata = mysql(sql);

if length(rundata) == 0,
   fprintf('batch number %d not found!\n', batch);
   return
end


for nn = 1:length(rundata)
    cellid = rundata(nn).cellid;
    
    fprintf('request_celldb_batch.m: CELLID=%s BATCH=%d\n', cellid, batch);

    train_set={};
    test_set={};
    file_code={};
    rawid=[];
    
    if ismember(rundata(nn).batch,[241 244 251 252 253 254 255 258])
        
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
                res=spn_tuning_match(cellid,batch,0);
                baphyparms=dbReadData(rawid(ii));
                if isfield(baphyparms,'Tar_Frequencies'),
                    maxtar=length(baphyparms.Tar_Frequencies);
                else
                    maxtar=length(baphyparms.Trial_TargetIdxFreq);
                end
                snr=baphyparms.Trial_RelativeTarRefdB;
                targetfreq=baphyparms.Trial_TargetIdxFreq(1:maxtar);
                commontargetid=min(find(targetfreq==max(targetfreq)));
                if ismember(rundata(nn).batch,[251 253 254])
                    % batch 251, 253, 254:  A1 hard, A2 easy 
                    % for in/out data, A1: attend in
                    % for high/low data. A1: attemd low
                    commonHz=baphyparms.Tar_Frequencies(commontargetid);
                    uncommontarid=min(find(targetfreq==min(targetfreq)));
                    uncommonHz=baphyparms.Tar_Frequencies(uncommontarid);
                    commontarinref=length(baphyparms.Ref_LowFreq)==1 && ...
                        commonHz>=baphyparms.Ref_LowFreq &&....
                        commonHz<=baphyparms.Ref_HighFreq;
                    uncommontarinref=length(baphyparms.Ref_LowFreq)==1 && ...
                        uncommonHz>=baphyparms.Ref_LowFreq &&....
                        uncommonHz<=baphyparms.Ref_HighFreq;
                    if commontarinref && ~uncommontarinref,
                        Aidx=1;
                    elseif uncommontarinref && ~commontarinref,
                        Aidx=2;
                    elseif snr(commontargetid)==max(snr),
                        Aidx=2;
                    else
                        Aidx=1;
                    end
                    
                elseif ismember(rundata(nn).batch,[ 244  255])
                    % lr  batches
                    % A1: attend contra; A2: attend ipsi
                    if isfield(baphyparms,'Trial_TargetChannel'),
                        Aidx=baphyparms.Trial_TargetChannel(commontargetid);
                    elseif isfield(baphyparms,'Tar_SplitChannels') &&...
                                strcmpi(baphyparms.Tar_SplitChannels,'Yes')
                        Aidx=commontargetid
                    else
                        disp('unknown lr condition');
                        keyboard
                    end
                    
                elseif ismember(rundata(nn).batch,[241 252 255 258])
                    % lr/hl batches
                    % A1: attend band near BF, A2: attend away
                    
                    if max(res)>0 && res(commontargetid)==max(res),
                        Aidx=1;
                    else
                        Aidx=2;
                    end
                end
                
                file_code{ii}=['A',num2str(Aidx)];
                acounter=acounter+1;
                pcounter=0;
            elseif ~pcounter,
                file_code{ii}=['P',num2str(acounter)];
                pcounter=1;
            else
                % disabled lettering for multiple passives.  now
                % just collapse into a single set
                %file_code{ii}=['P',num2str(acounter),char('a'-1+pcounter)];
                file_code{ii}=['P',num2str(acounter)];
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

