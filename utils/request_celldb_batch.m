function ret = request_celldb_batch(batch)
% Exploiting all the code written by Stephen, find the cellids, training,
% and test sets that correspond with a particular batch number. Returns a
% cell array of structs with three fields per struct:
%      .cellid
%      .training_set
%      .test_set

ret = {};

if ~exist('batch'),
   disp('syntax error: request_celldb_batch(batch) parameters required');
   return
end

dbopen;
sql = ['SELECT * from sRunData WHERE batch=', num2str(batch)];
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
    
    if rundata(nn).batch==241,
        cellfiledata=dbgetscellfile('cellid',cellid,'runclassid',[42 8 103]);
        filecount=length(cellfiledata);
        
        % special case, find only matching behavior sets
        runclassid=cat(2,cellfiledata.runclassid);
        stimspeedid=cat(2,cellfiledata.stimspeedid);
        activefile=zeros(1,length(cellfiledata));
        stimprint=zeros(filecount,5);
        prestimsilence=zeros(filecount,1);
        duration=zeros(filecount,1);
        poststimsilence=zeros(filecount,1);
        tarband=zeros(filecount,1);
        ForcePreStimSilence=0.25;
        for ii=1:length(cellfiledata),
            [parms,perf]=dbReadData(cellfiledata(ii).rawid);
            
            if strcmpi(cellfiledata(ii).behavior,'active'),
                activefile(ii)=1;
            end
            if isfield(parms,'Trial_TargetIdxFreq'),
                tif=parms.Trial_TargetIdxFreq;
                tarband(ii)=min(find(tif==max(tif)));
            end
            
            % "fingerprint" is subset and frequency range.  Need to match
            % between passive and active
            thisprint=[parms.Ref_Subsets parms.Ref_LowFreq parms.Ref_HighFreq];
            stimprint(ii,1:length(thisprint))=thisprint;
            %prestimsilence(ii)=parms.Ref_PreStimSilence;
            prestimsilence(ii)=ForcePreStimSilence;
            
            duration(ii)=parms.Ref_Duration;
            poststimsilence(ii)=parms.Ref_PostStimSilence;
            poststimsilence(ii)=0; % force to zero below
        end
        firstactiveidx=max(find(activefile));  % use last
                                               % activefile!
        if isempty(firstactiveidx),
            error('no active TSP data for this cell');
        end
        printmatch=double(sum(abs(stimprint-repmat(stimprint(firstactiveidx,:),filecount,1)),2)==0);
        
        useidx=find(printmatch);
        train_set={cellfiledata(useidx).stimfile};
        
    else
        % figure out what files to use for what stage of the analysis
        [cellfiledata, times, params] = cellfiletimes(cellid, rundata(nn).batch);
    
        for ii=1:length(times(1).fileidx),
            if times(1).stop(ii)>times(1).start(ii),
                train_set{end+1}=basename(params.stimfiles{times(1).fileidx(ii)});
            end
        end
        
        for ii=1:length(times(3).fileidx),
            if times(3).stop(ii)>times(3).start(ii),
                test_set{end+1}=basename(params.stimfiles{times(3).fileidx(ii)});
            end
        end
    end
    
    ret{nn} = [];
    ret{nn}.cellid = cellid;
    ret{nn}.training_set = train_set;
    ret{nn}.test_set = test_set;
    
end

