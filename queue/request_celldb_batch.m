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

if batch==243 || batch==256,
    celllist=wehr_db;
    celllist=sort(celllist);
    cc=0;;
    
    if batch==243,
        ret=cell(length(celllist),1);
        for ii=1:length(celllist),
            if ~exist('cellid','var') || strcmp(cellid,celllist{ii}),
                cc=cc+1;
                ret{cc}=struct();
                ret{cc}.cellid=celllist{ii};
                ret{cc}.training_set={[celllist{ii} '_est']};
                ret{cc}.test_set={[celllist{ii} '_val']};
                ret{cc}.filecode={};
            end
        end
        
    elseif batch==256
        realcellids=celllist;
        for ii=1:length(celllist),
            tc=strsep(celllist{ii},'_',1);
            realcellids{ii}=[tc{1} '_' tc{2}];
        end
        ucellids=unique(realcellids);
        ret={};
        for ii=1:length(ucellids),
            filematches=find(strcmp(realcellids,ucellids{ii}));
            
            if (~exist('cellid','var') || strcmp(cellid,ucellids{ii}))...
                    && length(filematches)>1,
                cc=cc+1;
                ret{cc,1}=struct();
                ret{cc}.cellid=ucellids{ii};
                for ff=1:length(filematches),
                    fn=celllist{filematches(ff)};
                    ret{cc}.training_set{ff}=[fn '_est'];
                    ret{cc}.test_set{ff}=[fn '_val'];
                    if fn(end)=='E',
                        ret{cc}.filecode{ff}='E';
                    elseif fn(end)=='I',
                        ret{cc}.filecode{ff}='I';
                    else
                        ret{cc}.filecode{ff}='C';
                    end
                end
             end
        end
    end
    
    return
    
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

if isempty(rundata),
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
    
    if ismember(rundata(nn).batch,[241 244 251 252 253 254 255 258 268])
        
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
                        Aidx=commontargetid;
                    else
                        disp('unknown lr condition');
                        keyboard
                    end
                    
                elseif ismember(rundata(nn).batch,[241 252 255 258 268])
                    % lr/hl batches
                    % A1: attend band near BF, A2: attend away
                    if isfield(baphyparms,'Trial_TargetChannel'),
                        commontargetchan=baphyparms.Trial_TargetChannel(commontargetid);
                    else
                        commontargetchan=commontargetid;
                    end
                    if max(res)>0 && res(commontargetchan)==max(res),
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
    elseif ismember(rundata(nn).batch,[262])
        % RDT data
        cellfiledata=dbbatchcells(rundata(nn).batch,cellid);
        activefile=strcmpi({cellfiledata.behavior},'active');
           
        acounter=1;
        pcounter=0;
        for ii=1:length(cellfiledata),
            train_set{ii}=[cellfiledata(ii).stimfile,'_est'];
            test_set{ii}=[cellfiledata(ii).stimfile,'_val'];
            rawid(ii)=cellfiledata(ii).rawid;
            
            if activefile(ii),
                Aidx=1;
                file_code{ii}=['A',num2str(Aidx)];
             elseif ~pcounter,
                file_code{ii}=['P',num2str(acounter)];
                pcounter=1;
            else
                file_code{ii}=['P',num2str(acounter)];
                pcounter=pcounter+1;
            end                            
        end
    elseif ismember(rundata(nn).batch,[263])
        % noisy vocalization data
        cellfiledata=dbbatchcells(rundata(nn).batch,cellid);
        snr=cat(1,cellfiledata.stimsnr);
        repcount=cat(1,cellfiledata.reps);
        ff=find(snr>=100);
        if length(ff)==1,
           train_set={[cellfiledata(ff).stimfile,'_est']};
           test_set={[cellfiledata(ff).stimfile,'_val']};
           rawid=cellfiledata(ff).rawid;
           test_rawid=cellfiledata(ff).rawid;
           file_code={'C1'};
        else
           ffm=ff(min(find(repcount(ff)==max(repcount(ff)))));
           train_set={cellfiledata(setdiff(ff,ffm)).stimfile};
           test_set={cellfiledata(ffm).stimfile};
           rawid=cat(1,cellfiledata(setdiff(ff,ffm)).rawid);
           test_rawid=cellfiledata(ffm).rawid;
           file_code=repmat({'C1'},[1 length(ff)-1]);
        end
        test_file_code={'C1'};

        ff=find(snr<100);
        if length(ff)==1,
           train_set={train_set{:} [cellfiledata(ff).stimfile,'_est']};
           test_set={test_set{:} [cellfiledata(ff).stimfile,'_val']};
           rawid=cat(1,rawid,cellfiledata(ff).rawid);
           test_rawid=cat(1,test_rawid,cellfiledata(ff).rawid);
           file_code=cat(2,file_code,{'N1'});
        else
           ffm=ff(min(find(repcount(ff)==max(repcount(ff)))));
           train_set={train_set{:} cellfiledata(setdiff(ff,ffm)).stimfile};
           test_set={test_set{:} cellfiledata(ffm).stimfile};
           rawid=cat(1,rawid,cellfiledata(setdiff(ff,ffm)).rawid);
           test_rawid=cat(1,test_rawid,cellfiledata(ffm).rawid);
           file_code=cat(2,file_code,repmat({'N1'},[1 length(ff)-1]));
        end
        test_file_code={'C1','N1'};
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
            test_rawid=cellfiledata(end).rawid;
        end
    end
    
    if ~exist('test_file_code','var')
       test_file_code=file_code;
    end
    if ~exist('test_rawid','var')
       test_rawid=rawid;
    end

    ret{nn} = [];
    ret{nn}.cellid = cellid;
    ret{nn}.training_set = train_set;
    ret{nn}.test_set = test_set;
    ret{nn}.filecode = file_code;
    ret{nn}.test_filecode = test_file_code;
    ret{nn}.training_rawid = rawid;
    ret{nn}.test_rawid = test_rawid;
end

