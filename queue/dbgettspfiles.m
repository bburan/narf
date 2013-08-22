function cellfiledata=dbgettspfiles(cellid,batch)

    global NARF_DEBUG
    
% find a file that matches specific behavioral selection
% criteria for this batch
[batchactivefiles,~,params]=dbbatchcells(batch,cellid);
firstactiveidx=[];

cellfiledata=dbgetscellfile('cellid',cellid,'runclassid',[42 8 103]);
filecount=length(cellfiledata);

% special case, find only matching behavior sets
%runclassid=cat(2,cellfiledata.runclassid);
%stimspeedid=cat(2,cellfiledata.stimspeedid);
activefile=zeros(1,length(cellfiledata));
stimprint=zeros(filecount,5);
prestimsilence=zeros(filecount,1);
duration=zeros(filecount,1);
poststimsilence=zeros(filecount,1);
tarband=zeros(filecount,1);
ForcePreStimSilence=0.25;
for ii=1:length(cellfiledata),
    [parms, ~]=dbReadData(cellfiledata(ii).rawid);
    if cellfiledata(ii).rawid==batchactivefiles(1).rawid,
        % first active that matches criteria for this batch
        firstactiveidx=ii;
    end
    
    if strcmpi(cellfiledata(ii).behavior,'active'),
        activefile(ii)=1;
    end
    if isfield(parms,'Trial_TargetIdxFreq'),
        tif=parms.Trial_TargetIdxFreq;
        tarband(ii) = find(tif==max(tif), 1);
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

if isempty(firstactiveidx),
    error('no active TSP data for this cell');
end
printmatch=double(sum(abs(...
    stimprint-repmat(stimprint(firstactiveidx,:),filecount,1)),2)==0);
useidx=find(printmatch);
if NARF_DEBUG,
    iso=cat(1,cellfiledata.isolation);
    useidx=find(printmatch & iso>params.miniso);
end

cellfiledata=cellfiledata(useidx);

 