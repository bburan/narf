function cellfiledata=dbgettspfiles(cellid,batch)

    global NARF_DEBUG
    
% find a file that matches specific behavioral selection
% criteria for this batch
[batchactivefiles,~,params]=dbbatchcells(batch,cellid);
firstactiveidx=[];
params.activeonly=getparm(params,'activeonly',0);

cellfiledata=dbgetscellfile('cellid',cellid,'runclassid',[42 8 103]);
filecount=length(cellfiledata);

% special case, find only matching behavior sets
%runclassid=cat(2,cellfiledata.runclassid);
%stimspeedid=cat(2,cellfiledata.stimspeedid);
activefile=zeros(1,length(cellfiledata));
stimprint=zeros(filecount,5);
for ii=1:length(cellfiledata),
    [parms, ~]=dbReadData(cellfiledata(ii).rawid);
    if cellfiledata(ii).rawid==batchactivefiles(1).rawid,
        % first active that matches criteria for this batch
        firstactiveidx=ii;
    end
    
    if strcmpi(cellfiledata(ii).behavior,'active'),
        activefile(ii)=1;
    end
    
    % "fingerprint" is subset and frequency range.  Need to match
    % between passive and active
    if params.activeonly,
        thisprint=[parms.Ref_Subsets parms.Ref_LowFreq ...
                   parms.Ref_HighFreq activefile(ii)];
    else
        thisprint=[parms.Ref_Subsets parms.Ref_LowFreq parms.Ref_HighFreq];
    end
    stimprint(ii,1:length(thisprint))=thisprint;
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

 
