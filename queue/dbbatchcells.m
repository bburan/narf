function [cellfiledata,cellids,params]=dbbatchcells(batchid,cellid);

dbopen;

if isnumeric(batchid),
   sql=['SELECT * FROM sBatch WHERE id=',num2str(batchid)];
   params=mysql(sql);
   if length(params.parmstring)>0,
      eval(char(params.parmstring));
   end
else
   params=batchid;
end

cellargs={};

if exist('cellid','var'),
   cellargs=cat(2,cellargs,{'cellid'},{cellid});
end
params.miniso=getparm(params,'miniso',0);
params.activeonly=getparm(params,'activeonly',0);

if 1 || params.stimspeedid==0,
   % disp('disabling speed test');
elseif params.stimspeedid>=60,
   cellargs=cat(2,cellargs,{'speedgt'},{params.stimspeedid-1});
else,
   cellargs=cat(2,cellargs,{'speed'},{params.stimspeedid});
end

if isfield(params,'stimsnr') && ~isempty(params.stimsnr),
   cellargs=cat(2,cellargs,{'stimsnr'},{params.stimsnr});
end

if isfield(params,'area') && ~isempty(params.area),
   %fprintf('restricted to cells matching area="%s"\n',params.area);
   cellargs=cat(2,cellargs,{'area'},{params.area});
end

if isfield(params,'dataparm'),
   cellargs=cat(2,cellargs,params.dataparm);
end

if strcmp(params.resploadcmd,'loadgammaraster'),
   % lfp data must exist.  also only need one unit per channel!
   cellargs=cat(2,cellargs,{'lfp'},{1});
end

if ismember(batchid,[140 143]),
   cellargs=cat(2,cellargs,{'runclassid'},{32});
else
   cellargs=cat(2,cellargs,{'runclassid'},{params.runclassid});
end

cellargs=cat(2,cellargs,{'respfmtcode'},{params.respfmtcode});
%cellargs=cat(2,cellargs,{'stimfmtcode'},{params.stimfmtcode});

switch batchid,
 case {139,140,143,173,174,198,199,215,216,218,219,220,229,230},
  fprintf('batch %d: special, only cells with active gDataRaw\n',...
          batchid);
  cellargs=cat(2,cellargs,{'behavior'},{1});
  cellargs=cat(2,cellargs,{'speedgt'},{0.5});
  
end

[cellfiledata,cellids,cellfileids]=dbgetscellfile(cellargs{:});

switch batchid,
 case 230,
  % exclude Bom data
  keepcells=ones(size(cellids));
  
  for ii=1:length(keepcells),
     if strcmp(cellids{ii}(1:2),'b0'),
        keepcells(ii)=0;
     end
  end
  keepidx=find(keepcells);
  cellids={cellids{keepidx}};
  cellfileids=cellfileids(keepidx);
  cellfiledata=cellfiledata(keepidx);
  
 case {206,207,221,222,223,224,225,226,227,228,231,232},
   fprintf('batch %d: special, only one cell per site\n',...
           batchid);
   masterid=cat(1,cellfiledata.masterid);
   singleid=cat(1,cellfiledata.singleid);
   try,
       [masterlist,uidx]=unique(masterid,'first');
   catch
       [masterlist,uidx]=unique(masterid);  % older versions of matlab
   end
   singlelist=singleid(uidx);
   keepidx=find(ismember(singleid,singlelist));
   
   cellids=unique({cellfiledata(keepidx).cellid});
   cellfileids=cellfileids(keepidx);
   
   cellfiledata=cellfiledata(keepidx);
end

if isfield(params,'specialbatch'),
    switch lower(params.specialbatch),
      case 'voc',
        
        % spn / center-surround
        keepfiles=zeros(size(cellfiledata));
        cleanfiles=zeros(size(cellfiledata));
        keepcellids={};
        for ii=1:length(cellfiledata),
            parms=dbReadData(cellfiledata(ii).rawid);
            if isfield(parms,'Ref_Subsets') && parms.Ref_Subsets~=5,
                keepfiles(ii)=1;
                keepcellids=union(keepcellids,cellfiledata(ii).cellid);
            end
        end
        
        for ii=1:length(cellfiledata),
            parms=dbReadData(cellfiledata(ii).rawid);
            if ismember(cellfiledata(ii).cellid,keepcellids) &&...
                    ((isfield(parms,'Ref_SNR') && parms.Ref_SNR>=100) ...
                     || ~isfield(parms,'Ref_SNR')),
                cleanfiles(ii)=1;
            end
        end
        
        keepidx=find(cleanfiles);
        cellfileids=cellfileids(keepidx);
        cellfiledata=cellfiledata(keepidx);
        cellids=keepcellids;
        
      case 'cs',
        
        % spn / center-surround
        keepfiles=zeros(size(cellfiledata));
        keepcellids={};
        for ii=1:length(cellfiledata),
            parms=dbReadData(cellfiledata(ii).rawid);
            if (~isfield(parms,'Ref_SplitChannels') ||...
                    strcmpi(strtrim(parms.Ref_SplitChannels),'No')) &&...
                    length(parms.Ref_LowFreq)>1 &&...
                    ismember(parms.Ref_LowFreq(end),[125 250 500 1000]) &&...
                    ismember(parms.Ref_HighFreq(end),[4000 8000 16000 32000]),
                keepfiles(ii)=1;
                keepcellids=union(keepcellids,cellfiledata(ii).cellid);
            end
        end
        
        keepidx=find(keepfiles);
        cellfileids=cellfileids(keepidx);
        cellfiledata=cellfiledata(keepidx);
        cellids=keepcellids;
       case 'lr',
        
        % spn / left-right same spectral features batch
        keepfiles=zeros(size(cellfiledata));
        keepcellids={};
        for ii=1:length(cellfiledata),
            parms=dbReadData(cellfiledata(ii).rawid);
            if isfield(parms,'Ref_SplitChannels') &&...
                    strcmpi(strtrim(parms.Ref_SplitChannels),'Yes') &&...
                    length(parms.Ref_LowFreq)>1 &&...
                    diff(parms.Ref_LowFreq)==0 &&...
                    diff(parms.Ref_HighFreq)==0 && ...
                    ~strcmpi(cellfiledata(ii).behavior,'active') ,
                keepfiles(ii)=1;
                keepcellids=union(keepcellids,cellfiledata(ii).cellid);
            end
        end
        keepidx=find(keepfiles);
        cellfileids=cellfileids(keepidx);
        cellfiledata=cellfiledata(keepidx);
        cellids=keepcellids;
      case 'lr-behavior',
        
        % spn / left-right same spectral features batch
        keepfiles=zeros(size(cellfiledata));
        keepcellids={};
        for ii=1:length(cellfiledata),
            parms=dbReadData(cellfiledata(ii).rawid);
            if isfield(parms,'Ref_SplitChannels') &&...
                    strcmpi(strtrim(parms.Ref_SplitChannels),'Yes') &&...
                    length(parms.Ref_LowFreq)>1 &&...
                    diff(parms.Ref_LowFreq)==0 &&...
                    diff(parms.Ref_HighFreq)==0 && ...
                    strcmpi(cellfiledata(ii).behavior,'active') &&...
                    cellfiledata(ii).isolation>params.miniso,
                keepfiles(ii)=1;
                keepcellids=union(keepcellids,cellfiledata(ii).cellid);
            end
        end
        keepidx=find(keepfiles);
        cellfileids=cellfileids(keepidx);
        cellfiledata=cellfiledata(keepidx);
        cellids=keepcellids;
      case 'lrhl-behavior',
        
        % spn / left-right same spectral features batch
        keepfiles=zeros(size(cellfiledata));
        keepcellids={};
        for ii=1:length(cellfiledata),
            parms=dbReadData(cellfiledata(ii).rawid);
            if isfield(parms,'Ref_SplitChannels') &&...
                    strcmpi(strtrim(parms.Ref_SplitChannels),'Yes') &&...
                    length(parms.Ref_LowFreq)>1 &&...
                    (diff(parms.Ref_LowFreq)~=0 ||...
                     diff(parms.Ref_HighFreq)~=0) && ...
                    strcmpi(cellfiledata(ii).behavior,'active') &&...
                    cellfiledata(ii).isolation>params.miniso,
                keepfiles(ii)=1;
                keepcellids=union(keepcellids,cellfiledata(ii).cellid);
            end
        end
        keepidx=find(keepfiles);
        cellfileids=cellfileids(keepidx);
        cellfiledata=cellfiledata(keepidx);
        cellids=keepcellids;
      case 'hl-behavior',
        
        % spn / left-right same spectral features batch
        keepfiles=zeros(size(cellfiledata));
        keepcellids={};
        for ii=1:length(cellfiledata),
            parms=dbReadData(cellfiledata(ii).rawid);
            if (~isfield(parms,'Ref_SplitChannels') ||...
                strcmpi(strtrim(parms.Ref_SplitChannels),'No')) &&...
                    length(parms.Ref_LowFreq)>1 &&...
                    (diff(parms.Ref_LowFreq)~=0 ||...
                     diff(parms.Ref_HighFreq)~=0) && ...
                    strcmpi(cellfiledata(ii).behavior,'active') &&...
                    cellfiledata(ii).isolation>params.miniso,
                keepfiles(ii)=1;
                keepcellids=union(keepcellids,cellfiledata(ii).cellid);
            end
        end
        keepidx=find(keepfiles);
        cellfileids=cellfileids(keepidx);
        cellfiledata=cellfiledata(keepidx);
        cellids=keepcellids;
      case 'lev-behavior',
        
        % spn / left-right same spectral features batch
        keepfiles=zeros(size(cellfiledata));
        keepcellids={};
        for ii=1:length(cellfiledata),
            parms=dbReadData(cellfiledata(ii).rawid);
            if strcmpi(parms.ReferenceClass,'SpNoise') &&...
                    strcmpi(parms.TrialObjectClass,'MultiRefTar') &&...
                    length(parms.Ref_LowFreq)==1 &&...
                    length(parms.Trial_RelativeTarRefdB)>1 &&...
                    strcmpi(cellfiledata(ii).behavior,'active') &&...
                    cellfiledata(ii).isolation>params.miniso,
                keepfiles(ii)=1;
                keepcellids=union(keepcellids,cellfiledata(ii).cellid);
            end
        end
        keepidx=find(keepfiles);
        cellfileids=cellfileids(keepidx);
        cellfiledata=cellfiledata(keepidx);
        cellids=keepcellids;
      case 'behavior',
        % spn / left-right same spectral features batch
        keepfiles=zeros(size(cellfiledata));
        keepcellids={};
        for ii=1:length(cellfiledata),
            if strcmpi(cellfiledata(ii).behavior,'active'),
                keepfiles(ii)=1;
                keepcellids=union(keepcellids,cellfiledata(ii).cellid);
            end
        end
        keepidx=find(keepfiles);
        cellfileids=cellfileids(keepidx);
        cellfiledata=cellfiledata(keepidx);
        cellids=keepcellids;
    end
    
    if params.activeonly,
        % require at least two active files if activeonly==1
        uniquecellids=cellids;
        allcellids={cellfiledata.cellid};
        
        keepcells=zeros(size(uniquecellids));
        keepfiles=zeros(size(cellfiledata));
        for ii=1:length(uniquecellids),
            ff=find(strcmp(uniquecellids{ii},allcellids));
            if length(ff)>1,
                keepfiles(ff)=1;
                keepcells(ii)=1;
            end
        end
        cellids=cellids(find(keepcells));
        cellfiledata=cellfiledata(find(keepfiles));
        cellfileids=cellfileids(find(keepfiles));
    end
end


if isempty(cellfiledata),
   %keyboard
end
