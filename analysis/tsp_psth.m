% function tsp_psth(cellid)
%
% example:  tsp_psth('por023a-b1');
%
function tsp_psth(cellid)

dbopen;

detparmfile='';
if ~exist('INCLUDEINCORRECT','var'),
   INCLUDEINCORRECT=0;
end
if ~exist('MUA','var'),
   MUA=0;
end
if ~exist('rasterfs','var') || isempty(rasterfs),
   rasterfs=1000;
end
if ~exist('nanlicks','var'),
   nanlicks=0;
end

torcdata=dbgetscellfile('cellid',cellid,'runclassid',[1]);
cellfiledata=dbgetscellfile('cellid',cellid,'runclassid',[42 8 103]);
filecount=length(cellfiledata);

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
firstactiveidx=min(find(activefile));
if isempty(firstactiveidx),
    error('no active TSP data for this cell');
end
printmatch=double(sum(abs(stimprint-repmat(stimprint(firstactiveidx,:),filecount,1)),2)==0);

useidx=find(printmatch);

r={};
rasterfs=1000;
rtarg={};
t={};
targfs=40;

minactiveidx=min(find(activefile));
preidx=max(find(~activefile(1:minactiveidx) & printmatch(1:minactiveidx)'));
postidx=min(find(~activefile(minactiveidx:end) & printmatch(minactiveidx:end)'))+minactiveidx-1;
if isempty(preidx),
    preidx=min(find(~activefile(1:minactiveidx)));
end
rmax=round(rasterfs.*(prestimsilence(preidx)+duration(preidx)+poststimsilence(preidx)));

mresp=0;
for ii=1:length(useidx),
    parmfile=[cellfiledata(useidx(ii)).stimpath cellfiledata(useidx(ii)).stimfile];
    spikefile=[cellfiledata(useidx(ii)).path cellfiledata(useidx(ii)).respfile];
    
    options=[];
    options.channel=cellfiledata(useidx(ii)).channum;
    options.unit=cellfiledata(useidx(ii)).unit;
    options.rasterfs=rasterfs;
    options.psthfs=25;
    options.includeprestim=[prestimsilence(ii) 0];
    options.includeincorrect=INCLUDEINCORRECT;
    options.tag_masks={'SPECIAL-COLLAPSE-REFERENCE'};
    
    [tr,tags]=loadspikeraster(spikefile,options);
    if size(tr,1)<rmax,
        tr((end+1):rmax,:,:)=nan;
    end
    r{ii}=tr(1:rmax,:,:);
    
    smcount=round(options.rasterfs./options.psthfs);
    smfilt=ones(smcount,1)./smcount.*1000;
    mr=rconv2(squeeze(nanmean(r{ii},2)),smfilt);
    mr=max(mr(:));
    mresp=max(mresp,mr);
end

figure;
if length(torcdata)>0,
    ha1=subplot(4,3,1);
    ha2=subplot(4,3,2);
    parmfile=[torcdata(1).stimpath torcdata(1).stimfile];
    spikefile=[torcdata(1).path torcdata(1).respfile];
    
    tor_tuning(parmfile,spikefile,options.channel,options.unit,[ha1 ha2])
end

for ii=1:length(useidx),
    ha=subplot(4,3,ii+1);
    options.psth=1;
    options.psthmax=mresp.*0.95;
    options.PreStimSilence=prestimsilence(ii);
    parmfile=[cellfiledata(useidx(ii)).stimpath cellfiledata(useidx(ii)).stimfile];
    raster_plot(parmfile,r{ii},tags,ha,options);
    
    ht=title(sprintf('%s - %s',cellid,cellfiledata(useidx(ii)).stimfile),...
          'Interpreter','none');
    if activefile(useidx(ii));
        set(ht,'FontWeight','Bold');
    end
end

colormap default

fprintf('print command:  print -f1 -djpeg -r150 /auto/users/daniela/plots/%s.jpg\n',cellid);



