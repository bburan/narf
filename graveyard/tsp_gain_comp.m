function [setidx,linewidths]=tsp_gain_comp(cellid,ha)
% function tsp_gain_comp(cellid,ha[default=new figure])
%
% created SVD 2013-01 - compare differences in response for each
% tsp/spn file
%
baphy_set_path;
dbopen(1);
cellfiledata=dbgetscellfile('cellid',cellid,'runclassid',[42 8 103],...
                            'ReferenceClass','SpNoise');
filecount=length(cellfiledata);

rasterfs=100;
psthfs=20;
behavior_only=1;

runclassid=cat(2,cellfiledata.runclassid);
stimspeedid=cat(2,cellfiledata.stimspeedid);
activefile=zeros(1,length(cellfiledata));
stimprint=zeros(filecount,6);
stimrelatten=zeros(filecount,2);
stimlow=zeros(filecount,2);
stimhigh=zeros(filecount,2);
prestimsilence=zeros(filecount,1);
duration=zeros(filecount,1);
poststimsilence=zeros(filecount,1);
tarband=zeros(filecount,1);
splitchannels=zeros(filecount,1);
tarstr=cell(filecount,1);
for ii=1:filecount,
   [parms,perf]=dbReadData(cellfiledata(ii).rawid);
   
   if strcmpi(cellfiledata(ii).behavior,'active'),
      activefile(ii)=1;
   end
   if isfield(parms,'Ref_SplitChannels'),
       splitchannels(ii)=strcmpi(parms.Ref_SplitChannels,'Yes');
   end
   if isfield(parms,'Trial_TargetIdxFreq'),
       tif=parms.Trial_TargetIdxFreq;
       tarband(ii)=min(find(tif==max(tif)));
       
       if splitchannels(ii) && tarband(ii)==2,
           bc='i';
       else
           bc='c';
       end
       tarstr{ii}=...
           sprintf('%.1f-%s',parms.Tar_Frequencies(tarband(ii))/1000,bc);
   end
   bandcount=min([length(parms.Ref_LowFreq),...
                 length(parms.Ref_HighFreq)]);
   
   stimrelatten(ii,:)=[parms.Ref_RelAttenuatedB...
                       zeros(2-length(parms.Ref_RelAttenuatedB),1)];
   stimlow(ii,:)=[parms.Ref_LowFreq zeros(2-bandcount,1)];
   stimhigh(ii,:)=[parms.Ref_HighFreq zeros(2-bandcount,1)];
   
   % "fingerprint" is subset and frequency range.  Need to match
   % between passive and active
   thisprint=[parms.Ref_Subsets stimlow(ii,:) stimhigh(ii,:) splitchannels(ii)];
   stimprint(ii,1:length(thisprint))=thisprint;
   prestimsilence(ii)=parms.Ref_PreStimSilence;
   duration(ii)=parms.Ref_Duration;
   if activefile(ii),
       poststimsilence(ii)=0; % force to zero below
   else
       poststimsilence(ii)=parms.Ref_PostStimSilence;
   end
end

if behavior_only,
    activeidx=max(find(activefile));
    printmatch=double(sum(abs(stimprint-repmat(stimprint(activeidx,:),filecount,1)),2)==0);

    fprintf('found %d files with same parameters as active file %d\n',...
            sum(printmatch),activeidx);
    
    usepassidx=find(printmatch & ~activefile');
    passivemergeidx=filecount+1;
    cellfiledata(passivemergeidx).stimfile='Passive Merge';
    stimhigh(passivemergeidx,:)=stimhigh(usepassidx(1),:);
    stimlow(passivemergeidx,:)=stimlow(usepassidx(1),:);
    stimrelatten(passivemergeidx,:)=stimrelatten(usepassidx(1),:);
    splitchannels(passivemergeidx)=splitchannels(usepassidx(1));
    tarstr{passivemergeidx}=tarstr{usepassidx(1)};
               
    usefileidx=find(printmatch);
    
else
    usefileidx=1:filecount;
    passivemergeidx=0;
end

INCLUDEINCORRECT=0;
r=cell(filecount,1);
p=cell(filecount,1);
s=cell(filecount,1);
stimparam=cell(filecount,1);
for ii=1:length(usefileidx),
    fidx=usefileidx(ii);
    parmfile=[cellfiledata(fidx).stimpath cellfiledata(fidx).stimfile];
    spikefile=[cellfiledata(fidx).path cellfiledata(fidx).respfile];
    
    % Load spike PSTH at high sampling rate
    options=[];
    options.channel=cellfiledata(fidx).channum;
    options.unit=cellfiledata(fidx).unit;
    options.rasterfs=rasterfs;
    options.includeprestim=0; [prestimsilence(fidx) poststimsilence(fidx)];
    options.includeincorrect=INCLUDEINCORRECT;
    options.tag_masks={'Reference'};
    
    [tr,rtags]=loadspikeraster(spikefile,options);
    r{fidx}=tr;

    % load at low sampling rate for PSTH
    options.rasterfs=psthfs;
    [p{fidx},ptags]=loadspikeraster(spikefile,options);
    
    % load stimulus at rate matched to r{fidx}--rasterfs
    [s{fidx},stimparam{fidx}]=...
        loadstimfrombaphy(parmfile,[],[],'envelope',rasterfs,0,0,options.includeprestim);

    if size(r{fidx},3)<size(s{fidx},3),
        for jj=1:length(rtags),
            tr=strsep(rtags{jj},',');
            rtags{jj}=strtrim(tr{2});
        end
        
        ff=find(ismember(stimparam{fidx}.tags,rtags));
        s{fidx}=s{fidx}(:,:,ff);
    end
    
    binspertrial=size(r{fidx},1);
    s{fidx}=s{fidx}(:,1:binspertrial,:);
    s{fidx}=s{fidx}(:,:)';
    
    if passivemergeidx && ismember(fidx,usepassidx),
        if fidx==usepassidx(1) || ...
                size(r{fidx},3)>size(r{passivemergeidx},3) || ...
                size(r{fidx},1)~=size(r{passivemergeidx},1)
            
            s{passivemergeidx}=s{fidx};
            r{passivemergeidx}=r{fidx};
        else
            r{passivemergeidx}=cat(2,r{passivemergeidx},r{fidx});
        end
        if fidx==usepassidx(end),
            r{passivemergeidx}=nanmean(r{passivemergeidx},2);
            r{passivemergeidx}=r{passivemergeidx}(:);
        end
    end
    
    %reshape to make stim and response match
    r{fidx}=nanmean(r{fidx},2);
    r{fidx}=r{fidx}(:);
end

if passivemergeidx,
    usefileidx=[usefileidx;passivemergeidx];
end

tt=dbReadTuning(cellid);
if ~isfield(tt,'bf'), tt.bf=0; end
if ~isfield(tt,'bw'), tt.bw=0; end
if ~isfield(tt,'torbf'), tt.torbf=0; end
cellinfo=sprintf('%s - BF %.0f, BW %.2f',cellid,tt.bf,tt.bw);

outpath='/auto/users/svd/data/spn_mod/';

puse=1:(length(usefileidx)-1);
pmax=size(p{usefileidx(1)},1);
for ii=puse,
    pmax=min(size(p{usefileidx(ii)},1),pmax);
end

pmean=[];%zeros(pmax,size(p{usefileidx(1)},3),length(puse));
pset={};
filenames={};
colorset=[0 0 1; 1 0 0; 0.2 0.2 0.8; 0 1 0; 0.4 0.4 0.6];
linespecset={'-','o-','-','s-','-'};
actidxset=find(activefile);
setidx=zeros(size(puse));
linewidths=ones(size(puse));
for ii=puse,
    fidx=usefileidx(ii);
    pset{ii}=squeeze(nanmean(p{fidx}(1:pmax,:,:),2));
    pmean=cat(2,pmean,p{fidx}(1:pmax,:,:));
    
    band1=sprintf('%.1f-%.1f(%d)-%s',stimlow(fidx,1)./1000,...
                  stimhigh(fidx,1)/1000, stimrelatten(fidx,1),'c');
    if splitchannels(fidx),
        bc='i';
    else
        bc='c';
    end
    band2=sprintf('%.1f-%.1f(%d)-%s',stimlow(fidx,2)/1000,...
                  stimhigh(fidx,2)/1000,stimrelatten(fidx,2),bc);
    
    if ~isempty(tarstr{fidx}) && activefile(fidx),
        filenames{ii}=sprintf('%s(%.0f): %s',cellfiledata(fidx).stimfile(8:(end-4)),...
                              cellfiledata(fidx).isolation,tarstr{fidx});
        if tarstr{fidx}(end)=='c',
            setidx(ii)=2;
        else
            setidx(ii)=4;
        end
        if (stimlow(fidx,tarband(fidx))<tt.bf && ...
                stimhigh(fidx,tarband(fidx))>tt.bf) ||...
           (stimlow(fidx,tarband(fidx))<tt.torbf && ...
                stimhigh(fidx,tarband(fidx))>tt.torbf),
            linewidths(ii)=2;
        end
    else
        filenames{ii}=sprintf('%s (%.0f)',cellfiledata(fidx).stimfile(8:(end-4)),...
                              cellfiledata(fidx).isolation);
        if fidx<min(actidxset),
            setidx(ii)=1;
        elseif fidx>max(actidxset),
            setidx(ii)=5;
        else
            setidx(ii)=3;
        end
    end
    
end
pmean=squeeze(nanmean(pmean,2));


[ps,psi]=sort(pmean(:));
bincount=14;
p0=zeros(bincount,1);
pcurve=zeros(bincount,length(pset));
pse=zeros(bincount,length(pset));
stepsize=length(psi)./bincount;
for bb=1:bincount,
    ii=psi(round((bb-1)*stepsize+1):round(bb*stepsize));
    p0(bb)=nanmean(pmean(ii)) * psthfs;
    for jj=1:length(pset),
        pcurve(bb,jj)=nanmean(pset{jj}(ii)) * psthfs;
        pse(bb,jj)=nanstd(pset{jj}(ii))./sqrt(sum(~isnan(pset{jj}(ii)))) ...
            * psthfs;
    end
end

if ~exist('ha','var') || isempty(ha);
    figure;
    ha=gca;
end

cla
ht=zeros(size(puse));
for ii=puse,
    ht(ii)=errorbar(p0,pcurve(:,ii),pse(:,ii),linespecset{setidx(ii)},'Color', ...
                colorset(setidx(ii),:),'LineWidth',linewidths(ii));
    hold on
end
%axis equal

hl=legend(ht,filenames,'Location','northwest','FontSize',6);
set(hl,'interpreter','none');
title([cellinfo ' : ' band1 ' / ' band2]);
xlabel('Average response (all files)');
ylabel('File-specific response');
%print('-f2','-djpeg','-r150',[outpath 'spn_rawcomp.' cellid '.jpg']);

