% dmsscratch.m
%
% created SVD 7/25/03
%

dbopen;
BATQUEUEID=[];

cellid='v0063';

clear params

FULLTIME=0;

params.outfile=['/auto/k5/david/data/dms/',cellid,'.dms6.res.mat'];
params.zipoutfile=1;
params.cellid=cellid;
params.resploadcmd='resploaddms';
if FULLTIME,
   params.resploadparms={'psth',1,0};
   params.respfiltercmd='';
   params.respfilterparms={};
else
   params.resploadparms={'psth',1,1,{'AB','A','B'}};
   params.respfiltercmd='respvarbins';
   params.respfilterparms={3,15,0};
end

params.stimloadcmd='loaddmsstim';
scalepix=32;
if 0,
   params.stimloadparms={scalepix,8,16};
   params.stimfiltercmd='movphasesep';
   params.stimfilterparms={0,0,0,1,0};
   params.kernfmt='fft';
elseif 0,
   params.stimloadparms={16,0,16};
   params.stimfiltercmd='movlinpower';
   params.stimfilterparms={0,0,1,1,0};
   params.kernfmt='lin2';
elseif 1,
   params.stimloadparms={16,0,16};
   params.stimfiltercmd='movpower';
   params.stimfilterparms={0,0,1,1,0};
   params.kernfmt='pfft';
elseif 0,
   params.stimloadparms={scalepix,0,32};
   params.stimfiltercmd='mov2wavsparse';
   params.stimfilterparms={5,4,4};
   params.kernfmt='wav';
elseif 0
   addpath /auto/k1/jess/ssa/
   params.stimloadparms={scalepix,scalepix+1,16};
   params.stimfiltercmd='sca';
   params.stimfilterparms={};
   params.kernfmt='space';
else
   params.stimloadparms={scalepix,scalepix+1,16};
   params.stimfiltercmd='';
   params.stimfilterparms={};
   params.kernfmt='space';
end
params.resampcount=20;
params.resampfmt=1;
params.sfsstep=5;
params.sfscount=10;
params.sffiltsigma=5;
params.fitfrac=0.075;
params.predfrac=0.075;
params.fitboot=0;
params.docellfit2=0;
params.repexclude=0;
params.nloutparm=4;
params.decorrspace=2;
params.decorrtime=0;
params.meansub=1;
params.noisecount=200;

if FULLTIME,
   params.maxlag=[-6 15];
   params.smoothtime=2;
else
   params.maxlag=[-4 4];
   params.smoothtime=0;
end

cellfiledata=dbgetscellfile('cellid',cellid,'runclassid',10);

starttimes=[];
stoptimes=[];
attcodes={'A+','B+','a-b-','A+B+'};
fcount=0;
for fidx=1:length(cellfiledata),
   params.respfiles{fidx}=[cellfiledata(fidx).path,...
                    cellfiledata(fidx).respfile];
   params.stimfiles{fidx}=[cellfiledata(fidx).stimpath,...
                    cellfiledata(fidx).stimfile];
   
   fprintf('%d. %s\n',fidx,params.stimfiles{fidx});
   
   starttimes=cat(1,starttimes,0);
   stoptimes=cat(1,stoptimes,0);
   %starttimes=cat(1,starttimes,1);
   %stoptimes=cat(1,stoptimes,imfileinfo(params.stimfiles{fidx},1));
end

% load info about the targets for reference, this can be found
% in the A+ response file
r=load(params.respfiles{1});
[fmask,crop]=movfmask(size(r.target,1),0,size(r.target,1));
targs=movresize(floor(r.target),params.stimloadparms{1},fmask,crop);
targs=feval(params.stimfiltercmd,targs,params.stimfilterparms{:});

kerncomp4;

disp('just ran kerncomp4.m');

keyboard

attidx=1;
xcloadfiles;

bsSAuse=[1 2];

if FULLTIME,
   

   disp('FULLTIME=1 not VALID!!! NEED TO SWITCH OVER TO DB-DRIVEN STRUCTURE')
   return
   
   rlen=zeros(attcount,1);
   
   % compute corr mtx for each stim set
   stim=cat(1,xcdata.stim);
   bresp=cat(1,xcdata.resp);
   resp=ones(length(bresp),attcount).*nan;
   tt=0;
   for attidx=1:attcount,
      tresp=xcdata(attidx).resp;
      tnan=find(isnan(tresp));
      tnan=tnan(find(tnan<length(tresp)-params.maxlag(2)));
      tresp(tnan+round(params.maxlag(2)./2))=nan;
      tresp(tnan+params.maxlag(2))=nan;
      
      resp(tt+(1:length(xcdata(attidx).resp)),attidx)=tresp;
      tt=tt+length(tresp);
      
      tnan=find(isnan(tresp));
      rlen(attidx)=length(tnan);
      
      if ismember(attidx,bsSAuse),
         bsSA=bsSA+xcdata(attidx).stim'*xcdata(attidx).stim;
         bsSAcount=bsSAcount+rlen(attidx);
      end
   end
else
   MINLEN=max([10 flashtimes]);
   
   attcount=length(attcodes)
   flashcount=size(resp,1);
   flashlen=size(resp,2)/attcount;
   resp=reshape(resp,flashcount,flashlen,attcount);
   
   flashtimes=3:15;
   
   if size(resp,2)>1,
      reresp=squeeze(nanmean(permute(resp(:,flashtimes,1:3),[3 1 2])));
      mresp=nanmean(reresp);
      sresp=nanstd(reresp) ./ sqrt(sum(sum(~isnan(reresp(:,3:4)))));
      meanallresp=nanmean(resp(:));
      stdallresp=nanstd(resp(:)) ./ ...
          sqrt(length(resp(find(~isnan(resp)))));
      
      
      
      savemresp=mresp;
      mresp=mresp-meanallresp;
      if length(find(abs(mresp)./sresp>=4))>0,
         mresp(find(abs(mresp)./sresp<4))=0;
      elseif length(find(abs(mresp)./sresp>=2))>0,
         mresp(find(abs(mresp)./sresp<2))=0;
      end
      mresp=mresp./sqrt(sum(mresp.^2));
      
      figure(1);
      clf
      errorbar((flashtimes-1).*16,savemresp,sresp.*4);
      hold on
      nzidx=find(abs(mresp)>0);
      scatter((flashtimes(nzidx)-1).*16,savemresp(nzidx),'ro','filled');
      plot([0 max(flashtimes)*16-16],[1 1].*meanallresp,'k--');
      hold off
      
      disp('plotted mresp');
      
      fresp=zeros(flashcount,attcount);
      for attidx=1:attcount
         fresp(:,attidx)=resp(:,flashtimes,attidx)*mresp';
      end
   else
      fresp=squeeze(resp);
   end
   
   rlen=sum(~isnan(fresp),1);
   
   tbsa=find(sum(~isnan(fresp(:,bsSAuse)),2)>0);
   bsSA=stim(tbsa,:)'*stim(tbsa,:)./length(tbsa);
   
end

resp=fresp;

rsize=size(resp);
spacecount=size(stim,2);
respcount=attcount;

firstseg=1;

xccore;


%bsSA=mean(sSA2,4);
%bsSA=bsSA.*repmat(reshape(rlen,1,1,attcount),[spacecount spacecount 1]);
%bsSA=sum(bsSA,3)./sum(rlen);

[ua,sa,va]=svd(bsSA);

% ssa=u*s*u'
% ssa_pc=ua'* (u * s * u') * ua
% ssa^-1=u*s^-1*u'
% ssa_pc^-1 = ua'* u* s^-1 * u'* ua

He=zeros(size(SR));
for attidx=1:attcount,
   for resampidx=1:params.resampcount,
      [u,s,v]=svd(sSA2(:,:,attidx,resampidx));
      ssa_pc1 = ua'* u* s^-1 * u';
         
      if FULLTIME,
         He(:,:,attidx,resampidx)=ssa_pc1 * H(:,:,1,attidx,resampidx);
      else
         %ssa_pc1 = s^-1 * ua';
         He(:,:,attidx,resampidx)=ssa_pc1 * SR(:,:,attidx,resampidx);
      end
   end
end

mHe=mean(He,4);
eHe=std(He,1,4).*sqrt(params.resampcount-1);
if FULLTIME
   postlats=-params.maxlag(1)+1:size(mHe,2);
   maxframe=postlats(5);
else
   postlats=1:size(mHe,2);
   maxframe=mean(postlats);
end

erange=1:20;

snr=abs(mHe)./eHe;
snr((max(erange)+1):end,:,:)=0;

figure;
amax=max(abs(mHe(:)).* (abs(snr(:))>1));
amax2=max(max(max((abs(mHe(erange,:,:))+eHe(erange,:,:)).* ...
                  (abs(snr(erange,:,:))>1))));

latrange=-1:1;
latcount=length(latrange);
for attidx=1:attcount,
   for latidx=1:latcount,
      subplot(attcount,latcount,(attidx-1)*latcount+latidx);
      curframe=maxframe+latrange(latidx);
      tsnr=abs(mHe(erange,curframe,attidx))./eHe(erange,curframe,attidx);
      errorbar(mHe(erange,curframe,attidx),eHe(erange,curframe,attidx));
      hold on
      plot(erange,erange.*0,'k--');
      
      plot(erange,eHe(erange,curframe,attidx).*2,'k:');
      plot(erange,-eHe(erange,curframe,attidx).*2,'k:');
      
      hold off
      title(sprintf('mHe attcode=%s lat=%d',...
                    attcodes{attidx},latrange(latidx)));
      axis([0 max(erange)+1 -amax2 amax2]);
   end
   
   if attidx==attcount,
      xlabel('eigenvector num');
   end
   
   %subplot(attcount,3,attidx*3-2);
   %imagesc(snr(erange,postlats,attidx)'>1);
   %title(sprintf('snr>1 attcode=%s',attcodes{attidx}));
  % 
  % subplot(attcount,3,attidx*3-1);
  % imagesc(mHe(erange,postlats,attidx)' .* ...
  %         (snr(erange,postlats,attidx)'>0.5),[-amax amax]);
  % title(sprintf('mHe attcode=%s',attcodes{attidx}));
   
end

titlecodes=attcodes;
for attidx=1:attcount,
   titlecodes{attidx}=[cellid,' ',titlecodes{attidx}];
end

figure
smd=(1-(snr./1.5).^(-2));
smd=smd.*(smd>0);
smd(find(isnan(smd)))=0;

tkern=mHe(:,postlats,:).*(smd(:,postlats,:));
tkern=ua*tkern(:,:);
tkern=cat(3,reshape(tkern,spacecount,length(postlats),attcount),...
          ones(spacecount,length(postlats),2).*nan);
showkern(tkern,params.kernfmt,params.stimloadparms{3},titlecodes,0,16);

tkern=cat(3,ones([size(targs,1) 1 attcount]).*nan,...
          reshape(sqrt(targs),spacecount,1,size(targs,2)));
showkern(tkern,params.kernfmt,params.stimloadparms{3},...
         {titlecodes{:},'Target A','Target B'},0,16);
xlabel('<-- neg latency  (peak)  pos latency -->');

if 0,
showstim(ua,params.kernfmt,params.stimloadparms{3},5,5,...
         ['cell ',cellid,' PCs']);
end

figure
smd=(1-(abs(mH)./eH./1.0).^(-2));
smd=smd.*(smd>0);
smd(find(isnan(smd)))=0;
%tkern=squeeze(mH(:,maxframe,:,:).*(smd(:,maxframe,:,:)));
tkern=squeeze(mH(:,maxframe,:,:));
tkern=permute(tkern,[1 3 2]);
showkern(tkern,params.kernfmt);

for ii=1:attcount,
   subplot(size(tkern,3),size(tkern,2),ii);
   title(titlecodes{ii});
   subplot(size(tkern,3),size(tkern,2),(size(tkern,3)-1).*size(tkern,2)+ii);
   xlabel('');
end

