% function tunedata=kern2tune(strf)
%
% tundata.something(..., signidx)  where signidx=1 is fit to peak
% of step response and signidx=2 is fit to last frame of step
%
% created 2003 SVD
% tried to document SVD 8/23/04 
%
function tunedata=kern2tune(strf)

h=strf.h;
spacebincount=size(h,1);
tbincount=size(h,2);
obincount=15;
sfbincount=8;
skipfit=1;
brotate90=0;

if strcmp(strf.parms.kernfmt,'fft') | strcmp(strf.parms.kernfmt,'sfft') | ...
      strcmp(strf.parms.kernfmt,'pfft') | ...
      findstr(strf.parms.kernfmt,'pfft+') | ...
      strcmp(strf.parms.kernfmt,'lin2'),
   
   if strcmp(strf.parms.kernfmt,'fft') | strcmp(strf.parms.kernfmt,'lin2'),
      phasecount=4;
   elseif length(strf.parms.kernfmt)>4 & findstr(strf.parms.kernfmt,'pfft+'),
      phasecount=str2num(strf.parms.kernfmt(end));
   else
      phasecount=1;
   end
   
   chancount=spacebincount/phasecount;
   Xmax=sqrt(chancount*2);
   [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);
   
   tsf=zeros(Xmax*Xmax,tbincount);
   tsf(cfilt,:,:)=squeeze(sum(reshape(h,chancount,phasecount,...
                                      tbincount),2));
   tsf(cfiltconj,:,:)=tsf(cfilt,:,:);
   tsf=reshape(tsf,Xmax,Xmax,tbincount);
   
   phasecount=1; % fix to 1 since we summed over phase
elseif strcmp(strf.parms.kernfmt,'fft') | strcmp(strf.parms.kernfmt,'pfft'),
   if strcmp(strf.parms.kernfmt,'fft'),
      phasecount=4;
   else
      phasecount=1;
   end
   
   chancount=spacebincount/phasecount;
   Xmax=sqrt(chancount*2);
   [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);
   
   if 0 & kcount==1,
      tsf=zeros(Xmax*Xmax,tbincount,phasecount);
      for ii=1:phasecount,
         tsf(cfilt,:,ii)=h((1:chancount)+(ii-1)*chancount,:);
      end
   else
      tsf=zeros(Xmax*Xmax,tbincount);
      tsf(cfilt,:,:)=squeeze(sum(reshape(h,chancount,phasecount,...
                                         tbincount),2));
   end
   tsf(cfiltconj,:,:)=tsf(cfilt,:,:);
   tsf=reshape(tsf,Xmax,Xmax,tbincount);
   
else
   % deal with other kernfmts!
   disp('kern2tune.m requires strfs with fft or pfft kernfmt! sorry');
   return
end



Xmax=size(tsf,1);
Ymax=size(tsf,2);
tbinuse=min([tbincount 12]);

signcount=2;
parmcount=4;
tparmcount=6;
tunedata.orsf=zeros(obincount,sfbincount);
tunedata.sftime=zeros(sfbincount,tbinuse,signcount);
tunedata.ortime=zeros(obincount,tbinuse,signcount);
tunedata.por=zeros(obincount,signcount);
tunedata.psf=zeros(sfbincount,signcount);
tunedata.ptime=zeros(tbinuse,signcount);
tunedata.pph=zeros(phasecount,signcount);
tunedata.fitor=zeros(parmcount,signcount);
tunedata.fitsf=zeros(parmcount,signcount);
tunedata.fittime=zeros(tparmcount,signcount);

bincounts=[obincount sfbincount tbinuse];

% calculate mean across resamples in GR space
h=mean(sf2gr(tsf(:,:,1:tbinuse,:,:),obincount,sfbincount,0,brotate90),5);

% figure out time of peak positive step response, integrated over space
th=cumsum(h,3);
ttime=squeeze(mean(mean(mean(th,1),2),4));

ttmax=min(find(ttime==max(ttime)));
if ttmax==1,
   ttmax=min(find(abs(ttime)==max(abs(ttime))));
end

orsf=mean(th(:,:,ttmax,:),4);
orsf=mean(th(:,:,ttmax,:).*abs(th(:,:,ttmax,:)),4);
xx=-1:1;
g=exp(-xx.^2);
g=g./sum(g(:));
orsf=rconv2(orsf,g);  % reflected boundaries along sf dim
orsf=cconv2(orsf,g'); % cicular boundaries along sf dim
smaxidx=min(find(orsf(:,:,1)==max(max(orsf(:,1:(sfbincount-1),1)))));
[oo,ss]=ind2sub([obincount,sfbincount],smaxidx);
fprintf('[oo,ss]=[%d,%d]\n',oo,ss);

binslices=[oo,ss,ttmax];

[tunedata.por,tunedata.psf,tunedata.ptime,tunedata.pph,...
 tunedata.orsf,tunedata.ortime,tunedata.sftime,phtime]=...
   kern2tunecurves(tsf,brotate90,bincounts,binslices);


%
% separability...
%

% separability stats for each resampled est

%h=sf2gr(tsf(:,:,:,:,respidx),obincount,sfbincount,0);
%h=h(:,:,1:tbinuse,:);

if tbincount>1,
   % alpha: space time separability for full kernel
   t=permute(tsf(:,:,1:tbinuse,:),[1 2 4 3]);
   t=reshape(t,Xmax*Ymax*phasecount,tbinuse);
   [u,s,v]=svd(t);
   sd=diag(s).^2;
   tunedata.tssep=sd(1)./sum(sd);
   
   % alpha_p: mean space time sep within each phase channel
   t=permute(h,[1 2 4 3]);
   tsphsep=zeros(phasecount,1);
   for phaseidx=1:phasecount,
      t=reshape(h(:,:,:,phaseidx),obincount*sfbincount,tbinuse);
      [u,s,v]=svd(t);
      sd=diag(s).^2;
      tsphsep(phaseidx)=sd(1)./sum(sd);
   end
   tunedata.tsphmsep=mean(tsphsep);
   
else
   tunedata.tssep=1;
   tunedata.tsphmsep=ones(phasecount,1);
end


t=repmat(mean(h,4),[1 1 1 phasecount]);
tunedata.phmsep=1-var(t(:)-h(:))./var(h(:));


tunedata.obins=linspace(0,180,obincount+1);
tunedata.obins=tunedata.obins(1:end-1);

tunedata.sfbins=linspace(1,round((Xmax+1)/2),sfbincount+1);
%sfbins=linspace(1,round((Xmax-1)/2),sfbincount+1);
tunedata.sfbins=tunedata.sfbins(1:end-1);

% use center of time bins for fit
tunedata.tbins=((0:(tbinuse-1))+0.5)*strf.parms.tbinms;


% 
% fit explicit tuning models
%

disp('only fitting to signidx=1');
for signidx=1:1,
   tunedata.fitor(:,signidx)=...
       fitflatcos(tunedata.obins,tunedata.por(:,signidx));
   tunedata.fitsf(:,signidx)=...
       fitgauss1d(log2(tunedata.sfbins),tunedata.psf(:,signidx));
   if tbincount>1,
      tunedata.fittime(:,signidx)=...
          fitexpdecaydiff(tunedata.tbins,tunedata.ptime(:,signidx));
   end
end
return

figure(1)
clf

subplot(2,2,1);
errorbar(obins,mean(tunedata.por(:,1,:),3),...
         std(tunedata.por(:,1,:),0,3));
title([strf.name,' orientation']);
subplot(2,2,2);
errorbar(sfbins,mean(tunedata.psf(:,1,:),3),...
         std(tunedata.psf(:,1,:),0,3));
title([strf.name,' spatial freq']);
subplot(2,2,3);
errorbar(tbins,mean(tunedata.ptime(:,1,:),3),...
         std(tunedata.ptime(:,1,:),0,3));
title([strf.name,' time']);
subplot(2,2,4);
%morsf=max(abs(tunedata.orsf(:)));
%imagesc(obins,sfbins,tunedata.orsf,[-morsf morsf]);
imagesc(obins,sfbins,tunedata.orsf);
colormap(hot);

drawnow
