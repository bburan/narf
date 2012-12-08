% function [por,psf,ptime,oo,ss,ttmax]=xctuning(strf)
%
% extract tuning properties from strfs.  currently designed to work
% with pfft, sfft strfs
%
% created SVD 4/16/03
%
function [por,psf,ptime,oo,ss,ttmax]=xctuning(strf)

obincount=8;
sfbincount=8;
tbinuse=size(strf(1).h,2);
kerncount=length(strf(:));

kernfmt=strf(1).parms.kernfmt;
if strcmp(kernfmt,'fft') | strcmp(kernfmt,'fftgr'),
   phasecount=4;
elseif strcmp(kernfmt,'pfft') | strcmp(kernfmt,'pfftgr'),
   phasecount=1;
else
   disp('error, doesn''t support anything but sfft and pfft kernfmts');
   return
end

por=zeros(obincount,kerncount);
psf=zeros(sfbincount,kerncount);
pph=zeros(phasecount,kerncount);
ptime=zeros(tbinuse,kerncount);
ortime=zeros(obincount,tbinuse,kerncount);
sftime=zeros(sfbincount,tbinuse,kerncount);
phtime=zeros(phasecount,tbinuse,kerncount);
orsf=zeros(obincount,sfbincount,kerncount);
oo=zeros(kerncount,1);
ss=zeros(kerncount,1);
ttmax=zeros(kerncount,1);

for kernidx=1:kerncount,
   kernfmt=strf(kernidx).parms.kernfmt;
   
   if strcmp(kernfmt,'fft') | strcmp(kernfmt,'fftgr'),
      phasecount=4;
      tsf=squeeze(fft2gr(strf(kernidx).h,obincount,sfbincount,phasecount));
      
      tsf=flipdim(tsf,2);
      tsf=reshape(permute(tsf,[2 4 1 3 5]),sfbincount*phasecount,...
                  obincount,tbincount,kcount);
   elseif strcmp(kernfmt,'pfft') | strcmp(kernfmt,'pfftgr'),
      phasecount=1;
      tsf=squeeze(fft2gr(strf(kernidx).h,obincount,sfbincount,phasecount));
      tsf=flipdim(tsf,2);
      tsf=permute(tsf,[2 1 3 4]);
   else
      disp('error, doesn''t support anything but sfft and pfft kernfmts');
      return
   end   
   
   Xmax=sqrt(size(strf(kernidx).h,1)./phasecount.*2);
   tbincount=size(strf(kernidx).h,2);
   
   obins=linspace(0,180,obincount+1);
   obins=obins(1:end-1);
   sfbins=linspace(1,round((Xmax+1)/2),sfbincount+1);
   sfbins=sfbins(1:end-1);
   tbins=(0:(tbincount-1))*strf(kernidx).parms.tbinms;
   
   h=tsf;
   th=cumsum(h,3);
   
   % figure out time of peak positive step response, integrated over space
   ttime=squeeze(mean(mean(mean(th,1),2),4));
   ptime=squeeze(mean(mean(mean(h,1),2),4));
   
   % restrict ttmax to first 8 bins?
   ttime=ttime(1:8);
   
   if max(ttime)>0,
      ttmax(kernidx)=min(find(ttime==max(ttime)));
   else
      ttmax(kernidx)=min(find(abs(ttime)==max(abs(ttime))));
   end
   torsp=reshape(th(:,:,ttmax(kernidx),:),obincount,sfbincount*phasecount);
   tsfsp=permute(th(:,:,ttmax(kernidx),:),[2,1,4,3]);
   
   [u,s,v]=svd(torsp);
   tpor=u(:,1);
   tmor=mean(torsp,2);
   if sum(tpor.*tmor)<0,
      tpor=-tpor;
   end
   por(:,kernidx)=tpor;
      
   tsfsp=tsfsp(:,:);
   [u,s,v]=svd(tsfsp);
   tpsf=u(:,1);
   tmsf=mean(tsfsp,2);
   if sum(tpsf.*tmsf)<0,
      tpsf=-tpsf;
   end
   psf(:,kernidx)=tpsf;
   
   
end

keyboard


return



titles={strf(:).name};
if isfield(strf(1),'tbinms'),
   tbinms=strf(1).tbinms;
else
   tbinms=16;
end

for kernidx=1:kerncount,
   tsf=ones([size(strf(kernidx).h) kerncount])*nan;
   tsf(:,:,kernidx)=strf(kernidx).h;
   
   showkern(tsf,strf(kernidx).parms.kernfmt,...
            strf(kernidx).parms.iconside,titles,1,tbinms);
   
end
