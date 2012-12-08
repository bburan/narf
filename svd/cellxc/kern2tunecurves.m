% function [por,psf,ptime,pph,orsf,ortime,sftime,phtime,oo,ss,ttmax]=...
%     kern2tunecurves(tsfIR,brotate90,bincounts,binslices)
%
% given an STRF (sx X sy X t X phase), calculate orienation, sf,
% phase, temporal tuning curves ('p*') and fit or, sf, time ('fit*')
% 
% Created SVD 9/1/01
% 
function [por,psf,ptime,pph,orsf,ortime,sftime,phtime,oo,ss,ttmax]=...
    kern2tunecurves(tsfIR,brotate90,bincounts,binslices)

if ~exist('brotate90','var'),
   brotate90=0;
end
if ~exist('bincounts','var'),
   obincount=15;
   sfbincount=8;
   tbinuse=12;
else
   obincount=bincounts(1);
   sfbincount=bincounts(2);
   tbinuse=bincounts(3);
end

Xmax=size(tsfIR,1);
Ymax=size(tsfIR,2);
tbincount=size(tsfIR,3);
phasecount=size(tsfIR,4);

obins=linspace(0,180,obincount+1);
obins=obins(1:end-1);
sfbins=linspace(1,round((Xmax+1)/2),sfbincount+1);
sfbins=sfbins(1:end-1);
tbins=(0:(tbinuse-1))*14;
signcount=2;

por=zeros(obincount,signcount);
psf=zeros(sfbincount,signcount);
pph=zeros(phasecount,signcount);
ptime=zeros(tbinuse,signcount);
ortime=zeros(obincount,tbinuse,signcount);
sftime=zeros(sfbincount,tbinuse,signcount);
phtime=zeros(phasecount,tbinuse,signcount);
orsf=zeros(obincount,sfbincount,signcount);

h=sf2gr(tsfIR(:,:,1:tbinuse,:),obincount,sfbincount,0,brotate90);
th=cumsum(h,3);

tsf=cumsum(tsfIR(:,:,1:tbinuse,:),3);
tsfIRsum=squeeze(sum(tsfIR(:,:,1:tbinuse,:),4));

hsum=sf2gr(tsfIRsum,obincount,sfbincount,0,brotate90);
thsum=cumsum(hsum,3);

% figure out time of peak positive step response, integrated over space
ttime=squeeze(mean(mean(mean(th,1),2),4));
ttime2=squeeze(mean(mean(mean(h,1),2),4));

% restrict ttmax to first 8 bins?
if length(ttime)>8,
   ttime=ttime(1:8);
end

if max(ttime)>0,
   ttmax=min(find(ttime==max(ttime)));
else
   ttmax=min(find(abs(ttime)==max(abs(ttime))));
end

if max(ttime2)>0,
   ttmax2=min(find(ttime2==max(ttime2)));
else
   ttmax2=min(find(abs(ttime2)==max(abs(ttime2))));
end

% generate small kernel for smoothing
xx=-1:1;
g=exp(-(xx*1.5).^2);
g=g./sum(g(:));

% find wx,wy coordinates for phase tuning
%ttsf=tsfIRsum(:,:,ttmax2);
ttsf=sum(tsf(:,:,1:ttmax),3);
%ttsf=conv2(ttsf,g,'same');
%ttsf=conv2(ttsf,g');

if max(ttime2)>0,
   phmaxidx=min(find(ttsf==max(ttsf(:))));
else
   phmaxidx=min(find(ttsf==min(ttsf(:))));
end
[phwx,phwy]=ind2sub([Xmax,Ymax],phmaxidx);

%
% orsf signature is calculated at peak time in the IR
orsf(:,:,1)=hsum(:,:,ttmax,:);
%orsf(:,:,1)=mean(th(:,:,ttmax,:),4);
%orsf(:,:,1)=mean(th(:,:,ttmax,:).*abs(th(:,:,ttmax,:)),4);

% smooth orsf signature and find peak:
orsf(:,:,1)=conv2(orsf(:,:,1),g,'same');  % reflected boundaries along sf dim
orsf(:,:,1)=cconv2(orsf(:,:,1),g'); % cicular boundaries along or. dim

if max(ttime)>0,
   smaxidx=min(find(orsf(:,:,1)==max(max(orsf(:,1:(sfbincount-1),1)))));
else
   smaxidx=min(find(orsf(:,:,1)==min(min(orsf(:,1:(sfbincount-1),1)))));
end
[oo,ss]=ind2sub([obincount,sfbincount],smaxidx);

if exist('binslices'),
   if length(binslices)>0,
      oo=binslices(1);
   end
   if length(binslices)>1,
      ss=binslices(2);
   end
   if length(binslices)>2,
      ttmax=binslices(3);
   end
end
fprintf('[oo,ss,ttmax]=[%d,%d,%d] (k2tc.m)\n',oo,ss,ttmax);

if oo<obincount,
   if oo>1,
      oo2=(oo-1):(oo+1);
   else
      oo2=[oo oo+1 obincount];
   end
else
   oo2=[1 oo-1 oo];
end
if sum(mean(mean(th(oo2,:,ttmax,:),1),4),2)==0,
   oo2=1:obincount;
end

% 1/26/02 SVD : use only nearest bins in sf for or/temp tuning
% rather than all sf bins
if 1,
   if ss<sfbincount,
      if ss>1,
         ss2=(ss-1):(ss+1);
      else
         ss2=[ss ss+1];
      end
   else
      ss2=[ss-1 ss];
   end
else
   ss2=1:(sfbincount-1);
end

%th=th.*abs(th);
% old way... only slice over max or/sf -- pieces
por(:,1)=psqrt(squeeze(mean(mean(th(:,ss,ttmax,:),2),4)));
if max(abs(por(:,1)))>0,
   por(:,1)=por(:,1)./max(abs(por(:,1)));
end
psf(:,1)=psqrt(squeeze(mean(mean(th(oo,:,ttmax,:),1),4)'));
if max(abs(psf(:,1)))>0,
   psf(:,1)=psf(:,1)./max(abs(psf(:,1)));
end
%por(:,1)=psqrt(squeeze(mean(mean(h(:,ss2,ttmax2,:),2),4)));
%psf(:,1)=psqrt(squeeze(mean(mean(h(oo2,:,ttmax2,:),1),4)'));
por(:,2)=psqrt(squeeze(mean(mean(th(:,ss2,end,:),2),4)));
por(:,2)=por(:,2)./sqrt(sum(por(:,2).^2));
psf(:,2)=psqrt(squeeze(mean(mean(th(oo,:,end,:),1),4))');
psf(:,2)=psf(:,2)./sqrt(sum(psf(:,2).^2));
%pph(:,2)=psqrt(squeeze(mean(mean(th(oo,ss,end,:),1),2)));

%pph(:,1)=psqrt(squeeze(mean(mean(h(oo,ss,ttmax2,:),1),2)));
%pph(:,1)=squeeze(tsfIR(phwx,phwy,ttmax2,:));
pph(:,1)=squeeze(tsf(phwx,phwy,ttmax,:));
pph(:,2)=squeeze(tsf(phwx,phwy,end,:));

% less old way... only slice over max or/sf -- pieces
%pph(:,1)=psqrt(squeeze(mean(mean(th(oo,ss,ttmax,:),1),2)));

% new way... take first PC to get sf and or tuning
for signidx=0:-1, % 1:signcount,
   if signidx==1,
      torsp=reshape(th(:,:,ttmax,:),obincount,sfbincount*phasecount);
      tsfsp=permute(th(:,:,ttmax,:),[2,1,4,3]);
   else
      torsp=reshape(th(:,:,end,:),obincount,sfbincount*phasecount);
      tsfsp=permute(th(:,:,end,:),[2,1,4,3]);
   end
   
   
   [u,s,v]=svd(torsp);
   tpor=u(:,1);
   tmor=mean(torsp,2);
   if sum(tpor.*tmor)<0,
      tpor=-tpor;
   end
   por(:,signidx)=tpor;
   
   tsfsp=tsfsp(:,:);
   [u,s,v]=svd(tsfsp);
   tpsf=u(:,1);
   tmsf=mean(tsfsp,2);
   if sum(tpsf.*tmsf)<0,
      tpsf=-tpsf;
   end
   psf(:,signidx)=tpsf;
end

ptime(:,1)=psqrt(squeeze(mean(mean(mean(th(oo2,ss2,:,:),1),2),4)));
if max(abs(ptime(:,1)))>0,
   ptime(:,1)=ptime(:,1)./sqrt(sum(ptime(:,1).^2));
end
ortime(:,:,1)=psqrt(squeeze(mean(mean(th(:,ss,:,:),2),4)));
sftime(:,:,1)=psqrt(squeeze(mean(mean(th(oo2,:,:,:),1),4)));
phtime(:,:,1)=psqrt(squeeze(mean(mean(th(oo,ss,:,:),1),2)))';

% skip, already calculated above
%orsf(:,:,1)=mean(th(:,:,ttmax,:),4);
orsf(:,:,2)=psqrt(mean(th(:,:,end,:),4));

ptime(:,2)=psqrt(squeeze(mean(mean(mean(th,1),2),4)));
sftime(:,:,2)=psqrt(squeeze(mean(mean(th,1),4)));
ortime(:,:,2)=psqrt(squeeze(mean(mean(th,2),4)));
phtime(:,:,2)=psqrt(squeeze(mean(mean(th,1),2)))';


return


% function x=psqrt(x2)
% return sqrt of abs(x2) .* sign(x2)
function x=psqrt(x2)

%disable, since not using squared th
x=x2;
return

x=sqrt(abs(x2)).*sign(x2);
