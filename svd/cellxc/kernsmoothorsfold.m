% function hs=kernsmoothorsf(h,obincount,sfbincount,phasecount)
%
function hs=kernsmoothorsfold(h,obincount,sfbincount,phasecount)

if ~exist('phasecount'),
   phasecount=4;
end
if ~exist('obincount'),
   obincount=12;
end
if ~exist('sfbincount'),
   sfbincount=8;
end
s=size(h);
spacebincount=s(1);
tbincount=s(2);
kerncount=prod(s(3:end));
h=reshape(h,spacebincount,tbincount,kerncount);
chancount=spacebincount/phasecount;

Xmax=sqrt(chancount*2);
[cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);

tsf=zeros(Xmax*Xmax,tbincount,phasecount,kerncount);
for kidx=1:kerncount,
   for ii=1:phasecount,
      tsf(cfilt,:,ii,kidx)=h((1:chancount)+(ii-1)*chancount,:,kidx);
   end
end
tsf(cfiltconj,:,:,:)=tsf(cfilt,:,:,:);

tsf=reshape(tsf,Xmax,Xmax,tbincount,phasecount,kerncount);
tsf=sf2gr(tsf,obincount,sfbincount,0,0);
hs=gr2sf(tsf,(Xmax-2)/2,Xmax,Xmax);
hs=permute(reshape(hs,Xmax*Xmax,tbincount,phasecount,kerncount),[1 3 2 4]);
hs=reshape(hs(cfilt,:,:,:),s);

