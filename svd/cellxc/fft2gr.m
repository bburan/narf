% function hs=fft2gr(h,obincount,sfbincount,phasecount,stepcount)
%
function hs=fft2gr(h,obincount,sfbincount,phasecount,stepcount)

if ~exist('obincount'),
   obincount=12;
end
if ~exist('sfbincount'),
   sfbincount=8;
end
if ~exist('phasecount'),
   phasecount=4;
end
if ~exist('stepcount'),
   stepcount=7;
end

s=size(h);
spacebincount=s(1);
chancount=spacebincount/phasecount;
framecount=prod(s)./chancount;

h=reshape(h,chancount,framecount);
Xmax=sqrt(chancount*2);
[cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);

sf2grfilt=sf2gr0(Xmax,obincount,sfbincount,stepcount);

% take only relevant bins
sf2grfilt=sf2grfilt(cfilt,:);

hs=sf2grfilt'*h;

% smoothing was already done in kernsmoothorsf?????
SMOOTHON=0;
if SMOOTHON,
   % if SMOOTHON is selected, apply a 2D gaussian filter to the
   % OR-SF domain representation of the function being transformed
   xx=(-1:1)';
   gsf=exp(-(xx./0.6).^2/2);
   gsf=gsf./sum(gsf(:));
   gor=exp(-(xx./0.6).^2/2);
   gor=(gor./sum(gor(:)));
   
   hs=reshape(hs,obincount,prod([sfbincount,phasecount,s(2:end)]));
   hs=cconv2(hs,gor);
   
   hs=reshape(hs,obincount,sfbincount,prod([phasecount,s(2:end)]));
   hs=permute(hs,[2 1 3]);
   hs=reshape(hs,sfbincount,prod([obincount,phasecount,s(2:end)]));
   hs=conv2(hs,gsf,'same');
   
   hs=reshape(hs,[sfbincount,obincount,phasecount,s(2:end)]);
   hs=permute(hs,[2 1 4 3 5 6]);
   
else
   hs=reshape(hs,[obincount,sfbincount,phasecount,s(2:end)]);
   hs=permute(hs,[1 2 4 3 5 6]);
end

return


h=zeros(size(cfilt));
h(50)=1;

subplot(1,3,1);
%imagesc(h);
plot(h);

subplot(1,3,2);
hgr=sf2grfilt'*h(:);
%imagesc(reshape(hgr,obincount,sfbincount));
plot(hgr);

subplot(1,3,3);
hout=gr2sffilt'*hgr;
%imagesc(reshape(hout,Xmax,Xmax));
plot(hout);

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

