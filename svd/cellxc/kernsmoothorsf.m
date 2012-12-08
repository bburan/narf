% function hs=kernsmoothorsf(h,obincount,sfbincount,phasecount,stepcount,...
%                            sigor,sigsf)
%
function hs=kernsmoothorsf(h,obincount,sfbincount,phasecount,stepcount,...
                           sigor,sigsf)

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
if ~exist('sigor'),
   sigor=0.75;
end
if ~exist('sigsf'),
   sigsf=0.75;
end

s=size(h);
spacebincount=s(1);
chancount=spacebincount/phasecount;
framecount=prod(s)./chancount;

h=reshape(h,chancount,framecount);
Xmax=sqrt(chancount*2);
[cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);

sf2grfilt=sf2gr0(Xmax,obincount,sfbincount,stepcount);
gr2sffilt=gr2sf0(obincount,sfbincount,Xmax,stepcount);

% take only relevant bins
sf2grfilt=sf2grfilt(cfilt,:);
gr2sffilt=gr2sffilt(:,cfilt)*2;

hs=sf2grfilt'*h;

% smooth in orsf domain????
SMOOTHON=1;
if SMOOTHON,
   fprintf('S');
   xx=(-1:1);
   gor=exp(-(xx./sigor).^2/2);
   gor=(gor./sum(gor(:)))';
   gsf=exp(-(xx./sigsf).^2/2);
   gsf=gsf./sum(gsf(:));
   
   hs=reshape(hs,obincount,sfbincount*framecount);
   hs=cconv2(hs,gor);
   
   hs=reshape(hs,obincount,sfbincount,framecount);
   hs=convn(hs,gsf,'same');
   
   hs=reshape(hs,sfbincount*obincount,framecount);
end

hs=gr2sffilt'*hs;

hs=reshape(hs,s);

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

