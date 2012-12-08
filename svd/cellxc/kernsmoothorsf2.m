% function hs=kernsmoothorsf2(h,obincount,sfbincount,kernfmt)
%
function hs=kernsmoothorsf2(h,obincount,sfbincount,kernfmt)

if ~exist('obincount','var'),
   obincount=12;
end
if ~exist('sfbincount','var'),
   sfbincount=8;
end
if ~exist('kernfmt','var'),
   kernfmt=pfft;
end
if strcmp(kernfmt,'fft') | strcmp(kernfmt,'lin2'),
   phasecount=4;
elseif length(kernfmt)>4 & ...
      (kernfmt(end-1)=='+' | kernfmt(end-1)=='-'),
   phasecount=str2num(kernfmt(end));
else
   phasecount=1;
end

s=size(h);
spacebincount=s(1);
chancount=spacebincount/phasecount;
framecount=prod(s./chancount);

h=reshape(h,s(1)./phasecount,phasecount*prod(s(2:end)));

Xmax=sqrt(chancount*2);
[cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);

sf2grfilt=sf2gr0(Xmax,obincount,sfbincount);
gr2sffilt=gr2sf0(obincount,sfbincount,Xmax);

% take only relevant bins
sf2grfilt=sf2grfilt(cfilt,:);
gr2sffilt=gr2sffilt(:,cfilt)*2;

hs=sf2grfilt'*h;
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

