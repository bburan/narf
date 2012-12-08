% function hgr=fft2gr(hf,obincount,sfbincount,phasecount);
%
% generate an orientation/spatial frequency-binned gratrev kernel
% from a phase-sep fourier domain function
%
% parameters:
%
% created SVD 8/2/02 - hacked from sf2gr.m
%
function hgr=fft2gr(hf,obincount,sfbincount,phasecount);

%INTMETH='nearest';
INTMETH='linear';
%INTMETH='cubic'; 

s=size(hf);
spacebincount=s(1);
timebincount=s(2);
kerncount=prod(s(3:end));
chancount=spacebincount/phasecount;
Xmax=sqrt(chancount*2);

if ~exist('obincount'),
   obincount=12;
end
if ~exist('sfbincount'),
   sfbincount=Xmax/2-1;
end
if ~exist('phasecount'),
   phasecount=4;
end
if ~exist('brotate90'),
   brotate90=0;
end

% convert to 2D fourier domain
[cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);
hf=reshape(hf,spacebincount/phasecount,phasecount,timebincount,kerncount);
sfIR=zeros(Xmax*Xmax,timebincount,phasecount,kerncount);
sfIR(cfilt,:,:,:)=permute(hf,[1 3 2 4]);
sfIR(cfiltconj,:,:,:)=sfIR(cfilt,:,:,:);
sfIR=reshape(sfIR,Xmax,Xmax,timebincount,phasecount,kerncount);

xc=round((Xmax+1)/2);

obins=linspace(0,pi,obincount+1);
obins=obins(1:end-1);
orstep=obins(2)-obins(1);

sfbins=linspace(1,round((Xmax+1)/2),sfbincount+1);
sfbins=sfbins(1:end-1);

[SF,OR]=meshgrid(sfbins,obins);

XX=cos(OR).*SF;
YY=-sin(OR).*SF;
% use these to smooth in or as necessary
XXa=cos(OR-orstep/3).*SF;
YYa=-sin(OR-orstep/3).*SF;
XXb=cos(OR+orstep/3).*SF;
YYb=-sin(OR+orstep/3).*SF;

WX=(1:Xmax)-xc;
WY=(1:Xmax)-xc;

hgr=zeros(obincount,sfbincount,timebincount,phasecount,kerncount);

for tt=1:timebincount,
   for ridx=1:kerncount,
      for phidx=1:phasecount,
         hgr(:,:,tt,phidx,ridx)= SF.* ...
             (0.5*interp2(WX,WY,sfIR(:,:,tt,phidx,ridx),XX,YY,INTMETH) + ...
              0.25*interp2(WX,WY,sfIR(:,:,tt,phidx,ridx),XXa,YYa,INTMETH) + ...
              0.25*interp2(WX,WY,sfIR(:,:,tt,phidx,ridx),XXb,YYb,INTMETH));
      end
   end
end

if brotate90,
   hgr=fftshift(hgr,1);
end

hgr(find(isnan(hgr)))=0;
hgr=reshape(hgr,obincount,sfbincount,timebincount,phasecount,s(3:end));

return

if ~exist('sumstr'),
   sumstr=sprintf('sf2gr.m: obins=%d (down) sfbins=%d (left)',...
                  obincount,sfbincount);
end

movdecorrres(hgr,[0 90 180 270 0],sumstr,1,min([timebincount 8]),...
             {'0','90','180','270','sum'},0);

return

%
% skip these displays for the time being
%
movdecorrres(H,[0 90 180 270 0],sumstr,1,min([timebincount 8]),...
             {'0','90','180','270','sum'},0);

figure(2)
subplot(2,2,1);
sfor=sum(H(:,:,:,5),3);
kmax=max(abs(sfor(:)));
imagesc(sfor,[-kmax kmax]);
title('sf vs. or');
xlabel('sf');
ylabel('or');

subplot(2,2,2);
timeor=squeeze(sum(H(:,:,:,5),2));
kmax=max(abs(timeor(:)));
imagesc(timeor,[-kmax kmax]);
title('time vs. or');
xlabel('time');
ylabel('or');

subplot(2,2,3);
timesf=squeeze(sum(H(:,:,:,5),1));
kmax=max(abs(timesf(:)));
imagesc(timesf,[-kmax kmax]);
title('time vs. sf');
xlabel('time');
ylabel('sf');

subplot(2,2,4);
orphase=squeeze(sum(sum(H(:,:,:,1:4),2),3))';
kmax=max(abs(orphase(:)));
imagesc(orphase,[-kmax kmax]);
title('or vs. phase');
xlabel('or');
ylabel('phase');


