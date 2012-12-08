% function [H_out,obins,rbins]=car2pol_strf(H,obincount,rbincount,kernfmt);
%
% generate an orientation/frequecy plot from an STRF
%
%
% parameters:
% H: the strf. exprected format depends on value of kernfmt:
%    kernfmt='strf' (default): H should be in freq X time format.
%    kernfmt='spect': H is MTF (scale X rate) format.
% obincount: number of orientation bins (default [8])
% rbincount: number of frequency bins (default [size(H,1)./2])
%
% returns: 
% H_out in orientation X frequency format
%
% created SVD 2007-04-26, ripped off of sf2gr.m
%
function [H_out,obins,rbins]=car2pol_strf(H,obincount,rbincount,kernfmt,...
                                          bphasemean,brotate90);

%INTMETH='nearest';
INTMETH='linear';
%INTMETH='cubic'; 

if ~exist('kernfmt','var'),
   kernfmt='strf';
end
if strcmp(kernfmt,'fft') | strcmp(kernfmt,'sfft') | ...
      strcmp(kernfmt,'pfft') | findstr(kernfmt,'pfft+') | ...
      strcmp(kernfmt,'lin2'),
   
   sfsize=size(H);
   spacebincount=sfsize(1);
   
   if strcmp(kernfmt,'fft') | strcmp(kernfmt,'lin2'),
      phasecount=4;
   elseif length(kernfmt)>4 & findstr(kernfmt,'pfft+'),
      phasecount=str2num(strf.parms.kernfmt(end));
   else
      phasecount=1;
   end
   
   chancount=spacebincount/phasecount;
   Xmax=sqrt(chancount*2);
   [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);
   
   tsf=zeros([Xmax*Xmax,sfsize(2:end)]);
   tsf(cfilt,:,:)=squeeze(sum(reshape(H,[chancount,phasecount,...
                                      sfsize(2:end)]),2));
   tsf(cfiltconj,:,:)=tsf(cfilt,:,:);
   tsf=reshape(tsf,[Xmax,Xmax,sfsize(2:end)]);
   
   H=tsf;
   
   phasecount=1; % fix to 1 since we summed over phase
elseif strcmp(kernfmt,'strf'),
   
   % convert from raw STRF (freq X time) format to MTF
   H=abs(fftshift(fft2(H)));
   
else
   % assume already given a modulation transfer function
   
   % deal with other kernfmts!
   %%disp('kern2tune.m requires strfs with fft or pfft kernfmt! sorry');
   %return
end

[Xmax,Ymax,timebincount,phasecount,respcount]=size(H);

if Xmax~=Ymax,
   H=imresize(H,[1 1].*max(Xmax,Ymax));
   [Xmax,Ymax]=size(H);
end

xc=round((Xmax+1)/2);
yc=round((Ymax+1)/2);
Hsum=sum(H,4);

if ~exist('obincount','var'),
   obincount=8;
end
if ~exist('rbincount','var'),
   rbincount=min(Xmax,Ymax)/2-1;
end
if ~exist('bphasemean','var'),
   bphasemean=0;
end
if ~exist('brotate90','var'),
   brotate90=0;
end

obins=linspace(0,pi,obincount+1);
obins=obins(1:end-1);
orstep=obins(2)-obins(1);

rbins=linspace(1,round((Xmax+1)/2),rbincount+1);
rbins=rbins(1:end-1);

[SF,OR]=meshgrid(rbins,obins);

XX=cos(OR).*SF;
YY=-sin(OR).*SF;
% use these to smooth in or as necessary
XXa=cos(OR-orstep/3).*SF;
YYa=-sin(OR-orstep/3).*SF;
XXb=cos(OR+orstep/3).*SF;
YYb=-sin(OR+orstep/3).*SF;

WX=(1:Xmax)-xc;
WY=(1:Ymax)-xc;

XX(XX>max(WX))=-XX(XX>max(WX));
XXa(XXa>max(WX))=-XXa(XXa>max(WX));
XXb(XXb>max(WX))=-XXb(XXb>max(WX));

H_out=zeros(obincount,rbincount,timebincount,phasecount+bphasemean,respcount);

for tt=1:timebincount,
   for ridx=1:respcount,
      for phidx=1:phasecount,
         H_out(:,:,tt,phidx,ridx)= ...
            (0.5*interp2(WX,WY,H(:,:,tt,phidx,ridx),XX,YY,INTMETH) + ...
              0.25*interp2(WX,WY,H(:,:,tt,phidx,ridx),XXa,YYa,INTMETH) + ...
              0.25*interp2(WX,WY,H(:,:,tt,phidx,ridx),XXb,YYb,INTMETH));
      end
      if bphasemean,
         H_out(:,:,tt,phasecount+1,ridx)=...
             interp2(WX,WY,Hsum(:,:,tt,1,ridx),XX,YY).*SF;
      end
   end
end

if brotate90,
   H_out=fftshift(H_out,1);
end

% convert orientation bins to degrees
obins = obins .* 180/pi;

H_out(find(isnan(H_out)))=0;

