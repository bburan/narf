% function [mov,transmov]=kern2opt(kern,movformat,framecount,dorand);
%
% generate random frames "optimal" to kern, normalize to global
% mean and std
%
% created SVD 6/5/02
%
function [mov,transmov]=kern2opt(kern,movformat,framecount,dorand);

% randomize here:
rand('state',sum(100*clock));
randstate=rand('state');

kernfmt=movformat.kernfmt;
stimwindowsize=movformat.cellfiledata(1).stimwindowsize;

if ~exist('dorand'),
   dorand=0;
end

if strcmp(kernfmt,'wg') | strcmp(kernfmt,'wav'),
   spacebincount=movformat.stimiconside(1);
   obincount=movformat.stimiconside(3);
   sbincount=movformat.stimiconside(4);
   
   [idxset,scalefactor,iconside]=wavidxset(stimwindowsize,spacebincount,...
                                           obincount,sbincount);
   mov0=rand(stimwindowsize,stimwindowsize,framecount);
   mov=zeros(size(mov0));
   [tmov,pind]=buildSFpyr(mov0(:,:,1),sbincount,obincount-1);
   wmov=zeros(length(tmov),framecount);
   transmov=zeros(length(idxset),framecount);
   %wmov(:,1)=tmov;
   f=cat(1,idxset{:});
   
   for ii=1:framecount,
      tmov=buildSFpyr(mov0(:,:,ii),sbincount,obincount-1);
      wmov(f,ii)=tmov(f);
   end
   for ii=1:size(kern,1),
      wmov(idxset{ii},:)=wmov(idxset{ii},:).*kern(ii).*(kern(ii)>0);
      transmov(ii,:)=scalefactor(ii).*sqrt(mean(wmov(idxset{ii},:).^2,1));
   end
   for ii=1:framecount,
      mov(:,:,ii)=reconSFpyr(wmov(:,ii),pind);      
   end
   
   for ii=1:min([size(mov0,3) 5]),
      figure(1);
      subplot(5,2,ii*2-1);
      imagesc(mov0(:,:,ii));
      subplot(5,2,ii*2);
      imagesc(mov(:,:,ii));
   end
   
   
elseif strcmp(kernfmt,'fft') | strcmp(kernfmt,'fftgr') | ...
      strcmp(kernfmt,'pfft') | strcmp(kernfmt,'pfftgr'),
   if strcmp(kernfmt,'fft') | strcmp(kernfmt,'fftgr'),
      phasecount=4;
   else
      phasecount=1;
   end
   
   chancount=size(kern,1)/phasecount;
   Xmax=sqrt(chancount*2);
   [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);
   
   skern=repmat(sum(reshape(kern,chancount,phasecount),2),[4 1]);
   skern=skern-min(skern);  % avoid negative terms in skern
   skernp=repmat(sum(reshape(kern.*(kern>0),chancount,phasecount),2),[4 1]);
   
   dssize=Xmax;
   smov0=rand(dssize,dssize,framecount);
   fkern=zeros(Xmax,Xmax);
   fmov=zeros(Xmax*Xmax,framecount);
   fmovp=zeros(Xmax*Xmax,framecount);
   
   for ii=1:framecount,
      tf=fftshift(fft2(smov0(:,:,ii)));
      fmov(cfilt,ii)=tf(cfilt).*skern(1:chancount);
      fmov(cfiltconj,ii)=tf(cfiltconj).*skern(1:chancount);
      fmovp(cfilt,ii)=tf(cfilt).*skernp(1:chancount);
      fmovp(cfiltconj,ii)=tf(cfiltconj).*skernp(1:chancount);
   end
   
   % skip phase separated crap... just scale amplitudes according
   % to power of positive kernel
   % lame and hacky... i know.
   if 0,  
   [tmov,tmovpower]=movphasesep(smov0,movformat.stimfilterparms{:});
   tmovp=zeros(size(tmov));
   for ii=1:framecount,
      tmovp(:,ii)=tmov(:,ii).* skern; %  .* (kern > 0);
      tmov(:,ii)=tmov(:,ii) .* kern;
   end
   
   for ii=1:phasecount,
      fmov(cfilt,:)=fmov(cfilt,:)+ ...
          tmov((1:chancount)+(ii-1)*chancount,:) .* i.^(ii-1);
      fmovp(cfilt,:)=fmovp(cfilt,:)+ ...
          tmovp((1:chancount)+(ii-1)*chancount,:) .* i.^(ii-1);
      fkern(cfilt)=fkern(cfilt)+kern((1:chancount)+(ii-1)*chancount);
   end
   %fmov(cfiltconj,:)=fmov(cfilt,:);
   end
   
   fmov=reshape(fmov,Xmax,Xmax,framecount);
   fmovp=reshape(fmovp,Xmax,Xmax,framecount);
   
   smov=zeros(Xmax,Xmax,framecount);
   smovp=zeros(Xmax,Xmax,framecount);
   for ii=1:framecount,
      smov(:,:,ii)=real(ifft2(fftshift(fmov(:,:,ii))));
      smovp(:,:,ii)=real(ifft2(fftshift(fmovp(:,:,ii))));
   end
   
   mov=zeros(stimwindowsize,stimwindowsize,framecount);
   
   for ii=1:framecount,
      mov(:,:,ii)=imresize(smovp(:,:,ii),[stimwindowsize stimwindowsize],...
                           'bilinear');
   end
   transmov=reshape(fmovp,Xmax*Xmax,framecount);
end

%normalize mean and variance

for ii=1:framecount,
   tmov=mov(:,:,ii);
   mov(:,:,ii)=(tmov-mean(mean(tmov(:))))./std(tmov(:));
end

GLOBALMEAN=116;
GLOBALSTD=40;

mov=mov.*GLOBALSTD+GLOBALMEAN;
mov(find(mov>255))=255;
mov(find(mov<0))=0;

return

for ii=1:5,
   figure(1);
   subplot(5,3,ii*3-2);
   imagesc(mov0(:,:,ii));
   subplot(5,3,ii*3-1);
   imagesc(mov(:,:,ii));
   %subplot(5,3,ii*3);
   %imagesc(movp(:,:,ii));
end
