% function r=xcpredpic(cellid,batch,nlidx,cnfidx[=thisbatch])
%
% load results from resfile in sRunData and display fits, pred
% results
%
% r=0 if no entries found in db, =1 otherwise
%
function r=xcpredpic(runidx,batch,nlidx0,cnfidx)

dbopen;
if isnumeric(runidx),
   sql=['SELECT * from sRunData WHERE id=',num2str(runidx)];
   rundata=mysql(sql);
   cellid=rundata.cellid;
else
   cellid=runidx;
   goodbatch=zeros(1,length(batch));
   sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
        ' AND batch=',num2str(batch)];
   rundata=mysql(sql);
end

if length(rundata)==0,
   disp('no entries found in db!');
   if nargout>0,
      r=0;
   end
   return
end

resfile=[rundata(1).respath,rundata(1).resfile,'.gz'];
if strcmp(resfile(end-2:end),'.gz'),
   zload(resfile);
else
   load(resfile);
end

global BATQUEUEID
BATQUEUEID=[];

if ~exist('nlidx0','var'),
   nlidx=params.nlidxsave;
else
   nlidx=nlidx0;
end

if nlidx>length(strf),
   nlidx=1;
end

if ~exist('cnfidx','var'),
   cnfidx=find(cat(1,params.predbatch{:})==params.batch);
   if isempty(cnfidx),
      cnfidx=1;
   end
end

fprintf('chose to show nlidx=%d cnfidx=%d\n',nlidx,cnfidx);

[pcellfiledata,ptimes,pbatchdata]=...
    cellfiletimes(params.cellid,params.predbatch{cnfidx});
predparams=params;
predparams.stimfiles={};
predparams.respfiles={};
predparams.stimcrfs=[];
for ii=1:length(pcellfiledata),
   predparams.stimfiles{ii}=[pcellfiledata(ii).stimpath,...
                    pcellfiledata(ii).stimfile];
   predparams.respfiles{ii}=[pcellfiledata(ii).path,...
                    pcellfiledata(ii).respfile];
   if pcellfiledata(ii).stimfilecrf>0,
      predparams.stimcrfs(ii)=pcellfiledata(ii).stimfilecrf;
   else
      tstimpix=strsep(pcellfiledata(ii).stimiconside,',');
      if length(tstimpix)>0,
         tstimpix=tstimpix{1};
      end
      sql=['SELECT * FROM gCellMaster WHERE cellid="',...
           params.cellid,'"'];
      celldata=mysql(sql);
      predparams.stimcrfs(ii)=tstimpix./celldata.rfsize;
   end
end

tpredstartframe=ptimes(3).start;
tpredstopframe=ptimes(3).stop;
tpredfile=ptimes(3).fileidx;
predparams.times=ptimes;

% load the stim/resp data
[cdata.stim,cdata.resp]=xcloadstimresp(tpredfile,tpredstartframe,...
                                       tpredstopframe,predparams);

fparams=predparams;
fparams.stimloadparms={0,0,32,0,1};
fparams.stimfiltercmd='';
[fdata.stim,fdata.resp]=xcloadstimresp(tpredfile,tpredstartframe,...
                                       tpredstopframe,fparams);
ostim=reshape(fdata.stim',32,32,size(fdata.stim,1));
clear fdata

hf=strf(nlidx).h;
ir=strf(nlidx).tempresp;
ir=ir(1:8);
maxlag=min(find(max(ir)==ir))-1;

ir2=ir([1:maxlag maxlag+2:end]);
maxlag2=min(find(max(ir2)==ir2))-1;
if maxlag2>=maxlag,
   maxlag2=maxlag2+1;
end

rel1=ir(maxlag+1)./(ir(maxlag+1)+ir(maxlag2+1));
rel2=ir(maxlag2+1)./(ir(maxlag+1)+ir(maxlag2+1));
fprintf('maxlag=%d (%.2f) maxlag2=%d (%.2f)\n',...
        maxlag,rel1,maxlag2,rel2);

if isfield(strf(nlidx),'hspace'),
   rr=cdata.resp;
   
   ir=strf(nlidx).tempresp;
   
   decr1=shift(rr,-maxlag);
   decr1(end-maxlag+1:end)=0;
   decr2=shift(rr,-maxlag2);
   decr2(end-maxlag2+1:end)=0;
   decr=decr1.*rel1+decr2.*rel2;
   
   cdata.resp0=cdata.resp;
   cdata.resp=decr;
   
   global ESTIMATIONPHASE VALIDATIONPHASE
   VALIDATIONPHASE=0; ESTIMATIONPHASE=1;
   tstrf=repmat(strf(nlidx),[3 1]);
   hs=tstrf(2).hspace;
   tr=std(tstrf(2).tempresp);
   hsp=hs.*(hs>0);
   hsn=hs.*(hs<0);
   tstrf(1).h=hs.*tr;
   tstrf(1).nltype='none';
   tstrf(2).h=hsp.*tr;
   tstrf(2).nltype='none';
   tstrf(3).h=hsn.*tr;
   tstrf(3).nltype='none';
   predres=xcval(tstrf,predparams,cdata);
   
   % don't need to do any more temporal manipulation since ir has
   % already been removed
   pp=predres.mod_psth{1}(:,1);
   ppp=predres.mod_psth{1}(:,2);
   ppn=predres.mod_psth{1}(:,3);
   rr=predres.act_resp{1}(:,1);
   pp(find(isnan(pp)))=0;
   rr(find(isnan(rr)))=0;
   decp=pp;
   decr=rr;
   
   nzidx=find(ppp>min(ppp));
   [junk,rri]=sortrows([-decr(nzidx) -decp(nzidx)]);
   
   figure(5);
   imagesc(1:length(ppp),1:length(ppp),...
           genpredpic(ppp(nzidx),ppn(nzidx),ostim(:,:,nzidx),rri(1:40)));
   xlabel('pos pred');
   ylabel('neg pred');
   title(sprintf('cell %s pos vs. neg',cellid));
   axis image
   fullpage('portrait');
   
   % normalize to make plots look good
   decp=decp./std(decp).*std(decr);
   pp=pp./std(pp).*std(rr);
   
else
   predres=xcval(strf(nlidx),predparams,cdata);

   pp=predres.mod_psth{1}(:,1);
   rr=predres.act_resp{1}(:,1);
   
   bidx=find(~isnan(pp) & ~isnan(rr));
   bidx2=bidx+maxlag2-maxlag;
   bidx1=bidx(bidx2>0 & bidx2<length(pp));
   bidx2=bidx2(bidx2>0 & bidx2<length(pp));
   bidx1=bidx1(find(~isnan(pp(bidx2)+rr(bidx2))));
   bidx2=bidx2(find(~isnan(pp(bidx2)+rr(bidx2))));
   
   pp(find(isnan(pp)))=0;
   rr(find(isnan(rr)))=0;
   

   if 0 ,
      % deconvolve responses with inverse of temporal response
      
      ir=zeros(1,10);
      ir(maxlag+1)=strf(nlidx).tempresp(maxlag+1);
      ir(maxlag2+1)=strf(nlidx).tempresp(maxlag2+1);
      
      hir=fft(ir);
      ihir=hir./(abs(hir).^2+(abs(hir)==0));
      fir=real(ifft(ihir));
      
      % do full deconvolution
      decp=conv(pp,flipud(fir(:)));
      decr=conv(rr,flipud(fir(:)));
      
      %decp=decp(1:(end-length(fir)+1));
      %decr=decr(1:(end-length(fir)+1));
      decp=decp(length(fir):end);
      decr=decr(length(fir):end);
      
      decp=decp./mean(decp).*mean(pp);
      decr=decr./mean(decr).*mean(rr);
      
   else
      % simply "deconvolve" by shifting stimulus in time according to
      % peak of the temporal impulse response
      
      z1=zeros([size(ostim(:,:,1)),maxlag]);
      z2=zeros([size(ostim(:,:,1)),maxlag2]);
      
      % only peak frame
      ostim=cat(3,z1,ostim(:,:,1:(end-maxlag)));
      ostim=ostim(:,:,bidx1);
      
      pp=(pp(bidx1)+pp(bidx2))./2;
      rr=(rr(bidx1)+rr(bidx2))./2;
      %pp=pp(bidx);
      %rr=rr(bidx);
      
      decp=pp;
      decr=rr;
   end
end

[junk,ppi]=sortrows([-decp -decr]);
[junk,rri]=sortrows([-decr -decp]);
%[junk,ppi]=sort([-decp]);
%[junk,rri]=sort([-decr]);




%
% plot psth and sorted responses
%
figure(1);
clf

subplot(3,1,1);
plot(rr,'k--');
hold on
plot(pp,'r');
hold off
legend('obs','pred');
title(sprintf('%s optinat analysis. predcorr=%.2f',...
              params.cellid,predres.predxc));

subplot(3,1,2);
plot(decr,'k--');
hold on
plot(decp,'r');
hold off
legend('deconv obs','deconv pred');
ylabel(sprintf('cc=%.2f',xcov(decr,decp,0,'coeff')));

subplot(3,1,3);
plot(decp(ppi),'k--');
hold on
plot(decr(ppi),'r');
hold off
legend('pred ranked','obs');

%
% compute angle stuff
%
if 1,
   figure(2);
   clf
   
   h=strf(nlidx).hspace;
   showkern([h h],strf(nlidx).parms.kernfmt);
   
   subplot(1,2,2);
   cla
   plot(0:16:16*(length(strf(nlidx).tempresp)-1),strf(nlidx).tempresp);
   axis square
      
   %keyboard
else
   h=strf(nlidx).hspace;
   spacecount=length(h);
   mS=strf(nlidx).mS;
   tstim=cdata.stim'-repmat(mS,1,size(cdata.stim,1));
   
   linpred=tstim'*h;
   
   tstim(find(abs(h)==0),:)=0;
   
   snorm=sqrt(sum(tstim.^2,1))';
   snorm(find(snorm==0))=1;
   angle=acos(linpred./norm(h(:))./snorm) *180/pi;
   
   figure(2);
   clf
   
   scatter(angle,decp-decr);
   xlabel('angle');
   ylabel('pred-obs resp');
end

figure(3);
clf

[bigplot,rstr]=genpredpic(decp,decr,ostim);

if strcmp(rstr,'rank') | strcmp(rstr,'rank2'),
   bpsize=size(bigplot,1);
   showcount=length(decp);
   hc=imagesc(linspace(1,showcount,bpsize),...
              linspace(1,showcount,bpsize),bigplot);
else
   rrange=max([tp; tr])-rmin;
   hc=imagesc(linspace(rmin,rmin+rrange,bpsize),...
              linspace(rmin,rmin+rrange,bpsize),bigplot);
end
ha=get(hc,'Parent');
tl=get(ha,'XTickLabel');
set(ha,'XTickLabel',flipud(tl));

xlabel(['pred ' rstr]);
ylabel(['actual ' rstr]);
axis image
colormap(gray);

title(sprintf('cell %s batch %d nlidx %d cnfidx %d predxc %.2f',...
              params.cellid,batch,nlidx,cnfidx,predres.predxc(1)));
set(gcf,'PaperPosition',[0.25 0.25 8 10.5],...
        'PaperOrientation','portrait');

if 1,
   figure(4);
   clf
   showcount=20;
   for ii=1:20,
      
      subplot(4,showcount/2,ii);
      if ppi(ii)<size(ostim,3) & ppi(ii)>1,
         sstim=cat(1,ostim(:,:,ppi(ii)-1),ostim(:,:,ppi(ii)),...
                   ostim(:,:,ppi(ii)+1));
      elseif ppi(ii)<size(ostim,3)
         sstim=cat(1,ostim(:,:,ppi(ii)),ostim(:,:,ppi(ii)+1));
      else
         sstim=cat(1,ostim(:,:,ppi(ii)-1),ostim(:,:,ppi(ii)));
      end
      
      imagesc(sstim,[0 255]);
      axis image
      axis off
      title(sprintf('pred %d v %d',ii,find(ppi(ii)==rri)));
      
      subplot(4,showcount/2,ii+showcount);
      if rri(ii)<size(ostim,3) & rri(ii)>1,
         sstim=cat(1,ostim(:,:,rri(ii)-1),ostim(:,:,rri(ii)),...
                   ostim(:,:,rri(ii)+1));
      elseif rri(ii)<size(ostim,3)
         sstim=cat(1,ostim(:,:,rri(ii)),ostim(:,:,rri(ii)+1));
      else
         sstim=cat(1,ostim(:,:,rri(ii)-1),ostim(:,:,rri(ii)));
      end
      imagesc(sstim,[0 255]);
      %imagesc(ostim(:,:,rri(ii)),[0 255]);
      axis image
      axis off
      title(sprintf('act: %d v %d',ii,find(rri(ii)==ppi)));
   end
   colormap(gray)
   
   if 0,
      figure(5);
      clf
      subplot(2,2,1);
      hc=imagesc(bigplot(1:bpsize/2,1:bpsize/2,:));
      axis off
      axis image
      subplot(2,2,2);
      hc=imagesc(bigplot(1:bpsize/2,bpsize/2+1:bpsize,:));
      axis off
      axis image
      subplot(2,2,3);
      hc=imagesc(bigplot(bpsize/2+1:bpsize,1:bpsize/2,:));
      axis off
      axis image
      subplot(2,2,4);
      hc=imagesc(bigplot(bpsize/2+1:bpsize,bpsize/2+1:bpsize,:));
      axis off
      axis image
   end
end



function [bigplot,rstr]=genpredpic(r1,r2,stim,hiliteidx,kernfmt);

DOPREDERR=0;
DORANK=1;    % 0-plot response strength
             % 1-plot ranks 
             % 2-plot some weird error thing vs actual resp rank
DOTRANS=0;   % ie, do pfft of images for display

decp=r1;
decr=r2;

[junk,ppi]=sortrows([-decp -decr]);
[junk,rri]=sortrows([-decr -decp]);

prank=ones(size(ppi));
prank(ppi)=(length(ppi):-1:1);
rrank=ones(size(ppi));
rrank(rri)=(length(rri):-1:1);

decb=rrank; % (prank+rrank)./2;
dece=(rrank-prank)./2;
%decb=(decp+decr)./2;
%dece=(decp-decr)./(decb+1);
[junk,bbi]=sort(-decb);
[junk,eei]=sort(-dece);

if DOTRANS,
   psize=16;
   psize2=16;
else
   psize=size(stim,1);
   psize2=size(stim,2);
end
bpsize=30*psize;
pplen=length(ppi);
bigplot=uint8(ones(bpsize,bpsize,3).*255);
showcount=pplen;

if DOPREDERR,
   dece=decr-decp;
   rmin=min([dece; decr]);
   tp=(dece-rmin);
   tr=(decr-rmin);
   trmin=0;
else
   rmin=min([decp; decr]);
   tp=(decp-rmin);
   tr=(decr-rmin);
   trmin=0;
end

rrange=max([tp; tr])-rmin;

if ~exist('hiliteidx','var') | isempty(hiliteidx),
   iirange=showcount:-1:1;
   hiliteidx=[];
else
   iirange0=showcount:-1:1;
   iirange1=iirange0(find(~ismember(ppi(iirange0),hiliteidx(:)')));
   iirange2=iirange0(find(ismember(ppi(iirange0),hiliteidx(:)')));
   iirange=[iirange1 iirange2];
end

for ii=iirange,
   
   if DORANK==1,
      pidx=ppi(ii);
      xpos=round((bpsize-psize)./pplen*(find(rri==ppi(ii))-1));
      %ypos=round((bpsize-psize2)./pplen*(ii-1));
      ypos=round((bpsize-psize2)./pplen*(showcount-ii));
      rstr='rank';
   elseif DORANK==2,
      pidx=bbi(ii);
      %xpos=round((bpsize-psize) ./ rrange .* (tp(pidx)-trmin));
      %ypos=round((bpsize-psize2) ./ rrange .* (tr(pidx)-trmin));
      
      xpos=round((bpsize-psize)./pplen*(ii-1));
      ypos=round((bpsize-psize2)./pplen*(find(eei==bbi(ii))-1));
      rstr='rank2';
   else
      pidx=ppi(ii);
      xpos=round((bpsize-psize) ./ rrange .* (tp(pidx)-trmin));
      ypos=round((bpsize-psize2) ./ rrange .* (tr(pidx)-trmin));
      rstr='resp';
   end
   
   if DOTRANS & (strcmp(kernfmt,'fft') | strcmp(kernfmt,'pfft')),
      
      fim=cdata.stim(pidx,:);
      %fim=cdata.stim(pidx,:)'-mS;
      %fim=fim-h./norm(h).^2 .* decp(pidx).^2;
      
      if strcmp(kernfmt,'fft'),
         phasecount=4;
      else
         phasecount=1;
      end
      
      spacebincount=size(cdata.stim,2);
      chancount=spacebincount/phasecount;
      Xmax=sqrt(chancount*2);
      [cfilt,cfiltconj,cmask,cmaskconj]=gencfilt(Xmax,Xmax);
      tsf=zeros(Xmax,Xmax);
      tsf(cfilt)=sum(reshape(fim,chancount,phasecount),2);
      tsf(cfiltconj)=tsf(cfilt);
      tsf=reshape(tsf,Xmax,Xmax);
      psize=size(tsf,1);
      im=tsf;
      %im=sqrt(tsf-min(tsf(:)));
      im(find(im(:)<0))=0;
      im=im.^(0.33);
      if max(im(:))>0,
         im=im./max(im(:)).*512;
      end
   else
      im=stim(:,:,pidx);
   end
   
   if xpos<0 | ypos<0,
      keyboard;
   end
   
   bigplot((xpos+1):(xpos+psize),(ypos+1):(ypos+psize2),1)=im;
   bigplot((xpos+1):(xpos+psize),(ypos+1):(ypos+psize2),2)=im;
   bigplot((xpos+1):(xpos+psize),(ypos+1):(ypos+psize2),3)=im;
   
   if ismember(ppi(ii),hiliteidx),
      % highlight in red
      cc=[255 0 0];
      r0=-1:0;
      r1=0:1;
   else
      cc=[0 0 255];
      r0=0; r1=0;
   end
   
   bigplot((xpos+1):(xpos+psize),(ypos+1)+r0,1)=cc(1);
   bigplot((xpos+1):(xpos+psize),(ypos+psize2+r1),1)=cc(1);
   bigplot((xpos+1+r0),(ypos+1):(ypos+psize2),1)=cc(1);
   bigplot((xpos+psize+r1),(ypos+1):(ypos+psize2),1)=cc(1);
   bigplot((xpos+1):(xpos+psize),(ypos+1+r0),2)=cc(2);
   bigplot((xpos+1):(xpos+psize),ypos+psize2+r1,2)=cc(2);
   bigplot((xpos+1+r0),(ypos+1):(ypos+psize2),2)=cc(2);
   bigplot((xpos+psize+r1),(ypos+1):(ypos+psize2),2)=cc(2);
   bigplot((xpos+1):(xpos+psize),(ypos+1+r0),3)=cc(3);
   bigplot((xpos+1):(xpos+psize),(ypos+psize2+r1),3)=cc(3);
   bigplot((xpos+1+r0),(ypos+1):(ypos+psize2),3)=cc(3);
   bigplot((xpos+psize+r1),(ypos+1):(ypos+psize2),3)=cc(3);
end
