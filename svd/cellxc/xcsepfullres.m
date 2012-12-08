% function r=xcsepfullres(cellid,batch)
%
% load results from resfile in sRunData and display fits, pred
% results
%
% r=0 if no entries found in db, =1 otherwise
%
function r=xcsepfullres(cellid,batch)

dbopen;
goodbatch=zeros(1,length(batch));
sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
     ' AND batch=',num2str(batch)];
rundata=mysql(sql);

if length(rundata)==0,
   disp('no entry found in db!');
   if nargout>0,
      r=0;
   end
   return
end

resfile=[rundata.respath,rundata.resfile,'.gz'];
fprintf('%s: loading %s\n',mfilename,resfile);
zload(resfile);

global BATQUEUEID
BATQUEUEID=[];

fitcount=size(strf,1);
bootcount=size(strf,2);
segcount=size(strf,3);

SHOWTEMPDECORR=0;
SHOWCEILING=1;

SHOWSEGS=0;
SHOWBOOTS=0;
SHOWFITS=0;   
if segcount>1,
   SHOWSEGS=1;
   strfcount=segcount;
   disp('showing different seg strfs');
elseif bootcount>1,
   SHOWBOOTS=1;
   strfcount=bootcount;
   disp('showing different jackknife strfs');
else
   SHOWFITS=1;
   strfcount=fitcount;
   disp('showing different stages of fit');
end

hshow=[];
tpredxc=zeros(strfcount,1);
tstrfid=zeros(strfcount,1);
for bidx=1:strfcount,
   if SHOWSEGS,
      tstrf=strf(end,end,bidx);
      tpredxc(bidx)=predxc(bidx,end);
      tstrfid(bidx)=sub2ind([fitcount bootcount segcount],...
                            fitcount,bootcount,bidx);
   elseif SHOWBOOTS,
      tstrf=strf(end,bidx);
      tpredxc(bidx)=vpredxc(end,bidx,end);
      tstrfid(bidx)=sub2ind([fitcount bootcount segcount],...
                            fitcount,bidx,1);
   else
      tstrf=strf(bidx);
      tpredxc(bidx)=predxc(1,bidx);
      tstrfid(bidx)=sub2ind([fitcount bootcount segcount],...
                            bidx,1,1);
   end
   
   if ~isfield(tstrf,'hspace'),
      hspace=sum(tstrf.h,2);
   else
      hspace=tstrf.hspace;
   end
   
   if isfield(tstrf,'hspacebiased'),
      hsb=tstrf.hspacebiased;
   else
      hsb=[];
   end
   hsb=hsb./max(abs(hsb)).*max(abs(hspace));
   
   hshow=cat(3,hshow,[hspace hsb hspace]);
end

figure(1);
clf
showkern(hshow,strf(1).parms.kernfmt);
tbins=(1:length(strf(1).tempresp))*14-7;
colcount=size(hshow,2);
for bidx=1:strfcount,
   subplot(strfcount,colcount,(bidx-1)*colcount+2);
   title(sprintf('%s bidx=%d: predxc=%.3f',...
                 cellid,bidx,tpredxc(bidx)));
   
   subplot(strfcount,colcount,bidx*colcount);
   tr1=strf(tstrfid(bidx)).tempresp;
   plot(tbins,tr1,'-');
   
   title('temp resp');
   
   fprintf('bidx=%d. sfs=%d sig=%d predxc=%.3f\n',...
           bidx,strf(tstrfid(bidx)).parms.sfsfit,...
           strf(tstrfid(bidx)).parms.sigfit,tpredxc(bidx));
end
xlabel('time lag (ms)');

set(gcf,'PaperOrientation','Portrait',...
        'PaperPosition',[1.5 0.25 6 10.5]);


%
% plot range of decorrelated kernel possibilities
%

figure(2);
clf

if ~exist('sigcount','var'),
   sigcount=params.sffiltsigma;
end
if ~exist('sigrange','var'),
   sigrange=exp(linspace(log(0.9),log(1.8),sigcount));
end

strfidx=length(tstrfid); % length(tstrfid)-1;

tH=repmat(strf(tstrfid(strfidx)).mH,[1 sigcount]);
eH=strf(tstrfid(strfidx)).eH;
for ii=1:sigcount,
   tH(:,ii,:)=shrinkage(tH(:,ii,:),eH,sigrange(ii));   
end
tH=permute(tH,[1 3 2]);
for ii=1:size(tH,2)*size(tH,3),
   if max(abs(tH(:,ii)))>0,
      tH(:,ii)=tH(:,ii)./max(abs(tH(:,ii)));
   end
end
for sigidx=1:sigcount,
   titles{sigidx}=sprintf('cellid %s sig=%.1f   (used %d/%d)',...
                          cellid,sigrange(sigidx),...
                          strf(tstrfid(strfidx)).parms.sfsfit,...
                          strf(tstrfid(strfidx)).parms.sigfit);
end

showkern(tH,strf(tstrfid(strfidx)).parms.kernfmt,...
         strf(tstrfid(strfidx)).parms.iconside,titles);
if isfield(strf(tstrfid(strfidx)),'pctvar'),
   for cc=1:size(tH,2),
      subplot(size(tH,3),size(tH,2),(size(tH,3)-1)*size(tH,2)+cc);
      xlabel(sprintf('%.4f',strf(tstrfid(strfidx)).pctvar(cc)));
   end
end

set(gcf,'PaperOrientation','Portrait',...
        'PaperPosition',[1.5 0.25 6 6]);


figure(3);
clf

imagesc(xc(:,:,end));
colormap(hot);
hold on
plot(strf(tstrfid(end)).parms.sigfit,strf(tstrfid(strfidx)).parms.sfsfit,'kx');
hold off
title('sig/sfs fit for optspacebidx kernel');

colorbar


%predstart=1101;
%predstop=1600;
%fprintf('predstart-predstop: %d-%d\n',predstart,predstop);

if params.predfrac==0,
   disp('no predfrac, skipping psth plot');
   return
end

% skip if not review
if params.runclassid==0,
   
   checkpath='/auto/k2/share/data/movfixdat/';
   if strcmp(cellid(1:3),'93G'),
      movThresh=10;
   else
      movThresh=0.75;
   end
   
   cnfidx=find(cat(1,params.predbatch{:})==params.batch);
   if isempty(cnfidx),
      cnfidx=1;
   end
   tpred=predres(cnfidx).mod_psth{1}(:,1,tstrfid(1));
   tact=predres(cnfidx).act_resp{1}(:,1);
   bidx=find(~isnan(tpred) & ~isnan(tact));
   firstkeep=bidx(1);
   tpred=tpred(bidx);
   tact=tact(bidx);
   
   if ~exist('predstart','var') | predstart<1 | predstart>length(tpred),
      predstart=1;
   end
   if ~exist('predstop','var') | predstop<1 | predstop>length(tpred),
      predstop=length(tpred);
   end
   
   [pcellfiledata,ptimes,pbatchdata]=...
       cellfiletimes(params.cellid,params.predbatch{cnfidx});
   
   stimfile=[pcellfiledata(ptimes(3).fileidx).stimpath, ...
             pcellfiledata(ptimes(3).fileidx).stimfile];
   [FirstFr,LastFr]=movgetcheck(checkpath,cellid,stimfile,movThresh);
   
   sidx=ptimes(3).start-params.maxlag(2);
   if sidx<1,
      sidx=1;
   end
   sidx=sidx+firstkeep-1+predstart;
   eidx=sidx-predstart+predstop;
   goodfr=find(FirstFr>=sidx & LastFr<=eidx);
   FirstFr=FirstFr(goodfr);
   LastFr=LastFr(goodfr);
   
   mov=loadimframes(stimfile,FirstFr);
   
   % normalize each frame for cleaner display
   mbase=mov(1,1,1);
   mov=mov-mbase;
   mmax=255-mbase;
   mmin=-mbase;
   for ii=1:size(mov,3),
      tmax=max(max(mov(:,:,ii).*(mov(:,:,ii)>0)));
      tmin=min(min(mov(:,:,ii).*(mov(:,:,ii)<0)));
      if (tmax/mmax)<(tmin/mmin) & tmin<0,
         mov(:,:,ii)=round(mov(:,:,ii).*(mmin/tmin));
      else
         mov(:,:,ii)=round(mov(:,:,ii).*(mmax/tmax));
      end
   end
   mov=mov+mbase;
   
   minfixlen=min(LastFr-FirstFr+1);
   imagescale=size(mov,2)/minfixlen;
   sframes=round((FirstFr-sidx)*imagescale);
   bigimagesc=ones(size(mov,1),round((eidx-sidx+1)*imagescale)).*255;
   for ii=1:size(mov,3),
      bigimagesc(1:size(mov,1),sframes(ii):sframes(ii)+size(mov,2)-1)=...
          fliplr(flipud(mov(:,:,ii)));
   end
end

figure(4)
clf

subplot(strfcount+1,1,1);
if params.runclassid==0,
   imagesc(bigimagesc);
   
   hold on
   for ii=1:size(mov,3),
      plot(sframes(ii).*[1 1],[size(mov,1) size(mov,1).*1.1],'k-');
   end
   hold off
   colormap(gray);
   axis image
   axis off
   title(sprintf('cell %s validation stimulus',params.cellid));
end

for ii=1:strfcount,
   cnfidx=tstrfid(ii);
         
   % plot predictions
   hs=subplot(strfcount+1,1,ii+1);
   tpred=predres(1).mod_psth{1}(:,1,tstrfid(ii));
   
   tact=predres(1).act_resp{1}(:,1);
      
   if length(tact)~=length(tpred),
      disp('length(tact)~=length(tpred) ??');
      keyboard
   end
   
   bidx=find(~isnan(tpred) & ~isnan(tact));
   tpred=tpred(bidx);
   tact=tact(bidx);
      
   if ~exist('predstart','var') | predstart<1,
      predstart=1;
   end
   if ~exist('predstop','var') | predstop>length(tpred) | predstop<1,
      predstop=length(tpred);
   end
   
   tpred=tpred(predstart:predstop);
   tpred=tpred-nanmin(tpred);
   
   tact=tact(predstart:predstop);
   
   %tbinms=z{ii}.strf(1).parms.tbinms;
   tbinms=14;
   tt=0:tbinms:(length(tact)-1)*tbinms;
   
   %data_predict_plot(data, prediction);
   %data_predict_plot(tact'./std(tact), tpred'./std(tpred));
   %data_predict_plot(tact', tpred'./sqrt(mean(tpred.^2)/mean(tact.^2)));
   
   cla;
   plot(tt,tact,'b--','Linewidth',1);
   hold on
   plot(tt,tpred,'r-','Linewidth',1);
   
   a=axis;
   
   if params.runclassid==0,
      
      for jj=1:length(FirstFr),
         plot((FirstFr(jj)-sidx).*[14 14],[0 a(4)],'k:');
      end
   end
   
   hold off
   a=axis;
   axis([tt(1) tt(end) -a(4)*0.05 a(4)]);
   if ii<strfcount,
      set(gca,'XTickLabel',[]);
   end
   
   title(sprintf('%s/%d: tstrfid=%d (r=%.3f)',...
                 cellid,ii,cnfidx,tpredxc(ii)));
end

if 0,
   figure(3);
   subplot(length(goodbatch),1,ii);
   [cxy,f]=cohere(tact,tpred,64,72,[],32);
   plot(f,cxy);
   figure(2);
end

legend('act','pred')
fullpage('portrait');


if SHOWTEMPDECORR,
   figure(5);
   clf
   
   mtime=squeeze(strf(2,end,end).mH);
   plot(0:params.maxlag(2),mtime(-params.maxlag(1)+1:end,1),'k--');
   hold on
   plot(0:params.maxlag(2),mtime(-params.maxlag(1)+1:end,2:end-1),'k:');
   plot(0:params.maxlag(2),mtime(-params.maxlag(1)+1:end,end),'k-');
   %plot(0:params.maxlag(2),strf(2,1,end).tempresp,'b-');
   hold off
   
   title('temporal response function vs. sfs');
   xlabel('time lag (bins)');
   drawnow
end

%keyboard

if SHOWCEILING,
   figure(6);
   clf
   
   ro=respload(params.respfiles{tpredfile});
   
   tpredstartframe=max([times(3).start-params.maxlag(2) 1]);
   tpredstopframe=min([times(3).stop-params.maxlag(1) size(ro,1)]);
   tpredfile=times(3).fileidx;
   
   ro=ro(tpredstartframe:tpredstopframe,:);
   r=compact_raster_matrix3(ro(:,2:end));      
   
   if size(r,2)>4,
      disp('measuring validation noise...');
      strfidx=sub2ind([fitcount bootcount segcount],...
                      ones(1,bootcount).*fitcount,1:bootcount,...
                      ones(1,bootcount).*segcount);
      if 0,
         p=mean(predres.mod_psth{1}(:,strfidx),2);
         if length(p)>size(r,1),
            p=p(1:size(r,1));
         end
         [ccr,ccrM,yinf,n,ty,ttp]=valnoise(p,r);
         ttp=ttp(find(~isnan(ttp(:,1))),:);
         
         if exist('vpredxc','var'),
            pxc=mean(squeeze(vpredxc(end,end,:)).^2);
         else
            pxc=predxc(end,end).^2;
         end
         adjrat=(mean(ty(:,end))./pxc);
         ttp=ttp.*adjrat;
         yinf=mean(1./ttp(:,2));
         ty=ty./adjrat;
         ccrM=pxc./yinf;
         
      yinf=1./ttp(:,2);
      else
         ccrM=zeros(bootcount,1);
         yinf=zeros(bootcount,1);
         ty=[];
         tp=[];
         for bb=1:bootcount,
            p=predres.mod_psth{1}(:,strfidx(bb));
            if length(p)>size(r,1),
               p=p(1:size(r,1));
            end
            [ccr,ccrM(bb),yinf(bb),n,tty,ttp]=valnoise(p,r);
            ttp=ttp(find(~isnan(ttp(:,1)) & ttp(:,2)~=0),:);
            if exist('vpredxc','var'),
               pxc=mean(squeeze(vpredxc(end,end,:)).^2);
            else
            pxc=predxc(end,end).^2;
            end
            adjrat=(mean(tty(:,end))./pxc);
            ttp=ttp.*adjrat;
            yinf(bb)=mean(1./ttp(:,2));
            tty=tty./adjrat;
            ccrM(bb)=pxc./yinf(bb);
         
            ty=[ty; tty];
            tp=[tp;ttp];
         end
         
         yinf=1./ttp(:,2);
         ty=tty;
      end
      
      subplot(1,2,1);
      cla
      plot(0,0);
      aa=1:(size(r,2)+4);
      bb=1./(ttp(:,2)*ones(size(aa))+ttp(:,1)*(1./aa));
      errorshade(aa,mean(bb),std(bb),[0 0 0],[0.9 0.9 0.9]);
      hold on
      
      errorbar(n,mean(ty,1),std(ty,0,1).*sqrt(size(ty,1)-1),'k+');
      errorbar(size(r,2)+6,mean(yinf),std(yinf).*sqrt(length(yinf)-1),'ko');
      hold off
      title('val noise');
      xlabel('validation trials');
   
      ccrM=nanmedian(ccrM);
      if ccrM<0.01,
         ccrM=nan;
      end
   else
      subplot(1,2,1);
      cla
      title('not enough data for val noise');
      ccrM=nan;
   end
   
   if segcount>1,
      disp('measuring estimation noise...');
      pxc=squeeze(vpredxc(:,end,:))';
      x=params.bfracs(1:size(pxc,2)).*sum(~isnan(resp(:,1)));
      y=pxc.^2; %.*abs(pxc);
      
      tp=zeros(size(y,1),2);
      xcestinf=zeros(size(y,1),1);
      for ii=1:size(y,1),
         [xcestinf(ii),errestinf,tp(ii,:)]=...
             fitceiling(x,mean(y([1:ii-1 ii+1:end],:)));
         %[xcestinf(ii),errestinf,tp(ii,:)]=fitceiling(x,y(ii,:));
      end
      mxcestinf=mean(xcestinf);
      
      %[xcestinf,errestinf,tp]=fitceiling(x,mean(y));
      
      if isnan(ccrM),
         ccrM0=1;
      else
         ccrM0=ccrM;
      end
      
      subplot(1,2,2);
      cla
      plot(0,0);
      aa=linspace(0.05,max(x).*1.8);
      bb=1./(tp(:,2)*ones(size(aa))+tp(:,1)*(1./aa)) ./ ccrM0;
      
      errorshade(aa,mean(bb),std(bb).*sqrt(size(y,1)-1).*2,...
                 [0 0 0],[0.9 0.9 0.9]);
      hold on
   
      errorbar(x,mean(y) ./ ccrM0,std(y ./ ccrM0).*2,'k+');
      errorbar(max(x).*2.0,mean(xcestinf ./ ccrM0),...
               std(xcestinf ./ ccrM0).*sqrt(size(y,1)).*2,'ko');
      hold off
      title('est noise (val corr)');
      xlabel('est samples');
      ylabel('r2valmax');
      a=axis;
      axis([0 max(x).*2.4 a(3) a(4)]);
      axis square
      
      subplot(1,2,1);
      axis([0 size(r,2)+10 a(3) a(4)]);
      axis square
      
      fprintf('samples:  ');
      fprintf('%7d',round(x));
      fprintf('\npreds:    ');
      fprintf('%7.3f',sqrt([ mean(pxc.^2) mxcestinf mxcestinf./mean(ccrM)]));
      fprintf('\npreds2:   ');
      fprintf('%7.3f',[ mean(pxc.^2) mxcestinf mxcestinf./mean(ccrM)]);
      fprintf('\n');
      
      fid=fopen('/auto/k1/david/tmp/res.dat','a');
      fprintf(fid,' %10.4f',cellfiledata(1).masterid,params.batch);
      fprintf(fid,' %10.4f',round(x));
      fprintf(fid,'%10.4f',([ mean(pxc.^2) mxcestinf mxcestinf./mean(ccrM)]));
      fprintf(fid,'\n');
      fclose(fid);
   else
      pxc=predxc(:,end);
      xcestinf=nan;
      
      subplot(1,2,2);
      cla
      title('no est noise ceil');
   end
end

%[ mean(pxc.^2) xcestinf xcestinf./mean(ccrM)]

if nargout==0,
   clear r
end

