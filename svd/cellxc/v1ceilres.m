% function v1cellcomp(cellid/runidx,batch,predstart,predstop,outnlidx)
%
% load results from resfile in sRunData and display fits, pred
% results, given cellid and a range of batch ids.  optional
% predstart and predstop are first and last frames of prediction
% psth comparison to show
%
%
function v1ceilres(cellid,batch,outnlidx)

rcsetstrings;
dbopen;

sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
     ' AND batch=',num2str(batch)];

trd=mysql(sql);
if isempty(trd),
   disp('no entries found in db!');
   return
end

rundata=trd;
fprintf('loading %s\n',[rundata.respath,rundata.resfile,'.gz']);
zload([rundata.respath,rundata.resfile,'.gz']);


repcount=0;
onetrialresp=[];

for fidx=1:filecount,
   
   if cellfiledata(fidx).repcount>repcount,
      cnfidx=fidx;
      repcount=cellfiledata(fidx).repcount;
   end
   
   if ~isnan(stoptimes(fidx,attidx)),
      stimstart(fidx)=max([starttimes(fidx,attidx)-params.maxlag(2) 1]);
      
      framecount=imfileinfo(params.stimfiles{fidx},1);
      stimstop(fidx)=min([stoptimes(fidx,attidx)-params.maxlag(1) framecount]);
   else
      stimstop(fidx)=0;
   end
   
   fprintf('%d: %d-->%d now %d-->%d\n',fidx,starttimes(fidx,attidx),...
           stoptimes(fidx,attidx),stimstart(fidx),stimstop(fidx));
   
   bslen=imfileinfo(params.stimfiles{fidx});
   if bslen<stimstop(fidx)-stimstart(fidx)+1,
      stimstop(fidx)=bslen+stimstart(fidx)-1;
   end
   
   tresp=respload(params.respfiles{fidx},cellfiledata(fidx).respvarname, ...
                        cellfiledata(fidx).respfiletype,1,0);
   tresp=compact_raster_matrix3(tresp(:,2:end));
   
   
   rsize=size(tresp);
   
   % trim response to appropriate time range
   if ~isnan(stoptimes(fidx,attidx)),
      if rsize(1)<stimstop(fidx),
         tresp=cat(1,tresp(stimstart(fidx):end,:,:),...
                   ones([stimstop(fidx)-rsize(1) rsize(2:end)]));
      else
         tresp=tresp(stimstart(fidx):stimstop(fidx),:,:);
      end
   end
   
   respcount=size(tresp,2);
   tresp(1:params.maxlag(2),:)=nan;
   tresp(length(tresp)+params.maxlag(1):end,:)=nan;
   
   if fidx==cnfidx,
      cnfresp=tresp;
      cnfpred=predmtx(length(onetrialresp)+(1:length(tresp)),:);
   end
   onetrialresp=cat(1,onetrialresp,tresp(:,1));
end


nlidx=params.nlidxsave;
disp('overriding nlidxsave and using nlidx=3');
nlidx=3;

figure(1);

tsf=cat(3,strf(nlidx,1,:).h);
kernfmt=strf(nlidx,1,1).parms.kernfmt;
iconside=strf(nlidx,1,1).parms.iconside;
titles={sprintf('%s 10%%',strf(nlidx,1,1).name),...
        sprintf('%s 20%%',strf(nlidx,1,1).name),...
        sprintf('%s 60%%',strf(nlidx,1,1).name)};
showkern(tsf,kernfmt,iconside,titles);

disp('computing one-trial predxc');

if ~isfield(params,'bfracs'),
   params.bfracs=bfracs;
end

segcount=length(params.bfracs)-1;
jcount=20;
cxy=zeros(jcount,segcount);
for nn=1:size(predmtx,2),
   goodidx=find(~isnan(onetrialresp) & ~isnan(predmtx(:,nn)));
   
   jstep=length(goodidx)/jcount;
   for jackn=1:jcount,
      jidx=goodidx([1:round((jackn-1)*jstep) ...
                    round((jackn+1)*jstep):length(goodidx)]);
      
      pp=predmtx(jidx,nn);
      rr=onetrialresp(jidx);
      cxy(jackn,nn)=xcov(pp,rr,0,'coeff');
   end
end
cxy=cxy.*abs(cxy); % preserve sign in R^2
exy=std(cxy,1).*sqrt(jcount-1);

y=cxy;
x=params.bfracs(1:segcount);
[yinf,yinferr,p]=fitceiling(x,y);

if size(cnfresp,2)>1,
   [ccr,ccrM,ccre]=frednoise(cnfresp);
   fprintf('cnfidx=%d. reps=%d. fred noise factor: ccr=%.2f\n',...
           cnfidx,repcount,ccr);
   [ccr2]=valnoise(cnfpred(:,3),cnfresp);
else
   ccr=1;
   ccr2=1;
   fprintf('cnfidx=%d. reps=%d. no reps for fnf. fix ccr=1\n',...
           cnfidx,repcount);
end

figure(2);
clf
subplot(3,1,1);
hold on
x0=linspace(0,1,20);
for ii=1:jcount,
   if ~isnan(p(ii,2)),
      y0=x0./(p(ii,2).*x0+p(ii,1));
      plot(x0,y0,'k:');
   end
end
errorbar(x,mean(cxy),exy);
hold off

title(sprintf('%s yinf=%.3f',cellid,yinf));

subplot(3,1,2);
hold on
x0=linspace(0,1,20);
for ii=1:jcount,
   if ~isnan(p(ii,2)),
      y0=x0./(p(ii,2).*x0+p(ii,1))./ccr;
      plot(x0,y0,'k:');
   end
end
errorbar(x,mean(cxy)./ccr,exy);
hold off

title(sprintf('%s yinf+fredcorr=%.3f',cellid,yinf./ccr));

subplot(3,1,3);
hold on
x0=linspace(0,1,20);
for ii=1:jcount,
   if ~isnan(p(ii,2)),
      y0=x0./(p(ii,2).*x0+p(ii,1))./ccr2;
      plot(x0,y0,'k:');
   end
end
errorbar(x,mean(cxy)./ccr2,exy);
hold off

title(sprintf('%s yinf+cnfasymp=%.3f',cellid,yinf./ccr2));

%keyboard



return


%
% display kernels from each batch
%
figure(1);
clf

for ii=find(goodbatch(:)'),
   tsf=ones([size(z{ii}.strf(1).h(:,:)),batchcount])*nan;
   if ~isfield(z{ii}.params,'batch'),
      z{ii}.params.batch=z{ii}.params.id;
   end
   if 0 & ii==1,
      tsf(:,:,1)=z{ii}.mH(:,(1:size(z{ii}.hf,2))-z{ii}.maxlag(1),1);
      titles{1}=sprintf('%s: Raw STA (%s)',cellid,z{1}.skernfile);
   end
   tsf(:,:,ii)=z{ii}.strf(outnlidx2(ii)).h;
   
   cnfstr='';
   for cnfidx=1:length(z{ii}.params.predbatch),
      cnfstr=sprintf('%s %d: %.2f',cnfstr,...
                     z{ii}.params.predbatch{cnfidx},...
                     z{ii}.predxc(cnfidx,outnlidx(ii),1,1));
   end
   
   titles{ii}=sprintf('%s Bat %d/nl %d: STRF (%s) Pred:%s',...
                      cellid,batch(ii),outnlidx(ii),...
                      basename(z{ii}.params.outfile),cnfstr);
   
   showkern(tsf,z{ii}.params.kernfmt,z{ii}.params.iconside,titles);
end
set(gcf,'PaperOrientation','landscape','PaperPosition',[0.25 0.25 10.5 8]);

%
% space-time decomposition
%
figure(2)
clf

if isfield(z{1}.strf(outnlidx2(1),1),'powunbiased'),
   commonpowunbiased=z{1}.strf(outnlidx2(1),1).powunbiased;
   if max(commonpowunbiased)>1,
      commonpowunbiased=commonpowunbiased./max(commonpowunbiased);
   end
   
else
   commonpowunbiased=ones(size(z{1}.strf(outnlidx2(1),1).h,1),1);
end

for ii=find(goodbatch(:)'),
   
   if strcmp(z{ii}.params.stimfiltercmd,'movphasesep'),
      kernfmt='pfft+4';
   else
      kernfmt=z{ii}.params.kernfmt;
   end
   
   % plot summary kernel stuff - take svd of kernel to get spatial
   % and temporal marginals
   h=z{ii}.strf(outnlidx2(ii)).h;
   spacecount=size(h,1);
   tcount=size(h,2);
   
   if isfield(z{ii}.strf(outnlidx2(ii),1),'powunbiased'),
      powunbiased=z{ii}.strf(outnlidx2(ii)).powunbiased;
      if max(powunbiased)>1,
         powunbiased=powunbiased./max(powunbiased);
      end
   else
      powunbiased=ones(size(h,1),1);
   end
   
   if length(commonpowunbiased)==length(powunbiased),
      ratub=commonpowunbiased./powunbiased;
      ratub(find(ratub>1))=1;
      ratub(find(ratub<0))=0;
   else
      ratub=ones(size(powunbiased));
   end
   
   hub=h .* repmat(ratub,1,tcount);
   
   if isfield(z{ii}.strf(outnlidx2(ii),1),'hspace'),
      spvec=z{ii}.strf(outnlidx2(ii),1).hspace;
      spvecub=z{ii}.strf(outnlidx2(ii),1).hspacebiased;
      ttimesep=z{ii}.strf(outnlidx2(ii),1).tempresp;
      
   else
      ttime=sum(h,1)';
      ttimepos=sum(h.*(h>0),1)';
      ttimeneg=sum(h.*(h<0),1)';
      
      cumttime=cumsum(ttime);
      tmax=min(find(cumttime==max(cumttime)));
      if tmax==1,
         tmax=min(find(ttimepos==max(ttimepos)));
      end
      
      fprintf('tmax=%d\n',tmax);
      tspacemaxpos=sum(h(:,1:tmax),2);
      tspaceall=sum(h,2);
      
      [u,s,v]=svd(hub);
      if tspacemaxpos'*u(:,1) > 0;
         ttimesep=v(:,1);
         tspacesep=u(:,1);
      else
      ttimesep=-v(:,1);
      tspacesep=-u(:,1);
      end
      teigvals=diag(s);
      
      spvec=tspacemaxpos;
      spvecub=spvec .* ratub;
   end
   
   ssp=sort(spvec);
   sptopmean=mean(ssp(end-3:end));
   
   tH=ones(spacecount,4,length(goodbatch)).*nan;
   tH(:,1,ii)=spvec;
   tH(:,2,ii)=powunbiased./max(powunbiased).*sptopmean;
   tH(:,3,ii)=spvecub./max(spvecub).*sptopmean;
   tH(:,4,ii)=spvecub./max(spvecub).*sptopmean;
   
   showkern(tH,kernfmt,z{ii}.iconside,titles);
   
   tvec=ttimesep;
   %tvec=tvec-mean(tvec([1 end]));
   tvec=tvec./max(abs(tvec));
   
   if length(tvec)>1,
      subplot(length(goodbatch),4,ii*4);
      
      tt=(1:length(tvec))*14-14;
      h=plot(tt,tvec,'k-','LineWidth',2);
      hold on
      plot(tt,zeros(1,length(tvec)),'k--');
      hold off
      
      if max(abs(tvec))>0,
         axis([0 max(tt) -max(abs(tvec)*1.1) max(abs(tvec)*1.1)]);
      end         
      axis square
   end
end

fullpage('portrait');
drawnow

%
% example prediction vs actual PSTH
%

% skip if not review
if z{1}.params.runclassid~=0,
   disp('not review. skipping pred plot');
   return
end

checkpath='/auto/k2/share/data/movfixdat/';
if strcmp(cellid(1:3),'93G'),
   movThresh=10;
else
   movThresh=0.75;
end

cnfidx=find(cat(1,z{1}.params.predbatch{:})==z{1}.params.batch);
if isempty(cnfidx),
   cnfidx=find(cat(1,z{1}.params.predbatch{:})==z{1}.params.batch);
end
if isempty(cnfidx),
   cnfidx=1;
end
tpred=z{1}.predres(cnfidx).mod_psth{1}(:,1,outnlidx2(1));
tact=z{1}.predres(cnfidx).act_resp{1}(:,1);
bidx=find(~isnan(tpred) & ~isnan(tact));
firstkeep=bidx(1);
tpred=tpred(bidx);
tact=tact(bidx);

if ~exist('predstart','var') | predstart<1,
   predstart=1;
end
if ~exist('predstop','var') | predstop<1 | predstop>length(tpred),
   predstop=length(tpred);
end


[pcellfiledata,ptimes,pbatchdata]=...
    cellfiletimes(z{1}.params.cellid,z{1}.params.predbatch{cnfidx});

stimfile=[pcellfiledata(ptimes(3).fileidx).stimpath, ...
          pcellfiledata(ptimes(3).fileidx).stimfile];
[FirstFr,LastFr]=movgetcheck(checkpath,cellid,stimfile,movThresh);

sidx=ptimes(3).start-z{1}.params.maxlag(2);
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


figure(3)
clf

subplot(batchcount+1,1,1);
imagesc(bigimagesc);

hold on
for ii=1:size(mov,3),
   plot(sframes(ii).*[1 1],[size(mov,1) size(mov,1).*1.1],'k-');
end
hold off
colormap(gray);
axis image
axis off
title(sprintf('cell %s validation stimulus',z{1}.params.cellid));

for ii=find(goodbatch(:)'),
   
   cnfidx=find(cat(1,z{ii}.params.predbatch{:})==z{1}.params.batch);
   if isempty(cnfidx),
      cnfidx=find(cat(1,z{ii}.params.predbatch{:})==z{ii}.params.batch);
   end
   if isempty(cnfidx),
      cnfidx=1;
   end
   
   if ~isempty(cnfidx),
      
      % plot predictions
      hs=subplot(batchcount+1,1,ii+1);
      tpred=z{ii}.predres(cnfidx).mod_psth{1}(:,1,outnlidx2(ii));
      tact=z{ii}.predres(cnfidx).act_resp{1}(:,1);
      if length(tact)~=length(tpred),
         tact=tact(z{ii}.leadbincount+1:end);
         tact=tact(find(~isnan(tact)));
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
      
      % last used:
      %data_predict_plot(tact', tpred');
      cla;
      plot(tt,tact,'b--','Linewidth',1);
      hold on
      plot(tt,tpred,'r-','Linewidth',1);
      
      a=axis;
      
      for jj=1:length(FirstFr),
         plot((FirstFr(jj)-sidx).*[14 14],[0 a(4)],'k:');
         %plot((LastFr(jj)-sidx).*[14 14],[0 a(4)],'r--');
      end
      
      hold off
      a=axis;
      axis([tt(1) tt(end) -a(4)*0.05 a(4)]);
      if ii<max(find(goodbatch(:))),
         set(gca,'XTickLabel',[]);
      end
      
      title(sprintf('%s/%d/%d: pred batch %d (r=%.3f)',...
                    cellid,batch(ii),outnlidx(ii),...
                    z{ii}.params.predbatch{cnfidx},...
                    z{ii}.predxc(cnfidx,outnlidx(ii))));
   end
   
   if 0,
      figure(3);
      subplot(length(goodbatch),1,ii);
      [cxy,f]=cohere(tact,tpred,64,72,[],32);
      plot(f,cxy);
      figure(2);
   end
end
legend('act','pred')

fullpage('landscape');



return

   figure(3);
   
   if z{1}.respfmtcode==1,
      for ii=1:(size(z{1}.mH,2)*size(z{1}.mH,3)),
         z{1}.mH(:,ii)=z{1}.mH(:,ii)./mean(abs(z{1}.mH(:,ii)));
      end
   end
   
   sampcount=6;
   sampidx=[round(linspace(1,z{1}.sfscount,sampcount))];
   smm=z{1}.mH(:, (1:size(z{1}.hf,2))-z{1}.maxlag(1), sampidx);
   sms=z{1}.eH(:, (1:size(z{1}.hf,2))-z{1}.maxlag(1), sampidx) .* ...
       z{1}.sigrange(1);
   smd=abs(smm)./(sms+(sms==0));
   if z{1}.DOTHRESH, % old--shrinkage filter
      % new -- threshold by # of std errs
      tsf=smm.*(smd>1);
   else
      smd=(1-smd.^(-2));
      smd=smd.*(smd>0);
      smd(find(isnan(smd)))=0;
      tsf=smm.*smd;
   end
   
   for ii=1:sampcount,
      
      if isfield(z{1},'lambda'),
         tl=z{1}.lambda(sampidx(ii));
      else
         tl=ii;
      end
      titles{ii}=sprintf('%s: STRF sample %.2f xc L/NL=%0.3f/%0.3f',...
                         cellid,tl,z{1}.xc(sampidx(ii),1,1),...
                         z{1}.xc(sampidx(ii),1,3));
   end
   showkern(tsf,z{1}.kernfmt,z{1}.iconside,titles);
   %showkern(squeeze(z{ii}.mSR),z{ii}.kernfmt,z.iconside)
   %showkern(squeeze(mH(:,:,end,:,:)),kernfmt,iconside)
   
   
   if z{1}.respfmtcode==1 & z{1}.attcount>1,
      figure(4);
      clf
      minx=min(min(min(z{1}.predxc(:,:,:,1))));
      maxx=max(max(max(z{1}.predxc(:,:,:,1))));
      if minx<-maxx & maxx>0,
         minx=-maxx;
      end
      
      subplot(ceil(z{1}.latcount/3)+1,2,1);
      plot(nanmean(z{1}.resp(:,:,1)));
      title([z{1}.cellid,' mean resp']);
      subplot(ceil(z{1}.latcount/3)+1,2,2);
      plot(squeeze(z{1}.predxc(1,1,:,1)));
      title('full xc');
      
      for ii=1:z{1}.latcount,
         subplot(ceil(z{1}.latcount/3)+1,3,ii+3);
         if sum(sum(abs(z{1}.predxc(:,:,ii)-mean(mean(z{1}.predxc(:,:,ii))))))>0,
            imagesc(z{1}.predxc(:,:,ii),[minx maxx]);
         end
         axis image
         axis off
         title(sprintf('lat=%d',ii));
         
         if ii==z{1}.latcount,
            colorbar
         end
      end
      colormap(hot);
      
      if length(z{1}.respfilterparms)>0 & ~strcmp(z{1}.respfiltercmd,''),
         
         sb=unique(z{1}.respfilterparms{1});
         eb=unique(z{1}.respfilterparms{2});
         
         figure(5);
         clf
         for attidx=1:z{1}.attcount,
            txc=zeros(length(eb),length(sb));
            for latidx=1:z{1}.latcount,
               sbidx=find(z{1}.respfilterparms{1}(latidx)==sb);
               ebidx=find(z{1}.respfilterparms{2}(latidx)==eb);
               txc(ebidx,sbidx)=z{1}.predxc(attidx,attidx,latidx);
            end
            subplot(1,z{1}.attcount,attidx);
            if maxx-minx > 0,
               imagesc(eb,sb,txc,[minx maxx]);
            end
            
            if attidx==1,
               ylabel('bin end time');
            end
            xlabel('bin start time');
            title(sprintf('%s att=%d',z{1}.cellid,attidx));
         end
         
         colormap(hot);
      end
   else
      for ii=1:batchcount,
         predcount=size(z{ii}.predxc,1);
         fprintf('File: %s\n',z{ii}.skernfile);
         for nlidx=1:z{ii}.nlcount,
            fprintf('nl=%d: (sfs,sig)=(%d,%d)',...
                    nlidx,z{ii}.sfsfit(nlidx),z{ii}.sigfit(nlidx));
            for predidx=1:predcount,
               fprintf('%7.3f',z{ii}.predxc(predidx,:,:,nlidx));
            end
            fprintf('\n');
         end
      end
   end


