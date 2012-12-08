% function v1cellcomp(cellid/runidx,batch,predstart,predstop,outnlidx)
%
% load results from resfile in sRunData and display fits, pred
% results, given cellid and a range of batch ids.  optional
% predstart and predstop are first and last frames of prediction
% psth comparison to show
%
%
function v1cellcomp(runidx,batch,predstart,predstop,outnlidx,cnfshowidx)

if ~exist('outnlidx','var'),
   outnlidx=ones(size(batch)).*nan;
elseif length(outnlidx)==1,
   outnlidx=ones(size(batch)).*outnlidx;
end
outnlidx2=outnlidx;

rcsetstrings;
dbopen;
if isnumeric(runidx),
   srunid='(';
   for ii=1:length(rawids),
      srunid=[srunid,num2str(runidx(ii)),','];
   end
   srunid(end)=')';
   sql=['SELECT * from sRunData WHERE id in ',srunid];
   rundata=mysql(sql);
   batchcount=length(rundata);
   
   cellid=rundata.cellid;
   
else
   cellid=runidx;
   goodbatch=zeros(1,length(batch));
   batchcount=length(batch);
   for ii=1:batchcount,
      sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
           ' AND batch=',num2str(batch(ii))];
      trd=mysql(sql);
      if ~isempty(trd),
         rundata(ii)=trd;
         goodbatch(ii)=1;
         z{ii}=zload([rundata(ii).respath,rundata(ii).resfile,'.gz']);
         if isnan(outnlidx(ii)),
            outnlidx(ii)=z{ii}.params.nlidxsave;
            outnlidx2(ii)=z{ii}.params.nlidxsave;
         end
         if isfield(z{ii},'goodbatchrange'),
            if ~isempty(find(z{ii}.goodbatchrange==outnlidx(ii))),
               goodbatch(ii)=1;
               outnlidx2(ii)=find(z{ii}.goodbatchrange==outnlidx(ii));
            end
         else
            goodbatch(ii)=1;
         end
      end
   end
end

if sum(goodbatch)==0,
   disp('no entries found in db!');
   return
end

% convert to gr format?
%GREXT='gr';
GREXT='';

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
   
   titles{ii}=sprintf('%s Bat %d/nl %d: STRF Pred:%s',...
                      cellid,batch(ii),outnlidx(ii),cnfstr);
   
   if strcmp(z{ii}.params.stimfiltercmd,'movphasesep'),
      kernfmt=['fft',GREXT];
   elseif strcmp(z{ii}.params.stimfiltercmd,'movpower'),
      kernfmt=['pfft',GREXT];
   else
      kernfmt=z{ii}.params.kernfmt;
   end
   showkern(tsf,kernfmt,z{ii}.params.iconside,titles);
end
set(gcf,'PaperOrientation','landscape','PaperPosition',[0.25 0.25 10.5 8]);

%
% space-time decomposition
%
figure(2)
clf
if isfield(z{1}.params,'altcore') & ...
      strcmp(z{1}.params.altcore,'xccorefet'),
   

   for ii=find(goodbatch(:)'),
      subplot(length(goodbatch),2,ii*2-1);
      
      ss=z{ii}.strf(1).h;
      mss=max(abs(ss(:)));
      imagesc(ss,[-1 1].*mss);
      axis xy
      
      title(sprintf('batch %d',batch(ii)));
      
      subplot(length(goodbatch),2,ii*2);
      if ii==1,
         ss=z{ii}.strf(1).h;
         mss=max(abs(ss(:)));
         imagesc(ss,[-1 1].*mss);
         axis xy
      else
         % reverse time to match fet code
         CSR=fliplr(z{ii}.strf(1).h);
         nband=size(CSR,1);
         twindow=[-1 1] .* (size(CSR,2)-1);
         CSR=[CSR zeros(nband,twindow(2))];
         
         nstd_val = 0.5;
         [fstim, fstim_spike] = ...
             fft_AutoCrossCorr(z{1}.meanCS,CSR,[],twindow(2), nband, ...
                               nstd_val);
         
         nb = nband;
         nt = 2*twindow(2) +1;
         nJN = 1;
         stim_size = size(fstim);
         stim_spike_size = size(fstim_spike);
         stim_spike_JNsize = size(fstim_spike);
         
         tol=z{1}.Tol_val(z{1}.strf(1).parms.sfsfit);
         
         nf = (nt-1)/2 +1;
         fstim_spike_biased=zeros(size(fstim_spike));
         for iff=1:nf
            nc = 1;
            for ii=1:nb
               for jj=ii:nb
                  stim_mat(ii,jj) = fstim(nc,iff);
                  if ii ~= jj
                     stim_mat(jj,ii) = conj(fstim(nc,iff));
                  end
                  nc = nc +1;
               end
            end
            fstim_spike_biased(:,iff)=stim_mat*fstim_spike(:,iff);
            if iff>1,
               fstim_spike_biased(:,nt-iff+2)=...
                   stim_mat*fstim_spike(:,nt-iff+2);
            end
         end
         forward = ...
             cal_Strf(fstim,fstim_spike_biased, fstim_spike, stim_size, ...
                      stim_spike_size,stim_spike_JNsize, nb, nt, nJN, tol);
         ss=forward(:,twindow(2)+1:end);
         mss=max(abs(ss(:)));
         imagesc(ss,[-1 1].*mss);
         axis xy
         
         title(sprintf('residual biased tol=%.4f',tol));
         
      end
   end
   
   return
   
else
   
   
   
   if isfield(z{1}.strf(outnlidx2(1),1),'powunbiased') & ...
         sum(z{1}.strf(outnlidx2(1)).powunbiased)~=0,
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
      
      if isfield(z{ii}.strf(outnlidx2(ii),1),'powunbiased') & ...
            size(z{ii}.strf(outnlidx2(ii),1).powunbiased,1)==size(h,1) & ...
            sum(z{ii}.strf(outnlidx2(ii),1).powunbiased(:))~=0,
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
   
      if size(ratub,1)==size(h,1),
         hub=h .* repmat(ratub,1,tcount);
      else
         hub=h;
      end
      
      if isfield(z{ii}.strf(outnlidx2(ii),1),'hspace'),
         spvec=z{ii}.strf(outnlidx2(ii),1).hspace;
         spvecub=z{ii}.strf(outnlidx2(ii),1).hspacebiased;
         ttimesep=z{ii}.strf(outnlidx2(ii),1).tempresp;
         
         if  ii>1,
            disp('recomputing bias-normalized strf');
            %tt=z{ii}.mSA2;
            %tt=pinv(tt,0.000001);
            
            % re-compute bias normalized strf
            
            hspace=spvec;
            if isfield(z{1},'mSA2'),
               mSA2=z{1}.mSA2;
            else
               mSA2=z{1}.strf(outnlidx2(1),1).sSA2;
            end
            
            if size(mSA2,2)==size(hspace,1),
               hspace=mSA2*hspace;
               th=normalizereg(hspace,mSA2,[],z{1}.params.sfscount,...
                               z{1}.params.sfsstep);
               spvecub=th(:,:,z{1}.strf(outnlidx2(1),1).parms.sfsfit);
            else
               spvecub=hspace;
            end
            
            %hspace=mSA2*savestrf(bidx).hspace;
            
         end
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
         if tspaceall'*u(:,1) > 0;
            ttimesep=v(:,1);
            tspacesep=u(:,1);
         else
            ttimesep=-v(:,1);
            tspacesep=-u(:,1);
         end
         teigvals=diag(s);
         
         spvec=tspacesep;
         if size(ratub,1)==size(h,1),
            spvecub=spvec .* ratub;
         else
            spvecub=spvec;
         end
      end
      
      ssp=sort(spvec);
      sptopmean=mean(ssp(end-3:end));
      
      tH=ones(spacecount,4,length(goodbatch)).*nan;
      if 0,
         shrf=3; df=4;
         hidx=find(spvec>shrf*std(spvec));
         spvec(hidx)=shrf*std(spvec)+(spvec(hidx)-shrf*std(spvec))/df;
         hidx=find(spvec<-shrf*std(spvec));
         spvec(hidx)=-shrf*std(spvec)+(spvec(hidx)+shrf*std(spvec))/df;
         hidx=find(spvecub>shrf*std(spvecub));
         spvecub(hidx)=shrf*std(spvecub)+(spvecub(hidx)-shrf*std(spvecub))/df;
         hidx=find(spvecub<-shrf*std(spvecub));
         spvecub(hidx)=-shrf*std(spvecub)+(spvecub(hidx)+shrf*std(spvecub))/df;
      end
      tH(:,1,ii)=spvec./max(abs(spvec)).*sptopmean;
      tH(:,2,ii)=powunbiased./max(powunbiased).*sptopmean;
      tH(:,3,ii)=spvecub./max(abs(spvecub)).*sptopmean;
      tH(:,4,ii)=spvecub./max(abs(spvecub)).*sptopmean;
      
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
end

fullpage('portrait');
drawnow

%
% example prediction vs actual PSTH
%

% skip if not review
if z{1}.params.runclassid~=0,
   disp('not review. skipping stim for pred plot');
else
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
   sframes(sframes==0)=1;
   bigimagesc=ones(size(mov,1),round((eidx-sidx+1)*imagescale)).*255;
   for ii=1:size(mov,3),
      bigimagesc(1:size(mov,1),sframes(ii):sframes(ii)+size(mov,2)-1)=...
          fliplr(flipud(mov(:,:,ii)));
   end
end

figure(3)
clf

subplot(batchcount+1,1,1);
if z{1}.params.runclassid==0,
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
end

for ii=find(goodbatch(:)'),
   
   if exist('cnfshowidx','var'),
      cnfidx=cnfshowidx;
   else
      cnfidx=find(cat(1,z{ii}.params.predbatch{:})==z{1}.params.batch);
      if isempty(cnfidx),
         cnfidx=find(cat(1,z{ii}.params.predbatch{:})==z{ii}.params.batch);
      end
      if isempty(cnfidx),
         cnfidx=1;
      end
   end
   
   if ~isempty(cnfidx),
      
      tpred=z{ii}.predres(cnfidx).mod_psth{1}(:,1,outnlidx2(ii));
      
      if 0 & isfield(z{cnfidx},'times'),
         tact=feval(z{cnfidx}.params.resploadcmd,...
                    z{cnfidx}.params.respfiles{z{cnfidx}.times(3).fileidx},...
                    z{cnfidx}.params.resploadparms);
         tact=tact(z{cnfidx}.times(3).start: ...
                   z{cnfidx}.times(3).stop,1);
      else
         tact=z{ii}.predres(cnfidx).act_resp{1}(:,1);
      end
      
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
      
      % plot predictions
         
       if z{1}.params.runclassid==0,
         hs=subplot(batchcount+1,1,ii+1);
      
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
         
         if z{1}.params.runclassid==0,
            
            for jj=1:length(FirstFr),
               plot((FirstFr(jj)-sidx).*[14 14],[0 a(4)],'k:');
               %plot((LastFr(jj)-sidx).*[14 14],[0 a(4)],'r--');
            end
         end
         
         hold off
         a=axis;
         axis([tt(1) tt(end) -a(4)*0.05 a(4)]);
         if ii<max(find(goodbatch(:))),
            set(gca,'XTickLabel',[]);
         end
         
      else
         hs=subplot(batchcount,2,ii*2-1);
         
         tbinms=16;
         tpred=tpred./tbinms.*1000;
         tact=tact./tbinms.*1000;
         
         scatter(tpred,tact,'k.');
         a=axis;
         axis([a(3) a(4) a(3) a(4)]);
         axis square

         hs=subplot(batchcount,2,ii*2);
         scatter(tpred,tact,'k.');
         a=axis;
         axis([a(3) a(4) a(3) a(4)]);
         axis square
         
         mm=mean([tpred tact]);
         scf=tpred'*tact./(tpred'*tpred);
         
         hold on
         plot([0 a(4)],[a(3) a(4).*scf],'b--','LineWidth',2);
         hold off
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


