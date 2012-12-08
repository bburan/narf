% function r=cellres(cellid,batch)
%
% load results from resfile in sRunData and display fits, pred
% results
%
% r=0 if no entries found in db, =1 otherwise
%
function r=xcstsepres(cellid,batch)

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
fprintf('loading: %s\n',resfile);
zload(resfile);

global BATQUEUEID
BATQUEUEID=[];

batchcount=size(strf,1);

hshow=[];
for bidx=1:batchcount,
   if ~isfield(strf(bidx),'hspace'),
      hspace=sum(strf(bidx).h,2);
   else
      hspace=strf(bidx).hspace;
   end
   
   if isfield(strf(bidx),'hspacebiased'),
      hsb=strf(bidx).hspacebiased;
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
for bidx=1:batchcount,
   subplot(batchcount,colcount,(bidx-1)*colcount+1);
   title(sprintf('%s bat %d: predxc=%.3f',...
                 cellid,sparams.estbatch(goodbatchrange(bidx)),...
                 predxc(1,goodbatchrange(bidx))));
   
   subplot(batchcount,colcount,bidx*colcount);
   tr0=strf(bidx).tempresp0;
   tr1=strf(bidx).tempresp;
   tr0=tr0./std(tr0).*std(tr1);
   plot(tbins,tr1,'r-');
   hold on;
   plot(tbins,tr0,'--');
   hold off
   
   title('temp resp');
   
   fprintf('spbatch: %d. sfs=%d sig=%d\n',...
           sparams.estbatch(goodbatchrange(bidx)),...
           strf(bidx).parms.sfsfit,strf(bidx).parms.sigfit)
end
legend('new','old');
xlabel('time lag (ms)');

set(gcf,'PaperOrientation','Portrait','PaperPosition',[0.25 0.25 8 10.5]);


%
% plot range of decorrelated kernel possibilities
%


if ~exist('optspacebidx','var'),
   optspacebidx=length(strf);
end

figure(2);
clf

showkern(strf(optspacebidx).mH,strf(end).parms.kernfmt);

set(gcf,'PaperOrientation','Portrait','PaperPosition',[0.25 0.25 8 10.5]);

figure(3);
clf
imagesc(xc(:,:,params.nlidxsave));
colormap(hot);
hold on
plot(strf(optspacebidx).parms.sigfit,strf(optspacebidx).parms.sfsfit,'kx');
hold off
title('sig/sfs fit for optspacebidx kernel');

return


rcsetstrings;

%keyboard
if ~exist('postanal'),
   postanal=('cellxc');
end

if strcmp(postanal,'cellxc'),
   for ii=1:batchcount,
      z{ii}=zload([rundata(ii).respath,rundata(ii).resfile,'.gz']);
   end
   
   % display fit results for different sfs and sigsmooth values in
   % each attentional state with lin and NL outputs
   %
   % also compute and extract the corresponding kernel for each
   % attentional/linearity condition
   figure(1);
   clf
   for ii=1:batchcount,
      maxc=max(max(max(max(z{ii}.xc(:,:,:,:,1)))));
      minc=min(min(min(min(z{ii}.xc(:,:,:,:,1)))));
      if minc<-maxc & maxc>0,
         minc=-maxc;
      end
      nlidx=1;
      for latidx=1:z{ii}.latcount,
         for attidx=1:z{ii}.attcount,
            for nlidx=1:z{ii}.nlcount,
               if z{ii}.latcount>1,
                  subplot(z{ii}.attcount,z{ii}.latcount*batchcount,...
                          ((attidx-1)*z{ii}.latcount+latidx-1)*batchcount+ii);
               else
                  subplot(z{ii}.attcount,z{ii}.nlcount*batchcount,...
                          ((attidx-1)*z{ii}.nlcount+nlidx-1)*batchcount+ii);
               end
               txc=z{ii}.xc(:,:,latidx,nlidx,attidx);
               
               if abs(nansum(abs(txc(:)-nanmean(txc(:)))))>0,
                  imagesc(txc,[minc maxc]);
               end
               
               hold on
               plot(z{ii}.sigfit(latidx,nlidx,attidx),...
                    z{ii}.sfsfit(latidx,nlidx,attidx),'x');
               hold off
               title(sprintf('%s\nl=%d a=%d nl=%d',...
                             z{ii}.skernfile,latidx,attidx,nlidx),...
                     'FontSize',10);
            end
         end
      end
      colorbar
   end
   colormap(hot);
   
   figure(2);
   if z{1}.respfmtcode==1,
      for ii=1:(size(z{1}.hf,2)*size(z{1}.hf,3)),
         z{1}.hf(:,ii)=z{1}.hf(:,ii)./mean(abs(z{1}.hf(:,ii)));
      end
   end
   tsf=z{1}.mH(:,(1:size(z{1}.hf,2))-z{1}.maxlag(1),1);
   titles{1}=sprintf('%s: Raw STA (%s)',cellid,z{1}.skernfile);
   if batchcount==1 & z{1}.attcount==1,
      tsf=cat(3,tsf,z{1}.hf);
      for ii=1:z{1}.nlcount,
         if z{1}.sfsfit(ii)==0,
            z{1}.sfsfit(ii)=1;
            z{1}.sigfit(ii)=1;
         end
         
         titles{ii+1}=sprintf('%s: STRF nl=%d xc=%0.3f pred=%.3f',...
                              cellid,ii,z{1}.xc(z{1}.sfsfit(ii),...
                                                z{1}.sigfit(ii),1,ii),...
                              z{1}.predxc(1,1,1,ii*3-1));
      end
   elseif batchcount==1,
      tsf=cat(3,tsf,squeeze(z{1}.hf(:,:,1,:)));
      for ii=1:z{1}.attcount,
         titles{ii+1}=sprintf('%s: STRF attidx=%d xc=%0.3f',...
                              cellid,ii,...
                              z{1}.xc(z{1}.sfsfit(1,1,ii),...
                                      z{1}.sigfit(1,1,ii),1,ii));
      end
   else
      for ii=1:batchcount,
         tsf=cat(3,tsf,z{ii}.hf(:,:,1));
         titles{ii+1}=sprintf('%s: STRF (%s) Pred=%0.3f',...
                              cellid,z{ii}.skernfile,z{ii}.predxc(1));
      end
   end

   kernfmt=z{1}.kernfmt;
   if strcmp(kernfmt,'pfftgr'),
      kernfmt='pfft';
      z{1}.iconside=[sqrt(z{1}.spacecount*2) sqrt(z{1}.spacecount*2)];
   elseif strcmp(kernfmt,'pixel') & length(z{ii}.iconside)<2,
      z{ii}.iconside=[sqrt(z{1}.iconside) sqrt(z{1}.iconside)];
   end
   %keyboard
   showkern(tsf,kernfmt,z{1}.iconside,titles);
   %showkern(squeeze(z{ii}.mSR),z{ii}.kernfmt,z.iconside)
   %showkern(squeeze(mH(:,:,end,:,:)),kernfmt,iconside)
   
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
   if ~isfield(z{1},'DOTHRESH'),
      z{1}.DOTHRESH=0;
   end
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
   showkern(tsf,kernfmt,z{1}.iconside,titles);
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
      %keyboard
      for ii=1:batchcount,
         predcount=size(z{ii}.predxc,1);
         fprintf('File: %s\n',z{ii}.skernfile);
         if size(z{ii}.predxc,4)>=9 & size(z{ii}.predxc,4)>z{ii}.nlcount,
            nlrange=[2 5 8];
         else
            nlrange=1:z{ii}.nlcount;
         end
         for nlidx=1:length(nlrange),
            fprintf('nl=%d: (sfs,sig)=(%2d,%2d)',...
                    nlidx,z{ii}.sfsfit(nlidx),z{ii}.sigfit(nlidx));
            for predidx=1:predcount,
               fprintf('%7.3f',z{ii}.predxc(predidx,:,:,nlrange(nlidx)));
               if z{ii}.predp(predidx,:,:,nlrange(nlidx))<0.01,
                  fprintf('*');
               else
                  fprintf(' ');
               end
            end
            fprintf('\n');
         end
      end
   end
elseif strcmp(postanal,'kerncomp'),
   
   if length(rundata)>1,
      disp('cellres-kerncomp: only displaying first batch');
   end
   
   fprintf('Loading %s...\n',[rundata(1).respath,rundata(1).kernfile,'.gz']);
   z=zload([rundata(1).respath,rundata(1).kernfile,'.gz']);
   
   attcount=z.attcount;
   attuse=z.attuse;
   spacelim=z.spacelim;
   attdocount=z.attdocount;
   noisesets=z.noisesets;
   respcount=z.respcount;
   attpairs=z.attpairs;
   bincount=z.bincount;
   T2=zeros(attdocount,noisesets,spacelim,respcount);
   attmod=zeros(attdocount,spacelim,respcount);
   attdiff=zeros(attdocount,spacelim,respcount);
   stimmod=zeros(attuse,spacelim,respcount);
   p=ones(attdocount,spacelim,respcount);
   r2=reshape(z.rbinned(:,:,2:end,:),bincount,respcount,attuse-1,...
              noisesets,spacelim);
   r2std=reshape(z.rbinnedstd(:,:,2:end,:),bincount,respcount,attuse-1,...
              noisesets,spacelim);
   for respidx=1:respcount,
      for pidx=1:attdocount,
         a1=attpairs(pidx,1);
         a2=attpairs(pidx,2);
         for eigidx=1:spacelim,
            x1=squeeze(permute(r2(:,respidx,a1,:,eigidx),[1 5 2 3 4]));
            x1=x1-repmat(mean(x1,1),[bincount,1,1]);
            x1=reshape(x1,bincount,noisesets);
            x2=permute(r2(:,respidx,a2,:,eigidx),[1 5 2 3 4]);
            x2=x2-repmat(mean(x2,1),[bincount,1,1]);
            x2=reshape(x2,bincount,noisesets);
            x1std=permute(r2std(:,respidx,a1,:,eigidx),[1 5 2 3 4]);
            x1std=reshape(x1std,bincount,noisesets);
            x2std=permute(r2std(:,respidx,a2,:,eigidx),[1 5 2 3 4]);
            x2std=reshape(x2std,bincount,noisesets);
            xstd=(x1std+x2std)./2;
            xstd(find(xstd==0))=1;
            T2(pidx,:,eigidx,respidx)=sum(abs(x1-x2)./xstd,1)./ ...
                (sqrt(sum(x1.^2./xstd,1)).*sqrt(sum(x2.^2./xstd,1)));
            
            
            x1=z.rbinned(:,respidx,a1+1,eigidx);
            x1std=z.rbinnedstd(:,respidx,a1+1,eigidx);
            %x1std=z.rbinnedstd(:,respidx,a1+1,eigidx).* ...
            %      sqrt(z.rbincount(:,respidx,a1+1));
            x2=z.rbinned(:,respidx,a2+1,eigidx);
            x2std=z.rbinnedstd(:,respidx,a2+1,eigidx);
            %x2std=z.rbinnedstd(:,respidx,a2+1,eigidx).* ...
            %    sqrt(z.rbincount(:,respidx,a2+1));
            d=x2-x1;
            dstd=(x1std+x2std)./2;
            attmod(pidx,eigidx,respidx)=sqrt(mean(((d-mean(d))./dstd).^2));
            attdiff(pidx,eigidx,respidx)=sqrt(mean((d./dstd).^2));
         end
      end
      for eigidx=1:spacelim,
         for pidx=1:attdocount,
            tt2=T2(:,2:end,eigidx,respidx);
            ttest=T2(pidx,1,eigidx,respidx);
            tt2=sort([tt2(:);ttest]);
            p(pidx,eigidx,respidx)=1-(min(find(ttest<=tt2))-1)./length(tt2);
         end
      end
   end
   
   for respidx=1:respcount,
      for attidx=1:attuse,
         for eigidx=1:spacelim,
            x1=z.rbinned(:,respidx,attidx,eigidx);
            x1std=z.rbinnedstd(:,respidx,attidx,eigidx);
            %x1std=z.rbinnedstd(:,respidx,attidx,eigidx).* ...
            %      sqrt(z.rbincount(:,respidx,attidx));
            
            %stimmod=zeros(attuse,spacelim,respcount);
            stimmod(attidx,eigidx,respidx)=...
                sqrt(mean(((x1-mean(x1))./x1std).^2));
         end
      end
   end
   
   figure(1);
   clf
   respidx=7;
   smark={'k','g','c','b','r'};
   %rm=reshape(z.brmean(:,:,2:end),respcount,attuse-1,noisesets);
   %rmall=squeeze(mean(rm(respidx,:,1)));
   rmall=0;
   rmin=min(min(min(min(z.rbinned(:,respidx,2:attcount,:)))))-rmall;
   rmax=max(max(max(max(z.rbinned(:,respidx,2:attcount,:)))))-rmall;
   %spaceuse=spacelim;
   spaceuse=15;
   rowcount=floor(spaceuse/5);
   for eigidx=1:spaceuse,
      subplot(rowcount,ceil(spaceuse/rowcount),eigidx);
      
      if 0,
         imagesc(squeeze(z.rbinned(:,respidx,1:attuse,eigidx))',[rmin,rmax]);
         hold on
         for attidx=1:attuse,
            plot(z.mpbinned(attidx,eigidx),attidx,'kx');
            plot(z.mpbinned(attidx,eigidx),attidx,'wo');
         end
         hold off
         
         axis off;
      else
         hold on
         labels={'A','1','2','3','4'};
         ll=zeros(attuse,1);
         attshowset=[1 2 3 4 5];
         for attidx=attshowset,
            r=z.rbinned(:,respidx,attidx,eigidx);
            %r=z.rbinned(:,respidx,attidx,eigidx)-...
            %  mean(z.rbinned(:,respidx,attidx,eigidx));
            a=errorbar(z.sbin(:,eigidx),r,...
                       squeeze(z.rbinnedstd(:,respidx,attidx,eigidx)).*2,...
                       [smark{attidx},'-']);
            set(a,'LineWidth',2);
            ll(attidx)=a(1);
            plot(z.sbin(z.mpbinned(attidx,eigidx),eigidx),...
                 r(z.mpbinned(attidx,eigidx)),'kx');
         end
         hold off
         axis([z.sbin(1,eigidx),z.sbin(end,eigidx),rmin,rmax]);
         axis([z.sbin(1,eigidx),z.sbin(end,eigidx),rmin,rmax]);
         if eigidx==spaceuse,
            legend(ll(find(ll)),labels{find(ll)});
         end
      end
      title(sprintf('%s: pc%d',cellid,eigidx));
   end
   colormap(hot);
   %set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8]);
   set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 2 10.5 4]);
   
   spaceuse=10;
   spacecount=z.spacecount;
   eigmtx=z.ua(:,1:spaceuse);
   eigtargets=eigmtx'*(z.mpatches-repmat(z.bsmean',[1 5]));
   
   eigmtx=repmat(eigmtx,[1 1 attuse]);
   targmtx=reshape(z.mpatches,spacecount,1,attuse);
   for ii=1:attuse,
      for jj=1:spaceuse,
         eigmtx(:,jj,ii)=eigmtx(:,jj,ii).*eigtargets(jj,ii) ...
             ./ mean(abs(eigtargets(jj,2:end)));
      end
      %targmtx(:,:,ii)=targmtx(:,:,ii)-targmtx(:,:,1);
   end
   
   for ii=1:attuse,
      targmtx(:,:,ii)=targmtx(:,:,ii) ./ ...
          max(abs(targmtx(:,:,ii))) .* max(abs(eigmtx(:)));
   end
   eigmtx=cat(2,targmtx,eigmtx);
   %for ii=1:attuse,
   %   eigmtx(:,:,ii)=eigmtx(:,:,ii) ./ max(max(abs(eigmtx(:,:,ii))));
   %end
   
   figure(2);
   kernfmt=z.kernfmt;
   if strcmp(kernfmt,'pfftgr'),
      kernfmt='pfft';
      z.iconside=[sqrt(spacecount*2) sqrt(spacecount*2)];
   end
   showkern(eigmtx,kernfmt,z.iconside,{},0);
   colorbar
   for eigidx=0:spaceuse,
      subplot(attuse,spaceuse+1,eigidx+1);
      if eigidx==0,
         title(sprintf('cell %s',cellid));
      else
         title(sprintf('pc%d',eigidx));
      end
      xlabel('');
   end
   for attidx=1:attuse,
      subplot(attuse,spaceuse+1,attidx*(spaceuse+1)-spaceuse);
      ylabel(sprintf('attidx=%d',attidx));
   end
   %spacecount=z.spacecount;
   %eigmtx=z.ua(:,1:(ceil(spaceuse/rowcount)*rowcount));
   %eigmtx=reshape(eigmtx,spacecount,ceil(spaceuse/rowcount),rowcount);
   %figure(2);
   %kernfmt=z.kernfmt;
   %if strcmp(kernfmt,'pfftgr'),
   %   kernfmt='pfft';
   %end
   %showkern(eigmtx,kernfmt,z.iconside);
   %for eigidx=1:spaceuse,
   %   subplot(rowcount,ceil(spaceuse/rowcount),eigidx);
   %   title(sprintf('%s: pc%d',cellid,eigidx));
   %   xlabel('');
   %end
   
   figure(3);
   clf
   for ii=1:length(z.targlist),
      subplot(length(z.targlist),1,ii);
      imagesc(z.bigpatches(:,:,z.targlist(ii)));
      title(sprintf('target %d (id %d)',ii,z.targlist(ii)));
      axis off
      axis image
   end
   colormap(gray);
   
   figure(4);
   clf
   %semilogy(p(:,:,respidx)');
   %legend('1-2','1-3','1-4','2-3','2-4','3-4');
   
   subplot(1,3,1);
   imagesc(stimmod(:,:,respidx));
   colorbar
   title('stimmod');
   subplot(1,3,2);
   imagesc(attdiff(:,:,respidx));
   colorbar
   title('attdiff: 1-2,1-3,1-4,2-3,2-4,3-4');
   subplot(1,3,3);
   imagesc(attmod(:,:,respidx));
   colorbar
   title('attmod');
   
   colormap(hot);
   
   %keyboard
   
   return
   
   
   PCCOUNT=15;
   latbase=4;
   latidx=5;
   
   if 0
   keyboard
   
   noisecount=z.attcount-z.attuse;
   bootcount=5000;
   T2=zeros(z.attuse-1,z.attuse-1);
   ppair=zeros(z.attuse-1,z.attuse-1);
   nT2=zeros(bootcount,1);
   
   % compute variance matrix of noise kernels.
   X=squeeze(z.eigH(1:PCCOUNT,latidx,z.attuse+1:end));
   mX=mean(X,2);
   X0=X-repmat(mX,[1 noisecount]);
   V=X0*X0' ./noisecount;
   
   for bidx=1:bootcount,
      att1=ceil(rand*noisecount);
      att2=ceil(rand*noisecount);
      E=z.eigH(1:PCCOUNT,latidx,att1+z.attuse) - ...
        z.eigH(1:PCCOUNT,latidx,att2+z.attuse);
      nT2(bidx)=E' * V^-1 * E ./ PCCOUNT;
   end
   
   for att1=1:z.attuse-1,
      for att2=att1+1:z.attuse-1,
         
         % mini hoteling's T kernel pair analysis
         E=z.eigH(1:PCCOUNT,latidx,att1+1) - ...
           z.eigH(1:PCCOUNT,latidx,att2+1);
         T2(att1,att2)=E' * V^-1 * E ./ PCCOUNT;
         
         sT2=sort([nT2; T2(att1,att2)]);
         ppair(att1,att2)=1-max(find(T2(att1,att2)>sT2))./(bootcount+1);
      end
   end
   end
   
   pgood=zeros(z.attuse,1);
   
   figure(1);
   clf
   for attidx=1:z.attuse,
      
      % plot kernel in PC domain
      subplot(z.attuse,3,attidx*3-2);
      
      errorbar(1:z.spacelim,z.eigH(:,latidx,attidx),z.seigH(:,latidx)*2,'b-');
      %plot(1:z.spacelim,z.mH(:,attidx)./z.eH(:,attidx));
      
      hold on
      plot([1 z.spacelim],[0 0],'k:');
      %plot(1:z.spacelim,z.eigH(:,latidx,attidx),'b-');
      errorbar(1:z.spacelim,z.eigHout(:,latidx,attidx),...
               z.seigHout(:,latidx),'k--');
      hold off
      %T2=zeros(attcount,spacelim,resptouse);
      %p=zeros(attuse,spacelim,resptouse);
      title(sprintf('%s/%d T%d=%.2f p<%.3f',...
                    z.cellid,z.rundata.batch,PCCOUNT,...
                    z.T2(attidx,PCCOUNT,latidx-latbase),...
                    z.p(attidx,PCCOUNT,latidx-latbase)));
      pgood(attidx)=z.p(attidx,PCCOUNT,latidx-latbase);
      
      ylabel(sprintf('att %d',attidx));
      if attidx==z.attuse;
         xlabel('eigenvector');
      end
      a=axis;
      axis([0 z.spacelim+1 a(3) a(4)]);
      
      % plot target in PC domain
      subplot(z.attuse,3,attidx*3-1);
      %errorbar(1:z.spacelim,z.meigpatches,z.seigpatches*2,'k--');
      plot(1:z.spacelim,z.attcc(:,attidx,1),'b-');
      hold on
      plot([1 z.spacelim],[0 0],'k:');
      plot(1:z.spacelim,z.attcc(:,attidx,2),'r--');
      plot(1:z.spacelim,z.attcc(:,attidx,3),'k:');
      %plot(1:z.spacelim,z.teigpatches(:,attidx),'r')
      hold off
      title(sprintf('targcc:'));
      %title(sprintf('targxc: %.3f/%.3f p<%.3f',...
      %              z.attsimxc(attidx,attidx,PCCOUNT),...
      %              mean(z.attsimxc(:,attidx,PCCOUNT)),...
      %              z.ptarg(attidx,attidx,PCCOUNT)));
      if attidx==z.attuse;
         xlabel('eigenvector');
      end
      a=axis;
      axis([0 z.spacelim+1 a(3) a(4)]);
      
      % plot kernel & target in PC domain, SE units
      subplot(z.attuse,3,attidx*3);
      tdev=(z.teigpatches(:,attidx)-z.meigpatches)./z.seigpatches;
      %kdev=(z.mH(:,attidx)-z.meigH) ./z.seigH;
      kdev=(z.eigH(:,latidx,attidx)-z.eigHout(:,latidx,attidx)) ./ ...
           z.seigHdiff(:,latidx);
      if attidx>1,
         tkcc=xcorr(z.eigH(1:PCCOUNT,latidx,attidx),...
                    z.teigpatches(1:PCCOUNT,attidx),0,'coeff');
         tkccout=xcorr(z.eigHout(1:PCCOUNT,latidx,attidx),...
                       z.teigpatches(1:PCCOUNT,attidx),0,'coeff');
      else
         tkcc=0;
         tkccout=0;
      end
      tkccn=zeros(z.attcount,2);
      
      if 0
      for nidx=1:z.attcount,
         tkccn(nidx,1)=xcorr(z.eigH(1:PCCOUNT,latidx,attidx),...
                             z.teigpatches(1:PCCOUNT,attidx),0,'coeff');
         if attidx>1,
            tkccn(nidx,2)=xcorr(z.eigHout(1:PCCOUNT,latidx,attidx),...
                             z.teigpatches(1:PCCOUNT,attidx),0,'coeff');
         end
      end
      ttcc=sort([tkcc; tkccn(z.attuse+1:end)]);
      ttpp=1-min(find(tkcc<=ttcc))./length(ttcc);
      end
      
      plot(1:z.spacelim,tdev,'r');
      hold on
      plot(1:z.spacelim,kdev,'b')
      plot([1 z.spacelim],[0 0],'k:');
      hold off
      axis([0 z.spacelim+1 -5 5]);
      title(sprintf('cc=%.3f ccout=%.3f pc=%d',tkcc,tkccout,PCCOUNT));
      ylabel('stderr');
      if attidx==z.attuse;
         xlabel('eigenvector');
      end
   end
   legend('target','kernel');
   set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8]);

   figure(2);
   clf
   showstep=2;
   showdim=10*showstep;
   if showdim>z.spacelim,
      showstep=1; showdim=10*showstep;
   end
   rowcount=(z.attuse-1)*2;
   colcount=showdim/showstep+1;
   ut=zeros(z.spacecount,showdim,rowcount);
   for attidx=1:z.attuse-1,
      % scale each stimulus eigenvector by kernel coefficient
      ut(:,:,attidx*2-1) = z.uu(:,1:showdim) * ...
          diag(z.eigH(1:showdim,latidx,attidx+1));
      titles{attidx*2-1}=sprintf('Cell %s att %d in (p<%.3f)',...
                                 z.cellid,attidx,pgood(attidx+1));
      ut(:,:,attidx*2) = z.uu(:,1:showdim,1) * ...
          diag(z.eigHout(1:showdim,latidx,attidx+1));
      titles{attidx*2}=sprintf('Cell %s att %d out',z.cellid,attidx);
   end
   
   % normalize to match scale of targets
   ut=ut./max(max(z.eigH(1:showdim,latidx,2:z.attuse)));
   uts=cumsum(ut,2);
   uts=uts(:,showstep:showstep:showdim,:);
   %uts(1,1,:)=max(max(uts));
   
   %keyboard
   %upatch=z.mpatches(:,2:z.attuse);
   %fpatch=z.fpatches;
   
   upatch=z.uu(:,1:showdim)*z.teigpatches(1:showdim,2:z.attuse);
   fpatch=z.uu(:,1:showdim)*z.eigpatches(1:showdim,:);
   %fpatch=z.uu(:,1:showdim)*z.beigpatches(1:showdim,:);
   
   ipatch=zeros(size(upatch));
   mpatch=mean(fpatch,2);
   spatch=std(fpatch,1,2);
   spatch(find(spatch==0))=1;
   for attidx=1:(z.attuse-1),
      ipatch(:,attidx)=(upatch(:,attidx)-mpatch)./spatch;
      ipatch(:,attidx)=ipatch(:,attidx) ...
          ./max(abs(ipatch(:,attidx))) .* max(uts(:));
      upatch(:,attidx)=upatch(:,attidx) ./ ...
          max(abs(upatch(:,attidx))) .* max(uts(:));
   end
   upatch=reshape(cat(1,upatch,ipatch),z.spacecount,1,rowcount);
   
   %upatch=reshape(repmat(upatch,[2,1]),z.spacecount,1,rowcount);
   
   if strcmp(z.kernfmt,'gr'),
      iconside=[8 8];
   elseif strcmp(z.kernfmt,'orsf'),
      iconside=[16 1];
   else
      iconside=strsep(z.cellfiledata(1).stimiconside,',');
      iconside=cat(2,iconside{:});
   end
   
   %showkern(cat(2,upatch,ut,uts),z.kernfmt,iconside,titles);
   showkern(cat(2,upatch,uts),z.kernfmt,iconside,titles);
   subplot(rowcount,colcount,(rowcount-1)*colcount+1);
   xlabel('target');
   for showidx=1:showdim/showstep,
      subplot(rowcount,colcount,(rowcount-1)*colcount+showidx+1);
      xlabel(sprintf('eig %d',showidx*showstep));
   end
   for attidx=1:(z.attuse-1),
      subplot(rowcount,colcount,(attidx*2-2)*colcount+1);
      title('matched');
      subplot(rowcount,colcount,(attidx*2-1)*colcount+1);
      title('ideal');
   end
   
         
   set(gcf,'PaperOrientation','Landscape','PaperPosition',[0.25 0.25 10.5 8]);
   
   %subplot(2,1,1);
   %hl=plot(z.mmstim(:,1:z.attuse));
   %title([z.cellid,' mean stim']);
   %
   %subplot(2,1,2);
   %hl=plot(z.mmresp(:,1:z.attuse));
   %set(hl,'LineWidth',2);
   %title([z.cellid,' mean resp']);
   %
   %legend('A','1','2','3','4');
   
   % temporarily disable display of targets in image domain
   if 0,
   figure(3);
   clf
   for ii=1:length(z.targlist),
      subplot(length(z.targlist),1,ii);
      imagesc(z.bigpatches(:,:,z.targlist(ii)));
      title(sprintf('target %d (id %d)',ii,z.targlist(ii)));
      axis off
      axis image
   end
   colormap(gray);
   end
   
   if 0,
   figure(2);
   clf
   for xx=1:z.spacelim,
      subplot(ceil(z.spacelim/5),5,xx);
      imagesc(squeeze(z.ucorr(xx,:,:)),[0 1]);
      axis image
      axis off
      title(sprintf('%d',xx));
      if xx==1,
         ylabel(z.cellid);
      end
   end
   colormap(gray);
   
   figure(3);
   clf
   for xx=1:z.spacelim,
      subplot(ceil(z.spacelim/5),5,xx);
      imagesc(squeeze(z.ediff(xx,:,:)),[-3 3]);
      axis image
      axis off
      title(sprintf('%d',xx));
      if xx==1,
         ylabel(z.cellid);
      end
   end
   colormap(hot);
   end
   
end


drawnow;

if nargout>0,
   r=z{1};
end

