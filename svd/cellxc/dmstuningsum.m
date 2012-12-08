
if ~exist('recalc','var'),
   recalc=1;
end

if ismember(batchid,[95 96]),
   FVVS=1;
else
   FVVS=0;
end

% load pre-saved res data (from dmsattsum.m)
RESPATH='/auto/k5/david/tmp/kvaparms/';
resfile=sprintf('%sdmsattsum.batch%d.mat',RESPATH,batchid);
fprintf('loading %s\n',resfile);
load(resfile);

keepidx=[];
for ii=1:length(res),
   if ~isempty(res(ii).cellid),
      keepidx=[keepidx; ii];
   end
end
res=res(keepidx);
cellcount=length(res);
if recalc,
   disp('calculating tuning curves...');
   %acount=1;
   acount=size(res(1).Hset,2);
   hcount=size(res(1).Hset,3);
   
   % MTS data:
   % aidx=1: feature attention
   % aidx=2: spatial attention
   
   if FVVS,
      patchcount=acount;
   else
      patchcount=hcount;
   end
   obincount=16;
   sfbincount=8;
   
   pamp=zeros(acount,hcount,cellcount);
   namp=zeros(acount,hcount,cellcount);
   seprat=zeros(acount,hcount,cellcount);
   orpeak=zeros(acount,hcount,cellcount);
   orbw=zeros(acount,hcount,cellcount);
   porpeak=zeros(acount,hcount,cellcount);
   porbw=zeros(acount,hcount,cellcount);
   norpeak=zeros(acount,hcount,cellcount);
   norbw=zeros(acount,hcount,cellcount);
   sfpeak=zeros(acount,hcount,cellcount);
   sfbw=zeros(acount,hcount,cellcount);
   psfpeak=zeros(acount,hcount,cellcount);
   psfbw=zeros(acount,hcount,cellcount);
   nsfpeak=zeros(acount,hcount,cellcount);
   nsfbw=zeros(acount,hcount,cellcount);
   targorpeak=zeros(patchcount,cellcount);
   targsfpeak=zeros(patchcount,cellcount);
   
   sfmatchxc=zeros(acount,hcount,cellcount);
   sfdiffxc=zeros(acount,hcount,cellcount);
   ormatchxc=zeros(acount,hcount,cellcount);
   ordiffxc=zeros(acount,hcount,cellcount);
   hmatchxc=zeros(acount,hcount,cellcount);
   hdiffxc=zeros(acount,hcount,cellcount);
   
   tpatch=[];
   for ii=1:length(res),
      tpatch=cat(2,tpatch,res(ii).mpatches(:,1));
   end
   tpatch=mean(tpatch,2);
   targsf=zeros(sfbincount,patchcount);
   targor=zeros(obincount,patchcount);
   targsfgr=zeros(obincount,sfbincount,patchcount);
   bincount=size(res(ii).Hset(:,aidx,hidx),1);
   targfft=zeros(bincount,patchcount);
   
   for ii=1:length(res),
      
      fprintf('%s - computing tuning curves\n',res(ii).cellid);
      %smmat=[.1 .8 .1];
      smmat=[0 1 0];
      
      for hidx=1:min([patchcount,size(res(ii).mpatches,2)]),
         tp=res(ii).mpatches(:,hidx);
         tp=tp-res(ii).mS(:,1);
         tsf0=pfft2sf(tp);
         [tsfgr,obins,sfbins]=sf2gr(tsf0,obincount,sfbincount);
         
         targfft(:,hidx)=tp;
         targsfgr(:,:,hidx)=tsfgr;
         
         tsfgr=rconv2(tsfgr,smmat);
         tsfgr=cconv2(tsfgr,smmat');
         
         % use mean, svd seems to give worse answer
         targsf(:,hidx)=mean(tsfgr,1)';
         targor(:,hidx)=mean(tsfgr,2);
         
         mm=min(find(tsfgr==max(tsfgr(:))));
         [mor,msf]=ind2sub([obincount,sfbincount],mm);
         
         targorpeak(hidx,ii)=obins(mor);
         
         sfnorm=sum(tsfgr.*(tsfgr>0),1)';
         sfsum=sfnorm.*log2(sfbins(:));
         if sum(sfsum)==0,
            sfsum(:)=0.1;
         end
         sfsum(sfsum==0)=min(sfsum(sfsum>0))./1000;
         sfsum=sfsum./(sum(sfsum) + (sum(sfsum)==0));
         sfsum=cumsum(sfsum);
         sfpoints=2.^interp1(sfsum,log2(sfbins),[.125 .5 .875]);
         
         targsfpeak(hidx,ii)=sfpoints(2);
         
      end
      
      for hidx=1:min([hcount,size(res(ii).Hset,3)]),
         for aidx=1:min([acount,size(res(ii).Hset,2)]),
            tH=res(ii).Hset(:,aidx,hidx);
            tsfIR=pfft2sf(res(ii).Hset(:,aidx,hidx),'pfft');
            [tsfgr,obins,sfbins]=sf2gr(tsfIR,obincount,sfbincount);
            
            tsfgr_sm=rconv2(tsfgr,smmat);
            tsfgr_sm=cconv2(tsfgr_sm,smmat');
            
            thissf=mean(tsfgr_sm,1)';
            thisor=mean(tsfgr_sm,2);
            
            mm=min(find(tsfgr_sm==max(tsfgr_sm(:))));
            [mor,msf]=ind2sub([obincount,sfbincount],mm);
            
            orpeak(aidx,hidx,ii)=obins(mor);
            %sfpeak(aidx,hidx,ii)=sfbins(msf);
            
            sfnorm=sum(tsfgr_sm.*(tsfgr_sm>0),1)';
            sfsum=sfnorm.*log2(sfbins(:));
            if sum(sfsum)==0,
               sfsum(:)=0.1;
            end
            sfsum(sfsum==0)=min(sfsum(sfsum>0))./1000;
            sfsum=sfsum./(sum(sfsum) + (sum(sfsum)==0));
            sfsum=cumsum(sfsum);
            sfpoints=2.^interp1(sfsum,log2(sfbins),[.125 .5 .875]);
            
            sfpeak(aidx,hidx,ii)=sfpoints(2);
            
            if FVVS,
               
               if aidx==1,
                  aoutrange=1;
               else
                  aoutrange=[2:(aidx-1) (aidx+1):patchcount];
               end
               ormatchxc(aidx,hidx,ii)=xcov(thisor,targor(:,aidx),0,'coeff');
               sfmatchxc(aidx,hidx,ii)=xcov(thissf,targsf(:,aidx),0,'coeff');
               aa=tsfgr_sm;
               bb=targsfgr(:,:,aidx);
               hmatchxc(aidx,hidx,ii)=xcov(aa(1:(end)),bb(1:(end)),0,'coeff');
               
               for aout=aoutrange,
                  ordiffxc(aidx,hidx,ii)=ordiffxc(aidx,hidx,ii)+...
                      xcov(thisor,targor(:,aout),0,'coeff')./length(aoutrange);
                  sfdiffxc(aidx,hidx,ii)=sfdiffxc(aidx,hidx,ii)+...
                      xcov(thissf,targsf(:,aout),0,'coeff')./length(aoutrange);
                  bb=targsfgr(:,:,aout);
                  hdiffxc(aidx,hidx,ii)=xcov(aa(1:(end)),bb(1:(end)),0,'coeff');
               end
               
            else
               
               if hidx==1,
                  hout=1;
               else
                  hout=5-hidx;
               end
               % xcov or xcorr?
               ormatchxc(aidx,hidx,ii)=xcov(thisor,targor(:,hidx),0,'coeff');
               ordiffxc(aidx,hidx,ii)=xcov(thisor,targor(:,hout),0,'coeff');
               sfmatchxc(aidx,hidx,ii)=xcov(thissf,targsf(:,hidx),0,'coeff');
               sfdiffxc(aidx,hidx,ii)=xcov(thissf,targsf(:,hout),0,'coeff');
               
               aa=tsfgr_sm;
               bb=targsfgr(:,:,hidx);
               hmatchxc(aidx,hidx,ii)=xcov(aa(1:(end-16)),bb(1:(end-16)),0,'coeff');
               bb=targsfgr(:,:,hout);
               hdiffxc(aidx,hidx,ii)=xcov(aa(1:(end-16)),bb(1:(end-16)),0,'coeff');
            end
            
         end
      end
   end
   
   clear recalc
   save(resfile);
end

PTHRESH=0.05;

if batchid==96,

   predxc=cat(3,res.predxc);
   pxc=cat(2,res.pxc)';
   origlocalsig=pxc(:,2:3)<PTHRESH;
   
   ftsig=squeeze(pxc(:,[2 3 5])<PTHRESH);
   
   ftshiftidx=find(ftsig(:,3));
   
   
   figure(5);
   clf
   subplot(3,2,1);
   aa=mean(sfmatchxc(2:5,1,ftshiftidx)-sfdiffxc(2:5,1,ftshiftidx),1);
   hist(aa,linspace(-0.4,0.4,8));
   title(sprintf('sf xc shift, mean=%.3f',mean(aa)));
   
   subplot(3,2,3);
   aa=mean(ormatchxc(2:5,1,ftshiftidx)-ordiffxc(2:5,1,ftshiftidx),1);
   hist(aa,linspace(-0.8,0.8,8));
   title(sprintf('or xc shift, mean=%.3f',mean(aa)))
   
   subplot(3,2,5);
   aa=mean(hmatchxc(2:5,1,ftshiftidx)-hdiffxc(2:5,1,ftshiftidx),1);
   hist(aa,linspace(-0.6,0.6,12))
   title(sprintf('h xc shift, mean=%.3f',mean(aa)))
   for ii=1:length(ftshiftidx),
      jj=ftshiftidx(ii);
      fprintf('%s: xc0=%.3f xcdcg=%.3f xctun=%.3f tsi=%.3f\n',...
              res(jj).cellid,fullxc([1 3 4],2,jj),aa(ii));
   end
   
else
   
   predxc=cat(3,res.predxc);
   pxc=cat(2,res.pxc)';
   origlocalsig=pxc(:,2:3)<PTHRESH;
   
   fullxc=cat(3,res.fullxc);
   fullp=cat(3,res.fullp);
   
   ftsig=squeeze(fullp(2:4,2,:)<PTHRESH)';
   
   
   ftshiftidx=find(ftsig(:,3));
   %spshiftidx=find(spsig(:,3));
   
   
   figure(5);
   clf
   subplot(3,2,1);
   aa=sum(sfmatchxc(1,2:3,ftshiftidx)-...
           sfdiffxc(1,2:3,ftshiftidx),2);
      %./(2-sum(sfdiffxc(1,2:3,ftshiftidx),2));
   hist(aa,linspace(-0.6,0.6,12));
   title(sprintf('sf xc shift, mean=%.3f',mean(aa)));
   subplot(3,2,2);
   aa=sum(sfmatchxc(2,2:3,ftshiftidx)-...
           sfdiffxc(2,2:3,ftshiftidx),2);
      %./(2-sum(sfdiffxc(2,2:3,ftshiftidx),2));
   hist(aa,linspace(-0.6,0.6,12));
   title(sprintf('sf xc control, mean=%.3f',mean(aa)));

   subplot(3,2,3);
   aa=sum(ormatchxc(1,2:3,ftshiftidx)-...
           ordiffxc(1,2:3,ftshiftidx),2);
      %./(2-sum(ordiffxc(1,2:3,ftshiftidx),2));
   hist(aa,linspace(-0.6,0.6,12))
   title(sprintf('or xc shift, mean=%.3f',mean(aa)))
   subplot(3,2,4);
   aa=sum(ormatchxc(2,2:3,ftshiftidx)-ordiffxc(2,2:3,ftshiftidx),2);
      %./(2-sum(ordiffxc(2,2:3,ftshiftidx),2)) ;
   hist(aa,linspace(-0.6,0.6,12))
   title(sprintf('or xc control, mean=%.3f',mean(aa)))
   
   strfranka=cat(4,res.strfranka);
   sara=cat(3,squeeze(([strfranka(1,2,:,:)-strfranka(1,3,:,:)])), ...
            squeeze(([strfranka(2,3,:,:)-strfranka(2,2,:,:)])));
   sm=squeeze(mean(mean(sara,3)>0,1))';
   se=squeeze(std(mean(sara,3)>0,0,1))'*sqrt(size(sara,1)-1);
   %sm(sm<0.5 & sm>0)=0;
   
   aa=squeeze(sum(hmatchxc(1,2:3,:)-hdiffxc(1,2:3,:),2));
       %./ (2-sum(hdiffxc(1,2:3,:),2));

   sm=sm(ftshiftidx);
   se=se(ftshiftidx);
   aa=aa(ftshiftidx);
   
   xx=linspace(-0.6,0.6,10);
   nxc=zeros(length(xx),3);
   catidx=zeros(length(ftshiftidx),1);
   catidx(sm>0.5 & sm>2.*se)=1;
   catidx(sm<=2.*se)=3;
   catidx(sm<0.5 & sm<2.*se)=2;
   aa(catidx==1 & aa<0)=0.05;
   
   for jj=1:3,
      nxc(:,jj)=hist(aa(find(catidx==jj)),xx);
   end
   
   subplot(3,2,5);
   %hist(aa,linspace(-0.6,0.6,10))
   bar(xx,nxc,'stacked')
   hold on
   aaa=axis;
   plot([0 0],[aaa(3) aaa(4)],'k--');
   hold off
   axis([xx(1)-0.1 xx(end)+0.1 0 aaa(4)+1]);
   title(sprintf('h xc shift, mean=%.3f',mean(aa)))
   
   subplot(3,2,6);
   aa=sum(hmatchxc(2,2:3,ftshiftidx)-hdiffxc(2,2:3,ftshiftidx),2);
       %./ (2-sum(hdiffxc(2,2:3,ftshiftidx),2));
   for jj=1:3,
      nxc(:,jj)=hist(aa(find(catidx==jj)),xx);
   end
   bar(xx,nxc,'stacked')
   hold on
   aaa=axis;
   plot([0 0],[aaa(3) aaa(4)],'k--');
   hold off
   axis([xx(1)-0.1 xx(end)+0.1 0 aaa(4)+1]);
   %hist(aa,linspace(-0.6,0.6,12))
   title(sprintf('h xc control, mean=%.3f',mean(aa)))
   
   fullpage portrait
   
   for fs=1:2,
      
      [mean(mean(sfmatchxc(fs,2:3,ftshiftidx)-sfdiffxc(fs,2:3,ftshiftidx),2))...
       mean(mean(sfmatchxc(fs,2:3,ftshiftidx),2))...
       mean(mean(sfdiffxc(fs,2:3,ftshiftidx),2))...
       mean(mean(ormatchxc(fs,2:3,ftshiftidx)-ordiffxc(fs,2:3,ftshiftidx),2))...
       mean(mean(ormatchxc(fs,2:3,ftshiftidx),2))...
       mean(mean(ordiffxc(fs,2:3,ftshiftidx),2))]
   end
   
   aa=sum(hmatchxc(1,2:3,:)-hdiffxc(1,2:3,:),2); 
   squeeze(aa([12 66 80 83]))
   
   fullxcfrac=squeeze(fullxc([1 3 4],2,:))';
   fullxcfrac=diff([zeros(size(fullxcfrac,1),1) fullxcfrac].^2,[],2)./...
       repmat(fullxcfrac(:,end).^2,[1 3]);
   
   [saa,sid]=sort(-aa(ftshiftidx));
   
   bigsrf=[];
   bigtar=[];
   bigpatch=[];
   for ii=1:length(ftshiftidx),
      jj=ftshiftidx(sid(ii));
      %fprintf('%s: xc0=%.3f xcdcg=%.3f xctun=%.3f tsi=%.3f\n',...
      %        res(jj).cellid,fullxc([1 3 4],2,jj),aa(jj));
      if catidx(sid(ii))<3,
         fprintf('*');
      end
      fprintf('%s: xc0=%.3f xcdcg=%.3f xctun=%.3f tsi=%.3f\n',...
              res(jj).cellid,fullxcfrac(jj,:),aa(jj));
      
      bigsrf=cat(3,bigsrf,squeeze(res(jj).Hset(:,1,2:3)));
      bigtar=cat(3,bigtar,res(jj).mpatches(:,2:3)-...
                 repmat(res(jj).mS(:,1),[1 2]));
      bigpatch=cat(3,bigpatch,reshape(res(jj).bigpatches,32*32,2));
   end
   
   bigsrf=cat(2,bigsrf,bigsrf(:,1,:)-bigsrf(:,2,:));
   
   figure(1);
   clf
   showkern(bigsrf,'pfftgr');
   fullpage portrait
   figure(2);
   clf
   showkern(bigtar,'pfftgr');
   fullpage portrait
   figure(3);
   clf
   showkern(bigpatch-100,'space');
   colormap(gray);
   fullpage portrait
   
end






