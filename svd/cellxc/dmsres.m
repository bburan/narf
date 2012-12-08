% function res=dmsres(cellid,batch,showextra)
%
% load results from kernfile in sRunData and display attentional
% modulation info ... kernfile generated from kerncomp4
%
% res.resp(targ,attcond): actual response to each target for
%                         attcond=(baseline, A, B)
%                         targ is currently a or b
% res.rank(targ,attcond): response in units of s.d. from mean
% res.mresp(attcond):     mean response in each attcond
% res.strfresp(targ,attcond): pred response to each target for
%                             attcond=(baseline, A, B)
% res.strfrank(targ,attcond): pred response in units of s.d. from mean
% res.strfmean(attcond):      fit mean response in each attcond
% res.predxc(controlstat,att): pred corr for att model(controlstat=1) 
%                              or rand att (controlstat=2)
% res.pxc(att):   p value associated with each att cond
%
% DUMP SOME INTERESTING KERNELS TO EPS:
%
%cellids={'v0102','v0116','v0121','v0137','v0146','v0147','v0164'};
%for ii=1:length(cellids),
%   dmsres(cellids{ii},99);
%   print('-f2','-depsc',[cellids{ii},'.kern.eps']);
%end
%
% created SVD 6/24/04 - hacked from kvares.m
%
function res=dmsres(cellid,batch,showextra)

if ~exist('showextra','var'),
   showextra='';
end

dbopen;

if exist('batch','var') & batch>0,
   fprintf('dmsres.m: (cell %s bat %d)\n',cellid,batch);
   goodbatch=zeros(1,length(batch));
   batchcount=0;
   for ii=1:length(batch),
      sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
           ' AND batch=',num2str(batch(ii))];
      trd=mysql(sql);
      if ~isempty(trd),
         goodbatch(ii)=1;
         batchcount=batchcount+1;
         rundata(batchcount)=trd;
      end
   end
   batch=batch(find(goodbatch));
   resfile=[rundata(1).respath,rundata(1).kernfile,'.gz'];
else
   % no batch specified, assume it's a file.
   resfile=cellid;
   if ~exist(resfile,'file'),
      batchcount=0;
   else
      batchcount=1;
   end
end


if batchcount==0,
   disp('no entries found in db!');
   if nargout>0,
      res=[];
   end
   return
end

global GCOLORMAP
GCOLORMAP=redblue;

rcsetstrings;

fprintf('Loading: %s\n',resfile);
zload(resfile);

cellid=params.cellid;
batch=params.batch;


% flag for format to display kernels. if showSFOR=1, show in
% spatial freq/orientation map. otherwise, show in native kernfmt
if strcmp(params.kernfmt,'space'),
   showSFOR=0;
else
   showSFOR=1;
end

if ~isfield(params,'altkva') | ~strcmp(params.altkva,'xcdms'),
   disp('need data analyzed by xcdms to work.');
   return
end
%if batch==82 & ~strcmp(params.kernfmt,'pfft+4'),
%   params.kernfmt='pfft+4';
%   params.kernfmt
%end

r=[];
for fidx=1:length(params.respfiles),
   
   [tr,picid,frameidx]=resploaddmsmatch(params.respfiles{fidx},...
                                        params.resploadparms{:});
   r=cat(1,r,squeeze(respvarbins(tr,params.respfilterparms{:})));
end

% raw kernels
H=zeros([size(vstrf(1).h,1),params.bootcount,attcount]);
mS=zeros([size(vstrf(1).mS,1),params.bootcount,attcount]);
nlparms=zeros(2,params.bootcount,attcount);
for attidx=1:attcount,
   H(:,:,attidx)=cat(2,vstrf(attidx,1,1,:).h);
   mS(:,:,attidx)=cat(2,vstrf(attidx,1,1,:).mS);
   nlparms(:,:,attidx)=cat(2,vstrf(attidx,1,1,:).nlparms);
   
   %for bootidx=1:params.bootcount,
   %   nlparms(1,bootidx,attidx)=nlparms(1,bootidx,attidx)-...
   %       mS(:,bootidx,attidx)'*H(:,bootidx,attidx);
   %   mS(:,bootidx,attidx)=0;
   %end
end

fprintf('sfsidx(%d)/sigidx(%d)=',params.sfscount,length(params.sffiltsigma));
for ii=1:params.bootcount,
   fprintf('%d/%d ',vstrf(1,1,1,ii).parms.sfsfit,vstrf(1,1,1,ii).parms.sigfit);
   if mod(ii,10)==0 | ii==params.bootcount,
      fprintf('\n');
   end
end
Hall=H; % save jackknifed kernels for display

% scale dc term to Hz
globaldc=globaldc.*60;

% H2 translate raw kernels into local kernels for Aa and Bb
H2=repmat(H(:,:,1),[1 1 3]);
H2(:,:,2)=H(:,:,1).*(1+repmat(globalgain(2,:),[spacecount 1])) + ...
          H(:,:,2);
H2(:,:,3)=H(:,:,1).*(1-repmat(globalgain(2,:),[spacecount 1])) - ...
          H(:,:,2);

Hm=H2(:,:,1);
Hd=H2(:,:,2)-H2(:,:,3);
DCd=globaldc(2,:);

H3=repmat(H(:,:,1),[1 1 3]);
H3(:,:,2)=H(:,:,1).*(1+repmat(globalgain(3,:),[spacecount 1])) + ...
          H(:,:,3);
H3(:,:,3)=H(:,:,1).*(1-repmat(globalgain(3,:),[spacecount 1])) - ...
          H(:,:,3);

% feature kernels for space in only
H4=repmat(H3(:,:,2),[1 1 3]);
H4(:,:,2)=H3(:,:,2).*(1+repmat(globalgain(4,:),[spacecount 1])) + ...
          H(:,:,4);
H4(:,:,3)=H3(:,:,2).*(1-repmat(globalgain(4,:),[spacecount 1])) - ...
          H(:,:,4);

% feature kernels for space out
H5=H2.*2-H4;

nlparms2a=zeros(size(nlparms,1),params.bootcount,3);
nlparms2a(:,:,1)=nlparms(:,:,1);
nlparms2a(:,:,2)=nlparms(:,:,1)+nlparms(:,:,2);
nlparms2a(:,:,3)=nlparms(:,:,1)-nlparms(:,:,2);
nlparms3a=zeros(size(nlparms,1),params.bootcount,3);
nlparms3a(:,:,1)=nlparms(:,:,1);
nlparms3a(:,:,2)=nlparms(:,:,1)+nlparms(:,:,3);
nlparms3a(:,:,3)=nlparms(:,:,1)-nlparms(:,:,3);
nlparms4a=zeros(size(nlparms,1),params.bootcount,3);
nlparms4a(:,:,1)=nlparms3a(:,:,1);
nlparms4a(:,:,2)=nlparms3a(:,:,1)+nlparms(:,:,4);
nlparms4a(:,:,3)=nlparms3a(:,:,1)-nlparms(:,:,4);


% how similar each target is to the strf
targpredsa=zeros(size(mpatches,2)-1,size(H2,3),size(H2,2));
difftargpredsa=zeros(size(mpatches,2)-1,size(H,3),size(H,2));

% how similar each stimulus is to the strf
targsima=zeros(size(mpatches,2)-1,size(H2,3),size(H2,2));
strfsima=zeros(blen,size(H2,3),size(H2,2));
diffstrfsima=zeros(blen,size(H,3),size(H,2));

targprojHm=zeros(size(mpatches,2)-1,size(Hm,2));
targprojHd=zeros(size(mpatches,2)-1,size(Hd,2));
stimprojHm=zeros(blen,size(Hm,2));
stimprojHd=zeros(blen,size(Hd,2));

stimlen=size(bstim,1);
chancount=size(bstim,2);
normstim=bstim'-repmat(mean(mS(:,:),2),[1 stimlen]);
normtarg=mpatches(:,2:end)-repmat(mean(mS(:,:),2),[1 size(mpatches,2)-1]);

%sS=std(normstim,0,2);
%normstim=normstim./repmat(sS,[1 stimlen]);
%normtarg=normtarg./repmat(sS,[1 size(mpatches,2)-1]);

%normstim=normstim-repmat(mean(normstim),[chancount 1]);
%normtarg=normtarg-repmat(mean(normtarg),[chancount 1]);
%normstim=normstim./repmat(std(normstim),[chancount 1]);
%normtarg=normtarg./repmat(std(normtarg),[chancount 1]);

for ii=1:size(H2,2),  % for each jackknife
   targpredsa(:,:,ii)=(squeeze(H2(:,ii,:))'* ...
                       (mpatches(:,2:end)-...
                        repmat(mS(:,ii,1),1,size(mpatches,2)-1)))' + ...
       squeeze(repmat(nlparms2a(1,ii,:),[size(mpatches,2)-1 1]));
   difftargpredsa(:,:,ii)=(squeeze(H(:,ii,:))'* ...
                       (mpatches(:,2:end)-...
                        repmat(mS(:,ii,1),1,size(mpatches,2)-1)))' + ...
       squeeze(repmat(nlparms(1,ii,:),[size(mpatches,2)-1 1]));
   
   normH=squeeze(H2(:,ii,:));
   %normH=normH-repmat(mean(normH),[size(normH,1) 1]);
   %normH=normH./repmat(sS,[1 size(normH,2)]);
   
   targsima(:,:,ii)=(normH' * normtarg)';
   strfsima(:,:,ii)=(normH' * normstim)';
   % ... + squeeze(repmat(nlparms2a(1,ii,:),[blen 1]));
   diffstrfsima(:,:,ii)=(squeeze(H(:,ii,:))' * normstim)';
   
   targprojHm(:,ii)=Hm(:,ii)'*normtarg;
   targprojHd(:,ii)=Hd(:,ii)'*normtarg;
   stimprojHm(:,ii)=Hm(:,ii)'*normstim;
   stimprojHd(:,ii)=Hd(:,ii)'*normstim;
end

% find outlier SRFs
dev=squeeze((targpredsa(1,1,:)-mean(targpredsa(1,1,:))) + ...
            (targpredsa(2,1,:)-mean(targpredsa(2,1,:))));
keepidx=setdiff(1:params.bootcount,...
                [min(find(dev==min(dev))) min(find(dev==max(dev)))]);


% average over non-outlier SRFs for display
targpreds=mean(targpredsa(:,:,keepidx),3);
difftargpreds=mean(difftargpredsa(:,:,keepidx),3);

targsim=mean(targsima(:,:,keepidx),3);
strfsim=mean(strfsima(:,:,keepidx),3);
diffstrfsim=mean(diffstrfsima(:,:,keepidx),3);


targprojHm=mean(targprojHm(:,keepidx),2);
targprojHd=mean(targprojHd(:,keepidx),2);
stimprojHm=mean(stimprojHm(:,keepidx),2);
stimprojHd=mean(stimprojHd(:,keepidx),2);

if 1,
   figure(3);
   clf
   %keyboard
   tempH=cat(3,squeeze(mean(H4(:,:,2:3),2)),...
             squeeze(mean(H5(:,:,2:3),2)));
   tempH=cat(2,tempH,-diff(tempH,[],2));
   
   showkern(tempH,'pfftgr');
elseif 0,
   figure(3);
   clf
   showkern(cat(3,Hm(:,1:2),Hd(:,1:2)),'pfftgr');
   subplot(1,3,3);
   cla
   sampidx=round(linspace(1,length(stimprojHm),200));
   plot(stimprojHm(sampidx),stimprojHd(sampidx),'.')
   hold on
   plot(targprojHm(1),targprojHd(1),'rx');
   plot(targprojHm(2),targprojHd(2),'kx');
   text(targprojHm(1),targprojHd(1),'A');
   text(targprojHm(2),targprojHd(2),'B');
   hold off
   axis equal
   axis square
end

shrf=0.5;
mH=mean(H,2);
eH=std(H,0,2).*sqrt(size(H,2)-1);
H=shrinkage(mH,eH,shrf);
mH=mean(H2,2);
eH=std(H2,0,2).*sqrt(size(H2,2)-1);
H2=shrinkage(mH,eH,shrf);
mH=mean(H3,2);
eH=std(H3,0,2).*sqrt(size(H3,2)-1);
H3=shrinkage(mH,eH,shrf);
mH=mean(H4,2);
eH=std(H4,0,2).*sqrt(size(H4,2)-1);
H4=shrinkage(mH,eH,shrf);
mS=mean(mS,2);

%H=(mean(H(:,keepidx,:),2));
%H2=(mean(H2(:,keepidx,:),2));
%H3=(mean(H3(:,keepidx,:),2));
%H4=(mean(H4(:,keepidx,:),2));
%mS=(mean(mS(:,keepidx,:),2));
nlparms=(mean(nlparms(:,keepidx,:),2));
nlparms2=(mean(nlparms2a(:,keepidx,:),2));
nlparms3=(mean(nlparms3a(:,keepidx,:),2));
nlparms4=(mean(nlparms4a(:,keepidx,:),2));

pstrfsim=((H2(:,1).*(H2(:,1)>0))'* (bstim'-repmat(mS(:,1),1,blen)))';
nstrfsim=((H2(:,1).*(H2(:,1)<0))'* (bstim'-repmat(mS(:,1),1,blen)))';

[r2,picid,frameidx,targmatch]=...
    resploaddms(params.respfiles{1},'',0,1,{'A','B','a','b'});
r2=squeeze(respvarbins(r2,params.respfilterparms{:}));
t1Aidx=find(targmatch'==1 & ~isnan(r2(:,1)));
t1Bidx=find(targmatch'==1 & ~isnan(r2(:,2)));
t1idx=[t1Aidx(1:min([length(t1Aidx) ]));  %20
       t1Bidx(1:min([length(t1Bidx) ]))];
t2Aidx=find(targmatch'==2 & ~isnan(r2(:,1)));
t2Bidx=find(targmatch'==2 & ~isnan(r2(:,2)));
t2idx=[t2Aidx(1:min([length(t2Aidx) ]));
       t2Bidx(1:min([length(t2Bidx) ]))];

% record response to each target in each attention state
% res.resp(targetnum,attstate)
clear res
res.cellid=cellid;
res.kernfile=[rundata(1).respath,rundata(1).kernfile];
res.resp(1,1)=nanmean([r2(t1idx,1);r2(t1idx,2)]);
res.resp(2,1)=nanmean([r2(t2idx,1);r2(t2idx,2)]);
res.resperr(1,1)=nanstd([r2(t1idx,1);r2(t1idx,2)])./sqrt(length(t1idx).*2);
res.resperr(2,1)=nanstd([r2(t2idx,1);r2(t2idx,2)])./sqrt(length(t1idx).*2);
res.mresp(1)=nanmean(nanmean(r2(find(targmatch==0),:)));
res.resp(1,2)=nanmean(r2(t1idx,1));
res.resp(2,2)=nanmean(r2(t2idx,1));
res.resperr(1,2)=nanstd(r2(t1idx,1))./sqrt(length(t1idx));
res.resperr(2,2)=nanstd(r2(t2idx,1))./sqrt(length(t1idx));
res.mresp(2)=nanmean(r2(find(targmatch==0),1));
res.resp(1,3)=nanmean(r2(t1idx,2));
res.resp(2,3)=nanmean(r2(t2idx,2));
res.resperr(1,3)=nanstd(r2(t1idx,2))./sqrt(length(t1idx));
res.resperr(2,3)=nanstd(r2(t2idx,2))./sqrt(length(t1idx));
res.mresp(3)=nanmean(r2(find(targmatch==0),2));

% figure out rank of these responses wrt response to all distractors.
tt=sort(([r2(find(targmatch'==0 & ~isnan(r2(:,1))),1); ...
         r2(find(targmatch'==0 & ~isnan(r2(:,2))),2);
         res.resp(:,2);res.resp(:,3)]));
tt=mean(tt);
%tt=tt(round(length(tt(:)).*0.95));
%res.rank(1,1)=max(find(res.resp(1,1)>tt))./length(tt);
%res.rank(2,1)=max(find(res.resp(2,1)>tt))./length(tt);
%res.rank(:,1)=(res.resp(:,1)-mean(tt))./std(tt);
res.rank(:,1)=res.resp(:,1)./max(tt(:));

% indices are (rank,targ)
tt=sort([r2(find(targmatch'==0 & ~isnan(r2(:,1))),1);
         res.resp(:,2)]);
tt=mean(tt);
%tt=tt(round(length(tt(:)).*0.95));
%res.rank(1,2)=max(find(res.resp(1,2)>tt))./length(tt);
%res.rank(2,2)=max(find(res.resp(2,2)>tt))./length(tt);
%res.rank(:,2)=(res.resp(:,2)-mean(tt))./std(tt);
res.rank(:,2)=res.resp(:,2)./max(tt(:));

tt=sort([r2(find(targmatch'==0 & ~isnan(r2(:,2))),2);
         res.resp(:,3)]);
tt=mean(tt);
%tt=tt(round(length(tt(:)).*0.95));
%res.rank(1,3)=max(find(res.resp(1,3)>tt))./length(tt);
%res.rank(2,3)=max(find(res.resp(2,3)>tt))./length(tt);
%res.rank(:,3)=(res.resp(:,3)-mean(tt))./std(tt);
res.rank(:,3)=res.resp(:,3)./max(tt(:));


% predicted baseline response to targets, mean response, rank of
% targets (rank calculated within attention condition. essentially
% a measure of dot product between stimulus and SRF
res.strfresp=targpreds;
res.strfmean=mean(strfsim);
tt=strfsim(:,[2 3]);

targpredsnorm=(H2(:,:)'* ...
    (mpatches(:,2:end)-repmat(mS(:,1),1,size(mpatches,2)-1)))';
for aa=1:3,
   tt=sort(strfsim(:,aa));
   %tt=tt.*60.*(tt>0)+1;
   %tt=sort([random('poiss',tt);random('poiss',tt);...
   %         random('poiss',tt);random('poiss',tt);...
   %         random('poiss',tt);random('poiss',tt)]-1)./60;
   %res.strfrank(1,aa)=max([0; find(targpreds(1,aa)>tt)])./length(tt);
   %res.strfrank(2,aa)=max([0; find(targpreds(2,aa)>tt)])./length(tt);
   %res.strfrank(:,aa)=(targpreds(:,aa)-mean(tt))./std(tt);
   res.strfrank(:,aa)=targsim(:,aa)./std(tt(:));
   %res.strfrank(:,aa)=targpredsnorm(:,aa)./norm(H2(:,aa));
   %res.strfrank(1,aa)=res.strfrank(1,aa)./...
   %    norm(mpatches(:,2)-mS(:,1));
   %res.strfrank(2,aa)=res.strfrank(2,aa)./...
   %    norm(mpatches(:,3)-mS(:,1));

end

% "distance" between targets A and B, normalized by mean of all
% stim.  maybe std of stim is better norm factor?  i guess it scales
res.targdist=mean(abs(mpatches(:,2)-mpatches(:,3))); %./mS(:,1));
%fprintf('targdist=%.2f\n',res.targdist);

res.globaldc=globaldc;
res.globalgain=globalgain;
res.baserange=[min(bresp(:,1)) max(bresp(:,1))];
if exist('globaldconly','var'),
   res.globaldconly=globaldconly;
else
   res.globaldconly=zeros(size(globaldc));
end

% prediction accuracy stuff.
res.predxc=predxc;
res.pxc=pxc;
res.randxc=randxc;

% predict WHOLE response with kernels, rather than just diffs
%
% these are now generated in xcdms.m:
%  fullxc  - xc normalized for addition of new att terms to model
%  fullp   - whether each additional term signif improves predictions
% rows: 1. no att over 0 (always sig), 2. dc over no att
%       3. dcg over dc, 4. local over dcg
% cols: 1. no att, 2. feature Aa vs Bb, 3. spatial AB vs ab

res.fullxc=fullxc;
res.fullp=fullp;

fprintf('             resp to:\n');
fprintf('atten.:     1    (rnk)     2     (rnk)   mean\n');
fprintf('both:     %5.2f (%5.2f)  %5.2f (%5.2f)  %.2f\n',...
        res.resp(1,1)*60,res.rank(1,1),...
        res.resp(2,1)*60,res.rank(2,1),res.mresp(1)*60);
fprintf('targ1:    %5.2f (%5.2f)  %5.2f (%5.2f)  %.2f\n',...
        res.resp(1,2)*60,res.rank(1,2),...
        res.resp(2,2)*60,res.rank(2,2),res.mresp(2)*60);
fprintf('targ2:    %5.2f (%5.2f)  %5.2f (%5.2f)  %.2f\n',...
        res.resp(1,3)*60,res.rank(1,3),...
        res.resp(2,3)*60,res.rank(2,3),res.mresp(3)*60);
fprintf('strfb:    %5.2f (%5.2f)  %5.2f (%5.2f)  %.2f\n',...
        res.strfresp(1,1)*60,res.strfrank(1,1),...
        res.strfresp(2,1)*60,res.strfrank(2,1),res.strfmean(1)*60);
fprintf('strf1:    %5.2f (%5.2f)  %5.2f (%5.2f)  %.2f\n',...
        res.strfresp(1,2)*60,res.strfrank(1,2),...
        res.strfresp(2,2)*60,res.strfrank(2,2),res.strfmean(2)*60);
fprintf('strf2:    %5.2f (%5.2f)  %5.2f (%5.2f)  %.2f\n',...
        res.strfresp(1,3)*60,res.strfrank(1,3),...
        res.strfresp(2,3)*60,res.strfrank(2,3),res.strfmean(3)*60);

obincount=16;
sfbincount=8;
if strcmp(showextra,'curves'),
   Hset=cat(2,H2,H3);
   acount=size(Hset,2);
   hcount=size(Hset,3);
   
   % show marginal tuning curves
   disp('computing tuning curves');
   ortuning=zeros(obincount,acount,hcount);
   sftuning=zeros(sfbincount,acount,hcount);
   orslice=zeros(obincount,acount,hcount);
   sfslice=zeros(sfbincount,acount,hcount);
   seprat=zeros(acount,hcount);
   orpeak=zeros(acount,hcount);
   orbw=zeros(acount,hcount);
   sfpeak=zeros(acount,hcount);
   sfbw=zeros(acount,hcount);
   
   for aidx=1:acount,
      
      % figure out mean for fixing svd output
      tsf0=pfft2sf(mean(Hset(:,1,1),2));
      tsf0=sf2gr(tsf0,obincount,sfbincount);
      mor=mean(tsf0,2);
      
      for hidx=1:hcount,
         tsfIR=pfft2sf(Hset(:,aidx,hidx),params.kernfmt);
         [tsfgr,obins,sfbins]=sf2gr(tsfIR,obincount,sfbincount);
         psfgr=tsfgr.*(tsfgr>0);
         
         [u,s,v]=svd(tsfgr);
         if sum(diag(s))>0,
            seprat(aidx,hidx)=s(1)./sum(diag(s));
         end
         if 0,
            por=u(:,1).*s(1);
            psf=v(:,1).*s(1);
            if por'*mor<0,
               por=-por;
               psf=-psf;
            end
            mor=mor+por;
         else
            por=mean(tsfgr,2);
            psf=mean(tsfgr',2);
         end
         
         peaksfidx=min(find(psf(:,1)==max(psf(:,1))));
         orslice(:,aidx,hidx)=tsfgr(:,peaksfidx);
         peakoridx=min(find(por(:,1)==max(por(:,1))));
         sfslice(:,aidx,hidx)=tsfgr(peakoridx,:)';
         
         % adjust circstats to deal with extra pi half of circle
         [orpeak(aidx,hidx),orbw(aidx,hidx)]=circstats(por(:,1));
         orbw(aidx,hidx)=orbw(aidx,hidx)./2 .*180/pi.*2; 
         orpeak(aidx,hidx)=mod(orpeak(aidx,hidx) ./2 .*180/pi,180);
         
         %beta=fitgauss1d(log2(sfbins),psf(:,1));
         
         % peak sf tuning is mean of sf tuning curve
         %sfpeak(aidx,hidx)=2.^beta(1);
         
         % sf bw defined by de valois is (log2 H2 - log2 L2)
         %sfbw(aidx,hidx)=(beta(2)+beta(2)) .* sqrt(2.*log(2));
         
         ortuning(:,aidx,hidx)=por(:,1);
         sftuning(:,aidx,hidx)=psf(:,1);
      end
   end
   res.Hset=Hset;
   res.obins=obins;
   res.sfbins=sfbins;
   res.ortuning=ortuning;
   res.sftuning=sftuning;
   res.sfpeak=sfpeak;
   res.sfbw=sfbw;
   res.orpeak=orpeak;
   res.orbw=orbw;
end

if 0,
   rr=[rprec(:,1) origresp(:,1:4)];
   rr=sortrows(rr(find(~isnan(sum(rr,2))),:));
   pp=rr(:,1);
   rr=rr(:,2:end);
   
   hingeparms=zeros(3,size(rr,2)-1);
   sigparms=zeros(4,2,size(rr,2)-1);
   sigrange=zeros(2,2,size(rr,2)-1);
   
   if 1,
      disp('fitting hinges!');
      hingeparms(:,1)=fithinge(rr(:,1),rr(:,2))';
      hingeparms(:,2)=fithinge(rr(:,1),rr(:,3))';
      hingeparms(:,3)=fithinge(rr(:,1),rr(:,4))';
   else
      % fit sigmoids!
      disp('fitting sigmoids!');
      
      for attidx=2:3,
         for segidx=1:2,
            pp=prednoatt((1:blen)+blen*(segidx-1),attidx);
            rr=respatt((1:blen)+blen*(segidx-1),attidx);
            gidx=find(~isnan(pp+rr));
         
            sigparms(:,segidx,attidx-1)=fitsigmoid(pp(gidx),rr(gidx))';
            sigrange(:,segidx,attidx-1)=[mean(pp)-3*std(pp); mean(pp)+3*std(pp)];
         end
      end
      %p0=linspace(sigrange(1,1,1),sigrange(2,1,1));
      %clf
      %plot(p0',[sigmoid(sigparms(:,1,1),p0') sigmoid(sigparms(:,2,1),p0')]);
      %drawnow
   end
   
   % define base range to span responses sampled on more than 1 trial
   res.baserange=[min(rr(:,1)) max(rr(:,1))];
   
   if 1,
      bincount=16;
      ppp=zeros(bincount,1);
      rrr=zeros(bincount,4);
      rrb=round(linspace(1,size(rr,1),bincount+1));
      for bb=1:bincount,
         ppp(bb)=mean(pp(rrb(bb):(rrb(bb+1)-1)));
         rrr(bb,:)=mean(rr(rrb(bb):(rrb(bb+1)-1),:));
      end
   else
      % average rows with same baseline response
      while (sum(diff(rr(:,1))==0)>0),
         sidx=min(find(diff(rr(:,1))==0));
         eidx=max(find(rr(:,1)==rr(sidx,1)));
         rr=[rr(1:sidx-1,:); mean(rr(sidx:eidx,:)); rr(eidx+1:end,:)];
      end
      sidx=max([1          max(find(rr(:,1)<mean(rr(:,1))-std(rr(:,1))*3))]);
      eidx=min([size(rr,1) min(find(rr(:,1)>mean(rr(:,1))+std(rr(:,1))*3))]);
      rr=rr(sidx:eidx,:);
      
      rrr=zeros(250,4);
      rrr(:,1)=linspace(res.baserange(1),res.baserange(2),size(rrr,1))';
      smrat=size(rr,1)./size(rrr,1)./2;
      
      if smrat<1,
         rrr(:,2)=interp1(rr(:,1),rr(:,2),rrr(:,1),'linear');
         rrr(:,3)=interp1(rr(:,1),rr(:,3),rrr(:,1),'linear');
         rrr(:,4)=interp1(rr(:,1),rr(:,4),rrr(:,1),'linear');
      else
         %rrr(:,2)=interp1(rr(:,1),gsmooth(rr(:,2),smrat),rrr(:,1),'linear');
         %rrr(:,3)=interp1(rr(:,1),gsmooth(rr(:,3),smrat),rrr(:,1),'linear');
         rrr(:,2)=interp1(rr(:,1),gsmooth(rr(:,2),smrat),rrr(:,1));
         rrr(:,3)=interp1(rr(:,1),gsmooth(rr(:,3),smrat),rrr(:,1));
         rrr(:,4)=interp1(rr(:,1),gsmooth(rr(:,4),smrat),rrr(:,1));
      end
   end
   
   clf
   plot(ppp,rrr);
   drawnow 
   
   res.pset=ppp;
   res.rset=rrr;
   
   res.hingeparms=hingeparms;
   res.sigparms=sigparms;
   res.sigrange=sigrange;
else
   % spit out responses under each discrete attention condition
   % binned according to linear prediction
   
   % "baseline prediction" is prediction of ABab STRF
   rr=[rprec(:,1) r(:,1:4)];
   rr=sortrows(rr(find(~isnan(sum(rr,2))),:));
   pp=rr(:,1);
   rr=rr(:,2:end);
   
   bincount=30;
   ppp=zeros(bincount,1);
   rrr=zeros(bincount,4);
   rrb=round(linspace(1,size(rr,1),bincount+1));
   for bb=1:bincount,
      ppp(bb)=mean(pp(rrb(bb):(rrb(bb+1)-1)));
      rrr(bb,1:4)=mean(rr(rrb(bb):(rrb(bb+1)-1),:));
   end
   
   % rset columns in resploaddmsmatch order: 'A' 'B' 'a' 'b'
   res.pset=ppp;
   res.rset=rrr;
   
end

if nargout>0,
   disp('nargout>0: skipping figures');
   return
end

SHOWALLKERNS=0;
if SHOWALLKERNS,
   % show jackknifed kernels
   figure(1);
   showkern(Hall,params.kernfmt);
   
   kshowcount=min([params.bootcount 12]);
   for ii=1:kshowcount,
      for attidx=1:attcount,
         subplot(attcount,kshowcount,(attidx-1).*kshowcount+ii),
         
         if ii==1,
            title(sprintf('%s %d/%d',cellid,...
                          vstrf(attidx,1,1,ii).parms.sfsfit,...
                          vstrf(attidx,1,1,ii).parms.sigfit));
         else
            title(sprintf('%d/%d',vstrf(attidx,1,1,ii).parms.sfsfit,...
                          vstrf(attidx,1,1,ii).parms.sigfit));
         end
      end
   end
   %keyboard
   
else
   
   MAXPTS=300;
   if size(bresp,1)>MAXPTS,
      ptidx=round(linspace(1,size(bresp,1),MAXPTS));
   else
      ptidx=1:size(bresp,1);
   end
   
   useidx=find(~isnan(r(:,1)) & ~isnan(r(:,2)));
   ptidx=useidx;
   rAB=mean(r(useidx,[1 2]),2).*60;
   rab=mean(r(useidx,[3 4]),2).*60;
   rABdiff=(r(useidx,1)-r(useidx,2))./2.*60;
   rA=r(ptidx,1).*60;
   rB=r(ptidx,2).*60;
   
   rAa=(r(ptidx,1)+r(ptidx,3))./2.*60;
   rBb=(r(ptidx,2)+r(ptidx,4))./2.*60;
   %rA=rAa;
   %rB=rBb;
   NN=10;
   pall=rprec(ptidx,1).*60;
   rall=bresp(ptidx,1).*60;
   pfdiff=rprec(ptidx,2).*60;
   rfdiff=bresp(ptidx,2).*60;
   spall=sort(rAB(~isnan(rAB+rab))+rab(~isnan(rAB+rab)))./2;
   spall=[spall(:); spall(end)+1];
   epall=spall(round(linspace(1,length(spall),NN+1)));
   pbinned=zeros(NN,1);
   rbinned=zeros(NN,4);
   rbinnederr=zeros(NN,4);
   for nn=1:NN,
      ff=find((rAB+rab)./2>=epall(nn) & (rAB+rab)./2<epall(nn+1));
      pbinned(nn)=nanmean(rAB(ff));
      rbinned(nn,1)=nanmean(rAa(ff));
      rbinnederr(nn,1)=nanstd(rAa(ff))./sqrt(length(ff));
      rbinned(nn,2)=nanmean(rBb(ff));
      rbinnederr(nn,2)=nanstd(rBb(ff))./sqrt(length(ff));
      rbinned(nn,3)=nanmean(rAB(ff));
      rbinnederr(nn,3)=nanstd(rAB(ff))./sqrt(length(ff));
      rbinned(nn,4)=nanmean(rab(ff));
      rbinnederr(nn,4)=nanstd(rab(ff))./sqrt(length(ff));
   end

   
   % show some basic stats of responses in different conditions
   figure(1);
   clf
   
   rmin=0;
   %rmax=nanmax([r(:);bresp(:,1);bresp(:,2)]).*60 .*2/3;
   rmax=nanmax([bresp(:,1)]).*60;
   
   subplot(2,3,1);
   
   %r=sqrt(r);
   %res.resp=sqrt(res.resp);
   %res.resperr=sqrt(res.resperr);
   
   scatter(rBb,rAa,'.');
   hold on
   if 1,
      %plot((res.resp(1,2)+res.resperr(1,2).*[-1 1]).*60,...
      %     res.resp(1,3).*[1 1].*60,'r-');
      %plot(res.resp(1,2)*[1 1].*60,...
      %     (res.resp(1,3)+res.resperr(1,3).*[-1 1]).*60,'r-');
      %plot((res.resp(2,2)+res.resperr(2,2).*[-1 1]).*60,...
      %     res.resp(2,3).*[1 1].*60,'k-');
      %plot(res.resp(2,2)*[1 1].*60,...
      %     (res.resp(2,3)+res.resperr(2,3).*[-1 1]).*60,'k-');
   else
      plot(res.resp(1,3).*60,r2(t1idx,1).*60,'r.');
      plot(r2(t1idx,2).*60,res.resp(1,2).*60,'r.');
      plot(res.resp(2,3).*60,r2(t2idx,1).*60,'k.');
      plot(r2(t2idx,2).*60,res.resp(2,2).*60,'k.');
   end
   plot([rmin rmax],[rmin rmax],'k--');
   
   p=polyfit(rAB,rABdiff,1);
   m=(1-p(1))./(1+p(1));
   b=-2.*p(2)./(1+p(1));
   
   %plot([0 rmax],[0 rmax].*m+b,'k-');
   hold off
   
   axis equal
   axis ([rmin rmax rmin rmax]);
   xlabel('r T2'); ylabel('r T1');
   
   title(cellid);
  
   figure(1);
   subplot(2,3,2);
   
   if 1,
      ebmax=max([abs(diff(rbinned(:,1:2),[],2))+...
                 abs(mean(rbinnederr(:,1:2),2));
                 abs(diff(rbinned(:,3:4),[],2))+...
                 abs(mean(rbinnederr(:,3:4),2))]).*1.1;
  
      errorbar(mean(rbinned(:,1:2),2),-diff(rbinned(:,1:2),[],2),...
               mean(rbinnederr(:,1:2),2));
      hold on
      plot([min(rbinned(:)) max(rbinned(:))],[0 0],'k--');
            
      %plot((res.resp(1,2)+res.resp(1,3))./2.*60,...
      %     (res.resp(1,2)-res.resp(1,3)).*60,'ro');
      %plot((res.resp(2,2)+res.resp(2,3))./2.*60,...
      %     (res.resp(2,2)-res.resp(2,3)).*60,'ro');
      
      hold off
      xlabel('average rA+rB');
      ylabel('diff r T1 -r T2');
      aa=axis;
      axis([0 max(rbinned(:)).*1.1 -ebmax ebmax])
      axis square
   else
      scatter(rBb,rAa,'.');
      hold on
      plot([rmin rmax],[rmin rmax],'k--');
      hold off
      axis equal
      axis ([rmin rmax rmin rmax]);
      xlabel('r T1'); ylabel('rBb');
   end
   
   subplot(2,3,3);
   rall=(r(ptidx,1)+r(ptidx,2))./2 .* 60;
   rfdiff=(r(ptidx,1)-r(ptidx,2))./2 .* 60;
   gidx=find(~isnan(rfdiff+rall));
      
   if 1,
      scatter((rAa+rBb)./2,rAa,'r.');
      hold on
      scatter((rAa+rBb)./2,rBb,'k.');
      plot([rmin rmax],[rmin rmax],'k--');
      hold off
      axis ([rmin rmax rmin rmax]);
      
   else
      
      scatter(rall,rfdiff,'.');
      hold on
      scatter(mean(res.resp(1,2:3).*60),...
              (res.resp(1,2)-res.resp(1,3)).*30,40,'r','filled');
      scatter(mean(res.resp(2,2:3).*60),...
              (res.resp(2,2)-res.resp(2,3)).*30,40,'k','filled');
      
      plot([rmin max(rall)],[0 0],'k--');
      plot([0 max(rAB)],[0 max(rAB)].*p(1)+p(2),'k-');
   
      hold off
      %axis equal
      %axis ([rmin rmax -rmax./2 rmax./2]);
      xlabel(sprintf('rAa+Bb s=%.2f',std(rall(gidx))));
      ylabel(sprintf('rAa-Bb s=%.2f',std(rfdiff(gidx))));
   end
   
   axis square
   
   subplot(2,3,4);
   
   scatter(rab,rAB,'.');
   hold on
   if 0
      plot((res.resp(1,2)+res.resperr(1,2).*[-1 1]).*60,...
           res.resp(1,3).*[1 1].*60,'r-');
      plot(res.resp(1,2)*[1 1].*60,...
           (res.resp(1,3)+res.resperr(1,3).*[-1 1]).*60,'r-');
      plot((res.resp(2,2)+res.resperr(2,2).*[-1 1]).*60,...
           res.resp(2,3).*[1 1].*60,'k-');
      plot(res.resp(2,2)*[1 1].*60,...
           (res.resp(2,3)+res.resperr(2,3).*[-1 1]).*60,'k-');
   end
   plot([rmin rmax],[rmin rmax],'k--');
   
   m=(1-p(1))./(1+p(1));
   b=-2.*p(2)./(1+p(1));
   
   %plot([0 rmax],[0 rmax].*m+b,'k-');
   hold off
   
   axis equal
   axis ([rmin rmax rmin rmax]);
   xlabel('r out'); ylabel('r in');
   
   
   subplot(2,3,5);
   errorbar(mean(rbinned(:,3:4),2),-diff(rbinned(:,3:4),[],2),...
                 mean(rbinnederr(:,3:4),2));
   hold on
   plot([min(rbinned(:)) max(rbinned(:))],[0 0],'k--');
   hold off
   axis square
   
   xlabel('average response (spatial in+out)');
   ylabel('spatial in - out')
   
   aa=axis;
   axis([0 max(rbinned(:)).*1.1 -ebmax ebmax])
   
   subplot(2,3,6);
   
   scatter(diffstrfsim(ptidx,1).*60,diffstrfsim(ptidx,2).*30,'.');
   hold on
   scatter(difftargpreds(1,1).*60,difftargpreds(1,2).*30,40,'r','filled')
   scatter(difftargpreds(2,1).*60,difftargpreds(2,2).*30,40,'k','filled')
   plot([rmin max(diffstrfsim(ptidx,1))*60],[0 0],'k--');
   hold off
   %axis equal
   %axis ([rmin rmax -rmax./2 rmax./2]);
   xlabel(sprintf('pAa+Bb s=%.2f',std(pall(gidx))));
   ylabel(sprintf('pAa-Bb s=%.2f',std(pfdiff(gidx))));
   axis square
   
   fullpage('landscape');
   
   % compare A to B
   r1=mean(r(:,[1 3]),2).*13;
   r2=mean(r(:,[2 4]),2).*13;
   useidx=find(~isnan(r1+r2));
   r1=r1(useidx);
   r2=r2(useidx);
   
   if 0
      figure(2)
      clf
      plot(r1,r2,'.');
      hold on
      plot([0 max([r1(:);r2(:)])],[0 max([r1(:);r2(:)])],'k--');
      hold off
      axis equal; axis square;
      
      % mean diff test
      mean([r1 r2])
      
      x=r1;y=r2;
      
      keyboard
   end
   
end

figure(2);
clf

%ap=fpatches;
tp=mpatches-repmat(mS(:,1),[1 size(mpatches,2)]);
%ap = ap ./ repmat(mean(ap,1),size(ap,1),1);
%tp = tp ./ repmat(mean(abs(tp),1),size(tp,1),1);
tp = tp ./ mean(abs(tp(:)));
%mp=repmat(mean(ap,2),1,size(tp,2));
%tp=tp-mp;

targsim=(tp(:,2:end)'* (bstim'-repmat(mS(:,1),1,blen)))';

kcount=size(H2,3);

if size(tp,1)>spacecount,
   targmtx=rand(spacecount,1,kcount);
elseif size(tp,2)<kcount,
   targmtx=reshape(tp,spacecount,1,size(tp,2));
   targmtx(:,1,size(tp,2)+1:kcount)=0;
   targmtx(1,1,size(tp,2)+1:kcount)=1;
   targmtx=abs(targmtx).^0.5 .* sign(targmtx);
else
   targmtx=reshape(tp(:,1:kcount),spacecount,1,kcount);
   targmtx=abs(targmtx).^0.5 .* sign(targmtx);
end
targmtx=targmtx./ max(abs(targmtx(:))) .* max(abs(H2(:)));
%for attidx=1:kcount,
%   targmtx(:,:,attidx)=targmtx(:,:,attidx) ./ ...
%       max(abs(targmtx(:,:,attidx))) .* max(max(abs(H2(:,:,attidx))));
%end
eigmtx=cat(2,targmtx,targmtx,H2,...
           cat(3,H(:,1),H2(:,2)-H2(:,3),H2(:,3)-H2(:,2)),...
           H3,cat(3,H(:,1),H3(:,2)-H3(:,3),H3(:,3)-H3(:,2)));
%           H4,cat(3,H(:,1),H4(:,2)-H4(:,3),H4(:,3)-H4(:,2)),...
eigmtx(:,:,size(eigmtx,3)+1)=nan;

colcount=size(eigmtx,2);
rowcount=size(eigmtx,3);

% depending on what flag is set, either show in sf/or map or native kernfmt
if showSFOR
   obins=linspace(0,180,obincount+1);
   obins=obins(1:end-1);
   sfbins=linspace(1,8,8);
   sformtx=sf2gr(eigmtx,length(obins),length(sfbins),0,0,params.kernfmt);
   sformtx=permute(sformtx,[2 1 3 4]);
elseif strcmp(params.kernfmt,'space'),
   Xmax=sqrt(size(eigmtx,1));
   sformtx=reshape(eigmtx,Xmax,Xmax,colcount,rowcount);
   sformtx=permute(sformtx,[2 1 3 4]);
   sformtx=flipdim(sformtx,1);
   obins=1:Xmax;
   sfbins=1:Xmax;
else
   sformtx=pfft2sf(eigmtx,params.kernfmt);
   Xmax=size(sformtx,1);
   xc=round((Xmax+1)/2);
   obins=(1:Xmax)-xc;
   sfbins=obins;
end

% find contour levels
ta=sformtx(:,:,3,1);
astd=std(ta(:));

% smoothing filters for contours
xx=-2:2;
gsf=exp(-(xx./0.8).^2/2);
gsf=gsf'./sum(gsf(:));
gor=exp(-(xx./1.0).^2/2);
gor=(gor./sum(gor(:)));

for rr=1:rowcount-1,
   if rr==1,
      crange=3;
   else
      crange=2:colcount;
   end
   for cc=crange,
      subplot(rowcount,colcount,(rr-1)*colcount+cc);
      ta=sformtx(:,:,cc,rr);
      
      smta=rconv2(ta,gsf); % ,'same');
      smta=cconv2(smta,gor);     %orientation is circular
      
      if cc<=3,
         mm=max(max(max(abs(sformtx(:,:,cc,rr)))));
      elseif cc>=3,
         mm=max(max(max(max(abs(sformtx(:,:,3:end,:))))));
      else
         mm=max(max(abs(sformtx(:,:,1,1))));
      end
      if mm==0,
         mm=1;
      end
      
      ta=imresize(ta,4,'bilinear');
      xx=-8:8;
      gsf2=exp(-(xx./2.4).^2/2);
      gsf2=gsf2'./sum(gsf2(:));
      gor2=exp(-(xx./4).^2/2);
      gor2=(gor2./sum(gor2(:)));
      ta=rconv2(ta,gsf2); % ,'same');
      ta=cconv2(ta,gor2);     %orientation is circular
      
      ta=repmat(ta,[1 1 3])./mm;
      %ta=sqrt(abs(ta)).*sign(ta);
      ta(:,:,1)=1+ta(:,:,1).*(ta(:,:,1)<0);
      ta(:,:,2)=1-abs(ta(:,:,2));
      ta(:,:,3)=1-ta(:,:,3).*(ta(:,:,3)>0);
      
      ta(find(ta>1))=1;
      ta(find(ta<0))=0;
      ta=ta.^1.1;
      
      if 1,
         imagesc(obins,sfbins,ta);
         
         hold on
         slev1=1.0;
         slev2=2.0;
         contour(obins,sfbins,smta,[-astd*slev2  -astd*slev2],'k-');
         contour(obins,sfbins,smta,[-astd*slev1  -astd*slev1],'k--');
         contour(obins,sfbins,smta,[ astd*slev1   astd*slev1],'k--');
         contour(obins,sfbins,smta,[ astd*slev2   astd*slev2],'k-');
         hold off
      else
         imagesc(obins,sfbins,smta);
         
      end
      
      if showSFOR,
         axis xy;
      else
         axis image;
      end
      
      if rr<rowcount-1,
         set(gca,'XTickLabel',[]);
      elseif showSFOR,
         xlabel('orientation');
      else
         xlabel('freq x');
      end
      if cc>1,
         set(gca,'YTickLabel',[]);
      elseif showSFOR,
         ylabel('spatial freq');
      else
         ylabel('freq y');
      end   
   end
end

z=load(params.respfiles{1});
targid1=unique(z.picid(find(z.targetmatch==1)));
targid2=unique(z.picid(find(z.targetmatch==2)));
otargpatches=cat(3,dmsorigframe(targid1),dmsorigframe(targid2));
bigpatches=otargpatches;

[xx,yy]=meshgrid(-64:63,-64:63);
dd=sqrt(xx.^2+yy.^2);
rr=(dd-61)./3;
rr(rr>1)=1;
rr(rr<0)=0;
rr=1-rr;
bigpatches(:,:,1)=((bigpatches(:,:,1)-46)./(rr+(rr==0)) + 46 ) .*rr;
bigpatches(:,:,2)=((bigpatches(:,:,2)-46)./(rr+(rr==0)) + 46 ) .*rr;


gam=2;
bigpatches=(bigpatches/255).^(1./gam);
mmax=max(max(max(bigpatches)));
for rr=2:rowcount-1,
   subplot(rowcount,colcount,rr.*colcount-colcount+1);
   
   ta=repmat(bigpatches(:,:,rr-1),[1 1 3])./mmax;
   imagesc(ta);
   axis off
   axis image
end


maxr=max(bresp(:,1));
xx=linspace(-maxr,maxr,100);

subplot(rowcount,colcount,1);
plot(targpreds(:,1),'k--');
hold on
plot(targpreds(:,2:end));
hold off
a=axis;
axis([0 kcount+2 a(3) a(4)]);
sleg={'all'};
for aa=2:kcount,
   sleg{aa}=sprintf('%s',num2str('A'+aa-2));
end
legend(sleg);


if 1, 
   pall=rprec(ptidx,1).*60;
   rA=r(ptidx,1).*60;
   rB=r(ptidx,2).*60;
   rAa=(r(ptidx,1)+r(ptidx,3))./2.*60;
   rBb=(r(ptidx,2)+r(ptidx,4))./2.*60;
   rAB=(r(ptidx,1)+r(ptidx,2))./2.*60;
   rab=(r(ptidx,3)+r(ptidx,4))./2.*60;
   NN=7;
   spall=sort(pall);
   spall=[spall(:); spall(end)+1];
   epall=spall(round(linspace(1,length(spall),NN+1)));
   pbinned=zeros(NN,1);
   rbinned=zeros(NN,4);
   rbinnederr=zeros(NN,4);
   for nn=1:NN,
      ff=find(pall>=epall(nn) & pall<epall(nn+1));
      pbinned(nn)=nanmean(pall(ff));
      rbinned(nn,1)=nanmean(rAa(ff));
      rbinnederr(nn,1)=nanstd(rAa(ff))./sqrt(length(ff));
      rbinned(nn,2)=nanmean(rBb(ff));
      rbinnederr(nn,2)=nanstd(rBb(ff))./sqrt(length(ff));
      rbinned(nn,3)=nanmean(rAB(ff));
      rbinnederr(nn,3)=nanstd(rAB(ff))./sqrt(length(ff));
      rbinned(nn,4)=nanmean(rab(ff));
      rbinnederr(nn,4)=nanstd(rab(ff))./sqrt(length(ff));
   end
   arange=[pbinned(1).*0.8 pbinned(end).*1.2 ...
           min(rbinned(:)-rbinnederr(:)).*0.9 ...
           max(rbinned(:)+rbinnederr(:)).*1.1];
   
   subplot(rowcount,colcount,colcount-1);
   
   plot((rAa+rBb)./2,rBb,'.','Color',[0.6 0.6 0.6]);
   hold on
   plot((rAa+rBb)./2,rAa,'.','Color',[1 0.5 0.5]);
   ht=errorbar(pbinned,rbinned(:,2),rbinnederr(:,2),'k');
   ht=[ht errorbar(pbinned,rbinned(:,1),rbinnederr(:,1),'r')];
   
   errorbar(res.strfresp(1,1).*60,res.resp(1,2).*60,...
            res.resperr(1,2).*60,'rx');
   errorbar(res.strfresp(1,1).*60,res.resp(1,3).*60,...
            res.resperr(1,3).*60,'kx');
   errorbar(res.strfresp(2,1).*60,res.resp(2,3).*60,...
            res.resperr(2,3).*60,'ro');
   errorbar(res.strfresp(2,1).*60,res.resp(2,2).*60,...
            res.resperr(2,2).*60,'ko');
   
   
   hold off
   axis(arange);
   
   subplot(rowcount,colcount,colcount);
   plot((rAa+rBb)./2,rab,'.','Color',[0.5 0.5 1]);
   hold on
   plot((rAa+rBb)./2,rAB,'.','Color',[0.5 0.8 0.5]);
   ht=[ht errorbar(pbinned,rbinned(:,4),rbinnederr(:,4),'b')];
   ht=[ht errorbar(pbinned,rbinned(:,3),rbinnederr(:,3),'g')];
   hold off
   axis(arange);
   
   set(ht,'LineWidth',2);
   
else
   
   subplot(rowcount,colcount,colcount-1);
   
   mdc=mean(globaldc(2:end,:),2);
   edc=std(globaldc(2:end,:),1,2) .* sqrt(params.resampcount-1);
   plot(linspace(0,kcount,params.bootcount),globaldc(2:end,:)','--');
   hold on
   errorbar(mdc,edc.*2,'ko');
   hold off
   title('global dc');
   
   subplot(rowcount,colcount,colcount);
   mdc=mean(globalgain(2:end,:),2);
   edc=std(globalgain(2:end,:),1,2) .* sqrt(params.resampcount-1);
   plot(linspace(0,kcount,params.bootcount),globalgain(2:end,:)','--');
   hold on
   errorbar(mdc,edc.*2,'ko');
   hold off
   title('global gain');
end

subplot(rowcount,colcount,2);
cla
axis off;
subplot(rowcount,colcount,3);
title(sprintf('%s no-att kernel',cellid));
for attidx=2:kcount,
   subplot(rowcount,colcount,attidx*colcount-colcount+1);
   title(sprintf('targ %s',num2str('A'-2+attidx)));
end

subplot(rowcount,colcount,colcount+3);
title('Aa');
subplot(rowcount,colcount,colcount+4);
title('Aa-Bb');
subplot(rowcount,colcount,colcount+5);
title('AB');
subplot(rowcount,colcount,colcount+6);
title('AB-ab');
subplot(rowcount,colcount,2*colcount+3);
title('Bb');
subplot(rowcount,colcount,2*colcount+4);
title('Bb-Aa');
subplot(rowcount,colcount,2*colcount+5);
title('ab');
subplot(rowcount,colcount,2*colcount+6);
title('ab-AB');


subplot(rowcount,1,rowcount);
axis off;
axis ij;

details={};
details{1}=sprintf('%s: pred xc :',cellid);
details{2}=[sprintf('attidx     :         '),...
            sprintf(' %11d',1:kcount)];
details{3}=[sprintf('framecount :         '),...
            sprintf(' %11d',ncount)];

for nlidx=1:attcount,
   details{3+nlidx}=[...
      sprintf('%-11s:',nlnames{nlidx}),...
      sprintf('%5.2f/%5.2f',predxc(1:2,nlidx)),...
      sprintf(' (p<%4.2f):',pxc(nlidx))];
end
ht=text(0,0,details,'VerticalAlignment','top','FontName','Courier');

fullpage landscape

%disp('skipping stim disp');
%return

if length(find(picid<0 & picid>=-12))>0,
   % ie, there are gratings in there
   
   rrtemp=load(params.respfiles{1});
   
   %ttidx=shift(-1:-1:-12,[1 6]);
   ttidx=fliplr(unique(picid(find(picid<0)))');
   
   gcount=zeros(length(ttidx),size(r2,2));
   gresp=zeros(length(ttidx),size(r2,2));
   gresperr=zeros(length(ttidx),size(r2,2));
   gframe=zeros(length(ttidx),1);
   
   for ii=1:length(ttidx),
      
      gframe(ii)=min([length(rrtemp.picid)+1 find(rrtemp.picid==ttidx(ii))]);
      
      for jj=1:size(gresp,2),
         gcount(ii,jj)=length(find(picid==ttidx(ii) & ~isnan(r2(:,jj))));
         gresp(ii,jj)=nanmean(r2(find(picid==ttidx(ii)),jj));
         gresperr(ii,jj)=nanstd(r2(find(picid==ttidx(ii)),jj));
      end
   end
   
   %keyboard
   ggmov=loadimframes(params.stimfiles{1},(gframe),params.stimloadparms{:});
   fgmov=feval(params.stimfiltercmd,ggmov,params.stimfilterparms{:});
   
   glen=size(fgmov,2);
   
   % how similar each grating is to the strf
   gratpreds=(H2(:,:)'* ...
              (fgmov-repmat(mS(:,1),1,glen)))' + ...
             repmat(nlparms2(1,:),[glen 1]);
   difftargpreds=(H(:,:)'* ...
       (mpatches(:,2:end)-repmat(mS(:,1),1,size(mpatches,2)-1)))' + ...
       repmat(nlparms(1,:),[size(mpatches,2)-1 1]);
   
   ggmov=reshape(ggmov,size(ggmov,1)*size(ggmov,2),size(ggmov,3));
   figure(4);
   tm=cat(3,zeros(size(ggmov)),ggmov);
   showkern(tm,'space');
   subplot(2,1,1);
   plot([mean(gresp(:,[1 ]),2) mean(gresp(:,[2 ]),2)].*60);
   hold on
   %plot(gratpreds(:,2:3).*60,'--');
   hold off
   legend('A','B')
   %legend('Aa','Bb','predA','predB')
end

drawnow

picid=[];
for fidx=1:filecount,
   [tr,tpicid]=resploaddmsmatch(params.respfiles{fidx});
   picid=cat(1,picid,tpicid);
end

bresp(bresp(:,1)<0,1)=0;

rresp=bresp(:,1);
rrespA=(r(:,1)+r(:,3));%-(r(:,2)+r(:,4));
rrespB=(r(:,2)+r(:,4));%-(r(:,1)+r(:,3));
%rrespA(r(:,1)<r(:,2))=0;
%rrespB(r(:,1)>r(:,2))=0;
rrespA(r(:,1)<r(:,2) | r(:,3)<r(:,4))=0;
rrespB(r(:,1)>r(:,2) | r(:,3)>r(:,4))=0;

%otargpatches=otargpatches(15:114,15:114,:);
%keyboard

%rrespA=rprec(:,1)+rprec(:,2);
%rrespB=rprec(:,1)-rprec(:,2);
%rrespA=rprec(:,2);
%rrespB=-rprec(:,2);
mresp=rprec(:,1);
avgresp=(rprec(:,1)+bresp(:,1))./2;

rresp(find(isnan(bresp(:,4))))=nan;
mresp(find(isnan(bresp(:,2))))=nan;
avgresp(find(isnan(bresp(:,2))))=nan;

[xx,bothidx]=sortrows(-avgresp);
[xx,strongidx]=sortrows(-mresp);
[xx,weakidx]=sortrows(mresp);
[xx,goodidx]=sortrows(-rresp);
[xx,badidx]=sortrows(rresp);
[xx,goodidxA]=sortrows(-rrespA);
[xx,goodidxB]=sortrows(-rrespB);

scount=20;

if 1,
   % used this for V4 attention paper
   fstimset=zeros(size(bstim,2),scount,5).*nan;;
   stimset=zeros(128,128,scount,5).*nan;
   
   mid=round(sum(~isnan(rresp))/2);
   
   %idset=[strongidx(10:-1:1) goodidx(10:-1:1) ...
   %       goodidx((10:-2:-8)+mid) badidx(10:-1:1)];
   idset=[goodidxA(1:scount) goodidxA((1:scount)+scount)  ...
          goodidxB(1:scount) goodidxB((1:scount)+scount)];
   
   disp('loading original stim patterns');
   
   [xx,yy]=meshgrid(-64:63,-64:63);
   dd=sqrt(xx.^2+yy.^2);
   rr=(dd-61)./3;
   rr(rr>1)=1;
   rr(rr<0)=0;
   rr=1-rr;
   %keyboard
   for ii=1:length(idset(:)),
      tts=dmsorigframe(picid(idset(ii)));
      
      % adjust ring to go to black rather than gray
      tts=((tts-46)./(rr+(rr==0)) + 46 ) .*rr;
      
      stimset(:,:,ii)=tts;
      fstimset(:,ii)=bstim(idset(ii),:)';
   end
   
   catstr={'gactA1','gactA2','gactB1','gactB2'};

elseif 1,
   % used this for V4 tuning paper
   
   fstimset=zeros(size(bstim,2),scount,7).*nan;;
   %stimset=zeros(100,100,scount,7).*nan;
   stimset=zeros(128,128,scount,7).*nan;
   
   mid=round(sum(~isnan(rresp))/2);
   
   %idset=[strongidx(10:-1:1) goodidx(10:-1:1) ...
   %       goodidx((10:-2:-8)+mid) badidx(10:-1:1)];
   idset=[strongidx(scount:-1:1) goodidx(scount:-1:1) ...
          goodidx(((scount./2):-1:(-scount./2+1))+mid) badidx(scount:-1:1) ...
          goodidxA(1:scount) goodidxB(1:scount)];
   
   disp('loading original stim patterns');
   for ii=1:length(idset(:)),
      tts=dmsorigframe(picid(idset(ii)));
      %tts=zeros(128);
      %stimset(:,:,ii)=tts(15:114,15:114);
      stimset(:,:,ii)=tts;
      
      fstimset(:,ii)=bstim(idset(ii),:)';
   end
   
   catstr={'gpred','gact','mact','bact','gactA','gactB'};
else
   
   if length(params.stimfilterparms)>=5,
      bgpix=params.stimfilterparms{5};
   else
      bgpix=46;
   end
   showstim=stim(:,:,goodidx(1:10))-bgpix;
   
   showstim=showstim.*repmat(hanning2(size(showstim,1)),...
                             [1 1 size(showstim,3)]);
   
   hmask=fftshift(pfft2sf(H2(:,1),params.kernfmt));
   pmask=hmask.*(hmask>0);
   nmask=hmask.*(hmask<0);
   fstimset=fft(fft(showstim,[],1),[],2);
   
   pstimset=fstimset.*repmat(pmask,[1 1 size(fstimset,3)]);
   nstimset=fstimset.*repmat(nmask,[1 1 size(fstimset,3)]);
   
   pstimset=real(ifft(ifft(pstimset,[],1),[],2));
   nstimset=real(ifft(ifft(nstimset,[],1),[],2));
   
   stimset(:,:,:,1)=stim(:,:,strongidx(1:10));
   stimset(:,:,:,2)=stim(:,:,goodidx(1:10));
   stimset(:,:,:,3)=pstimset;
   stimset(:,:,:,4)=nstimset;
   
   idset=[strongidx(1:10) goodidx(1:10) goodidx(1:10) goodidx(1:10)];
   catstr={'gpred','gact','gp','gn'};
end

srows=size(stimset,4);
stimset=reshape(stimset,size(stimset,1)*size(stimset,2),...
                 scount,srows);

if strcmp(showextra,'curves'),
   % curves generated above to be sent out to res.
   
   ormin=min(ortuning(:))-max(abs(ortuning(:))*0.1);
   ormax=max(ortuning(:))+max(abs(ortuning(:))*0.1);
   sfmin=min(sftuning(:))-max(abs(sftuning(:))*0.1);
   sfmax=max(sftuning(:))+max(abs(sftuning(:))*0.1);
   
   figure(3);
   clf
   
   showkern(cat(3,Hset,Hset(:,:,1:2)),'pfftgr',[],{},1);
   
   subplot(5,2,1);
   title(sprintf('%s baseline',cellid));
   subplot(5,2,2);
   cla; axis off
   mdc=mean(globaldc(2:end,:),2);
   edc=std(globaldc(2:end,:),1,2) .* sqrt(params.resampcount-1);
   plot(linspace(0,kcount,params.bootcount),globaldc(2:end,:)','--');
   hold on
   errorbar(mdc,edc.*2,'ko');
   hold off
   title('global dc');
   
   subplot(5,2,3);
   title(sprintf('Aa ft: %5.2f/%5.2f (p<%4.2f)',...
                 predxc(1:2,2),pxc(2)));
   subplot(5,2,4);
   title(sprintf('AB sp: %5.2f/%5.2f (p<%4.2f)',...
                 predxc(1:2,3),pxc(3)));
   subplot(5,2,5);
   title(sprintf('Bb'));
   subplot(5,2,6);
   title(sprintf('ab'));
   
   subplot(5,2,7);
   plot(obins,squeeze(ortuning(:,1,2:end)),'LineWidth',2);
   a=axis;
   axis([min(obins) max(obins) ormin ormax]);
   title(sprintf('Aa/Bb or0: %.1f/%.1f orbw: %.1f/%.1f',...
                 orpeak(1,2:end),orbw(1,2:end)));
   
   subplot(5,2,9);
   plot(sfbins,squeeze(sftuning(:,1,2:end)),'LineWidth',2);
   a=axis;
   axis([min(sfbins) max(sfbins) sfmin sfmax]);
   title(sprintf('Aa/Bb sf0: %.1f/%.1f sfbw: %.1f/%.1f',...
                 sfpeak(1,2:end),sfbw(1,2:end)));
   xlabel('spatial freq');
   legend('Aa','Bb');
   
   if 1,
      subplot(5,2,8);
      mstimor=squeeze(mean(sformtx(:,:,2,2:3),1));
      mstimsf=squeeze(mean(sformtx(:,:,2,2:3),2));
      
      plot(obins,mstimor,'LineWidth',2);
      %a=axis;
      %axis([min(obins) max(obins) ormin ormax]);
      title(sprintf('target or0: %.1f/%.1f orbw: %.1f/%.1f',...
                    orpeak(2,2:end),orbw(2,2:end)));
      
      subplot(5,2,10);
      plot(sfbins,mstimsf,'LineWidth',2);
      %a=axis;
      %axis([min(sfbins) max(sfbins) sfmin sfmax]);
      title(sprintf('targ sf0: %.1f/%.1f sfbw: %.1f/%.1f',...
                    sfpeak(2,2:end),sfbw(2,2:end)));
      xlabel('spatial freq');
      legend('AB','ab');
   
   else
      subplot(5,2,8);
      plot(obins,squeeze(ortuning(:,2,2:end)),'LineWidth',2);
      a=axis;
      axis([min(obins) max(obins) ormin ormax]);
      title(sprintf('AB/ab or0: %.1f/%.1f orbw: %.1f/%.1f',...
                    orpeak(2,2:end),orbw(2,2:end)));
      
      subplot(5,2,10);
      plot(sfbins,squeeze(sftuning(:,2,2:end)),'LineWidth',2);
      a=axis;
      axis([min(sfbins) max(sfbins) sfmin sfmax]);
      title(sprintf('AB/ab sf0: %.1f/%.1f sfbw: %.1f/%.1f',...
                    sfpeak(2,2:end),sfbw(2,2:end)));
      xlabel('spatial freq');
      legend('AB','ab');
   end
   
   fullpage portrait
   
elseif strcmp(showextra,'stim'),
   figure(3);
   % display with gamma-esque correction 
   gam=2;
   stimset=(stimset/255).^(1./gam) .*255-127;
   showkern(stimset,'space',0,{},0,1,scount);
   for ii=1:size(idset,1),
      for jj=1:size(idset,2);
         subplot(size(stimset,3),size(idset,1),ii+(jj-1)*size(idset,1));
         if strcmp(catstr{jj}(2:end),'act'),
            title(sprintf('%s %.1f',catstr{jj},bresp(idset(ii,jj)).*60));
         else
            title(sprintf('%s %.1f',catstr{jj},rprec(idset(ii,jj),1).*60));
         end
         axis off
      end
   end

   subplot(srows,2,srows*2-1);
   plot(sort(rresp(~isnan(rresp))).*60,'k-','LineWidth',2);
   title([cellid,' act response']);
   a=axis;
   axis([-10 sum(~isnan(rresp))+10 0 a(4)]);
   
   subplot(srows,2,srows*2);
   plot(sort(mresp(~isnan(rresp))).*60,'k-','LineWidth',2);
   title('pred response');
   a=axis;
   axis([-10 sum(~isnan(rresp))+10 0 a(4)]);
   
   fullpage landscape
   colormap(gray);
   
   figure(4);
   clf
   
   disp('showing pref/non-pref stimuli in fourier domain');
   if strcmp(params.kernfmt,'space'),
      showkern(fstimset(:,:,1:(srows-1))-repmat(mS(:,1),[1 scount srows-1]),...
               params.kernfmt,0,{},0,1,scount);
   else
      showkern(fstimset(:,:,1:(srows-1))-repmat(mS(:,1),[1 scount srows-1]),...
               'pfftgr',0,{},0,1,scount);
   end
   
   fullpage landscape
   
elseif strcmp(showextra,'pred'),
   figure(3);
   clf
   
   dist('showing pos/neg componeents of predictions');
   gidx=find(~isnan(bresp(:,1))& bresp(:,1)>=0) ;
   scatter(pstrfsim(gidx),nstrfsim(gidx),round(bresp(gidx,1).*250+1));
   hold on
   scatter(pstrfsim(goodidx(1:10)),nstrfsim(goodidx(1:10)),...
           round(bresp(goodidx(1:10),1).*300+1),'r','filled');
   hold off
   axis equal
end



fprintf('to print SRFs:   print -f2 -depsc output/%s.attsrf.eps\n',cellid);



return



