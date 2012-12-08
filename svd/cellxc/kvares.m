% function res=kvares(cellid,batch,showextra)
%
% load results from kernfile in sRunData and display attentional
% modulation info ... kernfile generated from kerncomp4
%
% created SVD 10/18/02 - hacked from kerncomp4res.m
%
function res=kvares(runidx,batch,showextra)

if ~exist('showextra','var'),
   showextra='';
end

dbopen;
if isnumeric(runidx),
   sql=['SELECT * from sRunData WHERE id=',num2str(runidx)];
   rundata=mysql(sql);
   batchcount=length(rundata);
   cellid=rundata.cellid;
else
   cellid=runidx;
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

if length(rundata)>1,
   disp('cellres-kerncomp: only displaying first batch');
end

fprintf('kvares.m: loading %s\n',...
        [rundata(1).respath,rundata(1).kernfile,'.gz']);
zload([rundata(1).respath,rundata(1).kernfile,'.gz']);

% show spatial kernels for each attention condition and jackknife

kernfmt=params.kernfmt;
if strcmp(kernfmt,'pfftgr'),
   kernfmt='pfft';
   iconside=[sqrt(spacecount*2) sqrt(spacecount*2)];
end

H=zeros([size(vstrf(1).h,1),params.bootcount,attcount]);
mS=zeros([size(vstrf(1).mS,1),params.bootcount,attcount]);
H(:,:,1)=cat(2,vstrf(1,1,1,:).h);
mS(:,:,1)=cat(2,vstrf(1,1,1,:).mS);

if isfield(params,'oddout') & params.oddout,
   for attidx=1:attcount-1,
      H(:,:,attidx+1)=cat(2,vstrf(4+attidx,attidx,1,:).h);
      mS(:,:,attidx+1)=cat(2,vstrf(4+attidx,attidx,1,:).mS);
   end
else
   for attidx=2:attcount,
      H(:,:,attidx)=cat(2,vstrf(end,attidx-1,1,:).h);
      mS(:,:,attidx)=cat(2,vstrf(end,attidx-1,1,:).mS);
   end
end

if nargout==0,
   
   figure(1);
   showkern(H,kernfmt,iconside,{},0);
   
   kshowcount=min([params.bootcount 12]);
   for ii=1:kshowcount,
      subplot(attcount,kshowcount,ii),
      if ii==1,
         title(sprintf('%s %d/%d',cellid,vstrf(1,1,1,ii).parms.sfsfit,...
                       vstrf(1,1,1,ii).parms.sigfit));
      else
         title(sprintf('%d/%d',vstrf(1,1,1,ii).parms.sfsfit,...
                       vstrf(1,1,1,ii).parms.sigfit));
      end
      
      for attidx=2:attcount,
         subplot(attcount,kshowcount,(attidx-1).*kshowcount+ii),
         title(sprintf('%d/%d',vstrf(end,attidx-1,1,ii).parms.sfsfit,...
                    vstrf(end,attidx-1,1,ii).parms.sigfit));
      end
   end
else
   disp('nargout>0: skipping figures');
end


H2=H;
dcgp=ones(2,attcount);
dcgpall=ones(2,params.bootcount,attcount);
for bootidx=1:params.bootcount,
   dcgpall(:,bootidx,1)=vstrf(1,1,1,bootidx).nlparms';
end
dcgp(:,1)=vstrf(1,1,1,1).nlparms;
for attidx=2:attcount,
   for bootidx=1:params.bootcount,
      % add local onto baseline kernel, scaled according to gain.
      %H2(:,1,attidx)=H2(:,1,attidx).*mean(globalgain(attidx-1,:))+...
      %    H(:,:,attidx);
      H2(:,bootidx,attidx)=H(:,bootidx,1).*globalgain(attidx-1,bootidx)+...
          H(:,bootidx,attidx);
      dcgpall(1,bootidx,attidx)=globaldc(attidx-1,bootidx);
   end
   
   
   
   dcgp(1,attidx)=mean(globaldc(attidx-1,:));
end

if 0,

   spatsim=zeros(size(H,3),size(H,3),size(H,2));
   pspatsim=zeros(size(H,3).*2,size(H,3).*2,size(H,2));
   nspatsim=zeros(size(H,3),size(H,3),size(H,2));
   targsim=zeros(size(H,3),size(H,3),size(H,2));
   
   tp=mpatches./repmat(mS(:,1),[1 size(mpatches,2)]);
   tp=mpatches;
   targsim=tp'*tp;
   targsim=targsim ./ sqrt(diag(targsim)*diag(targsim)');
   for bootidx=1:size(H,2),
      tsim=squeeze(H(:,bootidx,:))' * squeeze(H(:,bootidx,:));
      tsim=tsim+diag(diag(tsim)==0);
      spatsim(:,:,bootidx)=tsim ./ sqrt(diag(tsim)*diag(tsim)');
   
      th=[squeeze(H2(:,bootidx,:).*(H2(:,bootidx,:)>0)) ...
          squeeze(H2(:,bootidx,:).*(H2(:,bootidx,:)<0))];
      tsim=th'*th;
      tsim=tsim+diag(diag(tsim)==0);
      pspatsim(:,:,bootidx)=tsim ./ sqrt(diag(tsim)*diag(tsim)');
      
   end
   
   mspatsim=mean(pspatsim,3)
   espatsim=std(pspatsim,0,3).*sqrt(params.bootcount);
   
   figure(4);
   clf
   errorbar(mspatsim(1,2:attcount),espatsim(1,2:attcount),'r-');
   hold on
   errorbar(mspatsim(attcount+1,(2:attcount)+attcount),...
            espatsim(attcount+1,(2:attcount)+attcount),'b-');
   hold off
   keyboard
end

pcount=size(mpatches,2)-1;
strfsimboot=zeros(blen,params.bootcount,attcount);
targpredboot=zeros(pcount,params.bootcount,attcount);
for attidx=1:attcount
   strfsimboot(:,:,attidx)=...
       (H2(:,:,attidx)'* (bstim'-repmat(mS(:,1),1,blen)))' + ...
       repmat(dcgpall(1,:,attidx),[blen 1]);
   targpredboot(:,:,attidx)=...
       (mpatches(:,2:end)-repmat(mS(:,1),1,pcount))'*H2(:,:,attidx) + ...
       repmat(dcgpall(1,:,attidx),[pcount 1]);
end

eH=std(H,0,2).*sqrt(params.bootcount);
H=mean(H,2);
shH=shrinkage(H,eH,0.5);
mS=mean(mS,2);
H2=mean(H2,2);

% true local att difference, factoring in global gain
Hd=H2-repmat(H2(:,:,1),[1 1 attcount]);

% how similar each target is to the strf
targpreds=(H2(:,:)'* ...
    (mpatches(:,2:end)-repmat(mS(:,1),1,size(mpatches,2)-1)))' + ...
    repmat(dcgp(1,:),[size(mpatches,2)-1 1]);
patchpreds=(H2(:,:)'* ...
    (fpatches-repmat(mS(:,1),1,size(fpatches,2))))' + ...
    repmat(dcgp(1,:),[size(fpatches,2) 1]);
targpredsnodc=(H2(:,:)'* ...
    (mpatches(:,2:end)-repmat(mS(:,1),1,size(mpatches,2)-1)))';
targpredspos=((H2(:,:).*(H2(:,:)>0))'* ...
    (mpatches(:,2:end)-repmat(mS(:,1),1,size(mpatches,2)-1)))';
targpredsneg=((H2(:,:).*(H2(:,:)<0))'* ...
    (mpatches(:,2:end)-repmat(mS(:,1),1,size(mpatches,2)-1)))';


targpredcc=zeros(size(mpatches,2)-1,attcount);
for attidx=1:attcount,
   for targidx=1:size(mpatches,2)-1,
      if var(H2(:,attidx))>0 & var(mpatches(:,targidx+1)-mS(:,1))>0,
         targpredcc(targidx,attidx)=...
             xcorr(H2(:,attidx),mpatches(:,targidx+1)-mS(:,1),0,'coeff');
      else
         %keyboard
      end
   end
end

% how similar each stimulus is to the strf
strfsim=(H2(:,:)'* (bstim'-repmat(mS(:,1),1,blen)))' + ...
        repmat(dcgp(1,:),[blen 1]);
strfdiffsim=(shH(:,:)'* (bstim'-repmat(mS(:,1),1,blen)))' + ...
        repmat(dcgp(1,:),[blen 1]);

if size(mpatches,2)==5 & attcount==5,
   % matched filter test
   predin=zeros(attcount-1,1);
   predout=zeros(attcount-1,1);
   pp=squeeze(median(targpredboot,2));

   for attidx=1:size(pp,1),
      %predin(attidx)=targpreds(attidx,attidx+1);
      %predout(attidx)=mean(targpreds(attidx,[2:attidx attidx+1:end]));
      predin(attidx)=pp(attidx,attidx+1);
      predout(attidx)=mean(pp(attidx,[2:attidx attidx+1:end]));
   end
   predin=predin.*(predin>0);
   predout=predout.*(predout>0);
   
   % pull out targets in image domain
   ttpspace=bigpatches(:,:,targlist);
   ttpspace=ttpspace-mean(ttpspace(:));  %vs bigpatches 
   tpspace=zeros(32,32,length(targlist));
   for attidx=1:length(targlist),
      tpspace(:,:,attidx)=imresize(ttpspace(:,:,attidx),[32 32])-40;
   end
   tpspace=reshape(tpspace,size(tpspace,1)*size(tpspace,2),...
                   size(tpspace,3));
   
   tpm=mS(:,1);
   %tpm=mean(fpatches,2);
   tp=(mpatches(:,2:end)-repmat(tpm,[1 size(mpatches,2)-1]));
   
   tpred=strfsim(:,2:end)-repmat(strfsim(:,1),[1 size(strfsim,2)-1]);
   
   tdcpred=repmat(dcgp(1,2:end)-mean(dcgp(1,2:end)),[size(tpred,1) 1]);
   
   targsim=zeros(paircount,1);
   prefsim=zeros(paircount,1);
   predsim=zeros(paircount,1);
   dcsim=zeros(paircount,1);
   
   for pidx=1:paircount,
      p1=pairidx(pidx,1);
      p2=pairidx(pidx,2);
      
      % simrat 2/3 seems best
      %simrat=1/2;
      %simrat=2/3;
      simrat=0.85;
      %simrat=1; % all pfft domain
      
      % DON'T use xcov here
      targsim(pidx)=simrat*xcorr(tp(:,p1),tp(:,p2),0,'coeff') + ...
          (1-simrat)*xcorr(tpspace(:,p1),tpspace(:,p2),0,'coeff'); %,'coeff'
      
      % sim of normalized pred resp to each target in baseline
      %prefsim(pidx)=1-abs(targpredcc(p1,1)-targpredcc(p2,1))./...
      %    abs(targpredcc(p1,1)+targpredcc(p2,1));
      
      % sim of predicted resp to each target in baseline
      prefset=sort(1-abs(targpredboot(p1,:,1)-targpredboot(p2,:,1))./...
           abs(targpredboot(p1,:,1)+targpredboot(p2,:,1)));
      prefsim(pidx)=median(prefset(2:end-1));
      
      %prefsim(pidx)=1-abs(targpreds(p1,1)-targpreds(p2,1))./...
      %    abs(targpreds(p1,1)+targpreds(p2,1)-2*dcgp(1));
      %prefsim(pidx)=(targpreds(p1,p1+1)-mean(targpreds(:,1))) * ...
      %    (targpreds(p2,p2+1)-mean(targpreds(:,1)))./ ...
      %    var(targpreds(:,1));
      %prefsim(pidx)=(H(:,1))' * (tp(:,p1)-tp(:,p2)) ./ ...
      %    (norm(H(:,1)).*norm(tp(:,p1)-tp(:,p2)));
      
      if 0
         m1=mean(strfsimboot(:,:,p1+1)-strfsimboot(:,:,1),2);
         e1=std(strfsimboot(:,:,p1+1)-strfsimboot(:,:,1),0,2);
         tp1=shrinkage(m1-dcgp(1,p1+1),e1,0.5)+dcgp(1,p1+1);
         m2=mean(strfsimboot(:,:,p2+1)-strfsimboot(:,:,1),2);
         e2=std(strfsimboot(:,:,p2+1)-strfsimboot(:,:,1),0,2);
         tp2=shrinkage(m2-dcgp(1,p2+1),e2,0.5)+dcgp(1,p2+1);
         
         predsim(pidx)=xcorr(m1,m2,0,'coeff');
         %predsim(pidx)=xcorr(tp1,tp2,0,'coeff');
      elseif 1,
         predset=zeros(params.bootcount,1);
         for bootidx=1:params.bootcount
            predset(bootidx)=...
                xcorr(strfsimboot(:,bootidx,p1+1)-strfsimboot(:,bootidx,1),...
                      strfsimboot(:,bootidx,p2+1)-strfsimboot(:,bootidx,1),...
                      0,'coeff');    %  ,'coeff'
         end
         predsim(pidx)=median(predset);
      else
         predsim(pidx)=xcorr(tpred(:,p1),tpred(:,p2),0,'coeff');
      end
      
      dcsim(pidx)=((dcgp(1,1+p1)-mean(dcgp(1,2:end))) * ...
                   (dcgp(1,1+p2)-mean(dcgp(1,2:end))))./ ...
          var(dcgp(1,2:end));
   end
   
   
   figure(4);
   subplot(1,2,1);
   scatter(targsim,predsim);
   title(sprintf('%s targ v pred %.3f (psep<%.2f,ploc<%.2f)',cellid,...
                 xcov(targsim,predsim,0,'coeff'),...
                 pxc([4 5])));
   axis square
   
   subplot(1,2,2);
   scatter((predin+predout)/2,(predin-predout)/2);
   title(sprintf('predout v predin %.3f',...
                 xcov((predin+predout)/2,(predin-predout)/2,0,'coeff')));
   axis square
   
   
   %keyboard
else
   targsim=zeros(paircount,1);
   prefsim=zeros(paircount,1);
   predsim=zeros(paircount,1);
   dcsim=zeros(paircount,1);
   predin=zeros(4,1);
   predout=zeros(4,1);
end


if 0,
   
   iidx=[-(2:attcount)' zeros(attcount-1,1);
         (2:attcount)' zeros(attcount-1,1);
         (2:attcount)' ones(attcount-1,1)];
   % iidx=paircount+1;
   rowcount=size(iidx,1);
   colcount=11;
   
   PATCHONLY=0;
   
   stimextreme=zeros(spacecount,colcount,rowcount);
   if PATCHONLY,
      simdiff=zeros(length(patchpreds),rowcount);
   else
      simdiff=zeros(length(anyokidx),rowcount);
   end
   
   %keyboard
   
   for pidx=1:rowcount,
      i1=iidx(pidx,1);
      i2=iidx(pidx,2);
      if i1<0,
         i1=[2:(-i1-1) (-i1+1):attcount];
      elseif i2==1,
         i2=[2:i1-1 i1+1:attcount];
      end
      
      if PATCHONLY,
         if i2(1)==0,
            simdiff(:,pidx)=mean(patchpreds(:,i1),2);
         else
            simdiff(:,pidx)=mean(patchpreds(:,i1),2)-...
                mean(patchpreds(:,i2),2);
         end
         
      else
         if i2(1)==0,
            simdiff(:,pidx)=mean(strfsim(anyokidx,i1),2);
         else
            simdiff(:,pidx)=mean(strfsim(anyokidx,i1),2)-...
                mean(strfsim(anyokidx,i2),2);
         end
      end
      
      [crap,ttidx]=sort(-simdiff(:,pidx));
      
      if PATCHONLY,
         stimextreme(:,:,pidx)=...
             fpatches(:,ttidx([1 1:5 end-4:end])) -...
             repmat(mS(:,1),[1 colcount]);
      else
         stimextreme(:,:,pidx)=...
             bstim(anyokidx(ttidx([1 1:5 end-4:end])),:)' -...
             repmat(mS(:,1),[1 colcount]);
      end
   end
   
   if 0 & nargout==0,
      figure(4);
      clf
      showkern(stimextreme,kernfmt);
   elseif nargout==0
      figure(4);
      clf
      
      obins=linspace(0,180,15+1);
      obins=obins(1:end-1);
      sfbins=linspace(1,8,8);
      sformtx=sf2gr(stimextreme,length(obins),length(sfbins),0,0,kernfmt);
      sformtx=permute(sformtx,[2 1 3 4]);
      
      ttmax=max(abs(sformtx(:)));
      
      for ii=1:size(sformtx,3).*size(sformtx,4),
         subplot(size(sformtx,4),size(sformtx,3),ii);
         ta=repmat(sformtx(:,:,ii),[1 1 3])./ttmax;
         ta(:,:,1)=1+ta(:,:,1).*(ta(:,:,1)<0);
         ta(:,:,2)=1-abs(ta(:,:,2));
         ta(:,:,3)=1-ta(:,:,3).*(ta(:,:,3)>0);
         
         ta(find(ta>1))=1;
         ta(find(ta<0))=0;
         
         imagesc(obins,sfbins,ta);
         %imagesc(obins,sfbins,sformtx(:,:,ii),[-ttmax ttmax]);
         
         if ii<=size(sformtx,3).*(size(sformtx,4)-1),
            set(gca,'XTickLabel',[]);
         end
         axis xy
      end
   end
   
   colcount=size(stimextreme,2);
   chtargsim=zeros(rowcount,colcount-1);
   for attidx=1:rowcount,
      if iidx(attidx,1)>1,
         ttarg=mpatches(:,iidx(attidx,1))-mS(:,1);
      else
         ttarg=mpatches(:,attidx+1)-mS(:,1);
      end
      ttarg=ttarg./norm(ttarg);
      for sigidx=1:colcount-1,
         ttmatch=stimextreme(:,sigidx+1,attidx)./...
                 norm(stimextreme(:,sigidx+1,attidx));
         chtargsim(attidx,sigidx)=ttarg'*ttmatch;
      end
   end
   
   mchtargsim=[mean(chtargsim(:,1:(colcount-1)/2),2) ...
               mean(chtargsim(:,(colcount-1)/2+1:end),2)];
   
   if nargout==0,
      for pidx=1:rowcount,
         subplot(rowcount,colcount,(pidx-1).*colcount+1);
         plot(sort(simdiff(:,pidx)).*1000);
         title(sprintf('targ %d - %d',iidx(pidx,:)));
      
         subplot(rowcount,colcount,(pidx-1).*colcount+2);
         title(sprintf('max %d > %d',iidx(pidx,:)));
         
         subplot(rowcount,colcount,(pidx-1).*colcount+round(colcount/2+1));
         title(sprintf('max %d > %d',fliplr(iidx(pidx,:))));
      end 
      mchtargsim
   end
   
   %keyboard
   
end


if 0      
   nn=zeros(segcount,attcount);
   mm=zeros(segcount,attcount);
   ee=zeros(segcount,attcount);
   pp=zeros(segcount,attcount);
   for segidx=1:segcount,
      for attidx=1:attcount,
         segrange=anyokidx(ttidx(tedges(segidx):tedges(segidx+1)-1));
         segrange=segrange(~isnan(bresp(segrange,attidx)));
         nn(segidx,attidx)=length(segrange);
         mm(segidx,attidx)=mean(bresp(segrange,1));
         ee(segidx,attidx)=std(bresp(segrange,1))./sqrt(nn(segidx,attidx));
         pp(segidx,attidx)=mean(strfsim(segrange,1));
      end
   end
   
   figure(4);
   clf;
   sfmt={'k--','r-','g-','b-','k-'};
   for attidx=1:attcount,
      errorbar(pp(:,attidx),mm(:,attidx),ee(:,attidx),sfmt{attidx});
      hold on
   end
   hold off
   legend('a','1','2','3','4')
end   


hzscale=1000;
globaldc=globaldc.*hzscale; % or 60?
mdc=mean(globaldc,2);
edc=std(globaldc,1,2) .* sqrt(params.resampcount-1);
mgain=mean(globalgain,2);
egain=std(globalgain,1,2) .* sqrt(params.resampcount-1);

if nargout>0,
   res.kernfile=[rundata(1).respath,rundata(1).kernfile,'.gz'];
   res.cellid=cellid;
   
   % save mean response
   res.mresp=nanmean(squeeze(bresp));
   res.ncount=ncount;
   
   % save predxc
   res.predxc=predxc;
   res.pxc=pxc;
   res.predxccross=squeeze(predxccross(1,:,:));
   res.predxccrossrand=squeeze(predxccross(2,:,:));
   res.randxc=randxc;
   
   % predicted baseline response to targets, mean response, rank of targets
   res.strfresp=targpreds;
   res.strfrespnodc=targpredsnodc;
   res.strfresppos=targpredspos;
   res.strfrespneg=targpredsneg;
   res.strfmean=mean(strfsim);
   res.dcerr=(mdc'-mean(mdc))./(edc+(edc==0))';
   for aa=1:size(strfsim,2),
      %tt=sort(strfsim(:,aa));
      tt=sort(patchpreds(:,aa));
      for bb=1:size(targpreds,1),
         res.strfrank(bb,aa)=...
             max([0; find(targpreds(bb,aa)>tt)])./length(tt);
      end
      %res.strfrank(:,aa)=(targpreds(:,aa)-mean(tt))./std(tt);
   end
   
   % sample srf from each att condition
   res.Hset=squeeze(H2);
   
   res.targsim=targsim;
   res.prefsim=prefsim;
   res.predsim=predsim;
   res.dcsim=dcsim;
   res.predout=predout;
   res.predin=predin;
   
   return
end

ap=fpatches;

%tp=mpatches-repmat(mS(:,1),[1 size(mpatches,2)]);
tp=mpatches-repmat(mean(fpatches,2),[1 size(mpatches,2)]);
%tp=[mpatches(:,1) ...
%    mpatches(:,2:end)-repmat(mpatches(:,1),[1 size(mpatches,2)-1])];

% normalize for display
ap = ap ./ repmat(mean(ap,1),size(ap,1),1);
tp = tp ./ repmat(mean(abs(tp))+(mean(abs(tp))==0),size(tp,1),1);
mp=repmat(mean(ap,2),1,size(tp,2));

if size(tp,1)>spacecount,
   targmtx=rand(spacecount,1,attcount);
elseif size(tp,2)<attcount,
   targmtx=reshape(tp,spacecount,1,size(tp,2));
   targmtx(:,1,size(tp,2)+1:attcount)=0;
   targmtx(1,1,size(tp,2)+1:attcount)=1;
   targmtx=abs(targmtx).^0.5 .* sign(targmtx);
else
   targmtx=reshape(tp(:,1:attcount),spacecount,1,attcount);
   %targmtx=abs(targmtx).^0.5 .* sign(targmtx);
end

for attidx=1:attcount,
   targmtx(:,:,attidx)=targmtx(:,:,attidx) ./ ...
       max(abs(targmtx(:,:,attidx))) .* max(max(abs(H2(:,:,attidx))));
end

fprintf('pred xc : indiv att state MSE (att/rand):\n');
fprintf('attidx     :         ');
fprintf(' %11d',1:attcount);
fprintf('\n');

fprintf('framecount :         ');
fprintf(' %11d',ncount);
fprintf('\n');

for nlidx=1:nlcount,
   fprintf('%-11s:',nlnames{nlidx});
   fprintf('%5.2f/%5.2f',predxc(1:2,nlidx));
   fprintf(' (p<%4.2f):',pxc(nlidx));
   fprintf(' %5.2f/%5.2f',predxccross(1:2,nlidx,2:attcount));
   fprintf('\n');
end

%xcboot(find(isnan(xcboot)))=0;
%[mean(xcboot(1,:,:),3); mean(mean(xcboot(2:end,:,:),3))]

figure(2);
clf

eigmtx=cat(2,H2,targmtx,H2,Hd);
eigmtx(:,:,size(eigmtx,3)+1)=nan;
colcount=size(eigmtx,2);
rowcount=size(eigmtx,3);

showSFOR=1;
if showSFOR,
   obins=linspace(0,180,15+1);
   obins=obins(1:end-1);
   sfbins=linspace(1,8,8);
   sformtx=sf2gr(eigmtx,length(obins),length(sfbins),0,0,params.kernfmt);
   sformtx=permute(sformtx,[2 1 3 4]);
   %sformtx=flipdim(sformtx,1);
   
   ttmax=max(max(max(abs(sformtx(:,:,2,:)))));
   thmax=max(max(max(max((sformtx(:,:,[3 4],:))))));
   for rr=1:rowcount-1,
      for cc=1:colcount,
         subplot(rowcount,colcount,(rr-1)*colcount+cc);
         ta=sformtx(:,:,cc,rr);
         xx=-2:2;
         gsf=exp(-(xx./0.8).^2/2);
         gsf=gsf./sum(gsf(:))';
         gor=exp(-(xx./1.0).^2/2);
         gor=(gor./sum(gor(:)));
         
         astd=std(ta(:));
         
         smta=conv2(ta,gsf,'same');
         smta=cconv2(smta,gor);
         
         if cc==2,
            mm=ttmax;
         elseif ismember(cc,[1 3 4]),
            mm=thmax;
         else
            mm=max(max(abs(sformtx(:,:,1,1))));
            %mm=thmax;
         end
         if mm==0,
            mm=1;
         end
         
         ta=repmat(ta,[1 1 3])./mm;
         ta(:,:,1)=1+ta(:,:,1).*(ta(:,:,1)<0);
         ta(:,:,2)=1-abs(ta(:,:,2));
         ta(:,:,3)=1-ta(:,:,3).*(ta(:,:,3)>0);
         
         ta(find(ta>1))=1;
         ta(find(ta<0))=0;
         
         %keyboard
         imagesc(obins,sfbins,ta);
         %imagesc(obins,sfbins,ta,[-mm mm]);
         hold on
         slev1=1.0;
         slev2=2.0;
         contour(obins,sfbins,smta,[-astd*slev2  -astd*slev2],'k-');
         contour(obins,sfbins,smta,[-astd*slev1  -astd*slev1],'k--');
         contour(obins,sfbins,smta,[ astd*slev1   astd*slev1],'k--');
         contour(obins,sfbins,smta,[ astd*slev2   astd*slev2],'k-');
         hold off
         
         axis xy
         if rr<rowcount-1,
            set(gca,'XTickLabel',[]);
         else
            xlabel('orientation');
         end
         if cc>1,
            set(gca,'YTickLabel',[]);
         else
            ylabel('spatial freq');
         end   
      end
   end
   
   mmax=max(max(max(bigpatches(:,:,targlist))));
   for rr=2:rowcount-1,
      subplot(rowcount,colcount,rr.*colcount-colcount+1);
      
      ta=repmat(bigpatches(:,:,targlist(rr-1)),[1 1 3])./mmax;
      imagesc(ta);
      axis off
      axis image
      xlabel(sprintf('attidx %d',rr-1));
   end
   colormap(redblue);
   
   Ta=squeeze(sformtx(:,:,2,1:5));
   Ha=squeeze(sformtx(:,:,3,1:5));
   figure(3);
   clf
   
   %for pidx=1:length(pairidx),
   rc2=attcount-1;
   for pidx=1:attcount-1,
      
      if 0,
         p1idx=pairidx(pidx,1)+1;
         p2idx=pairidx(pidx,2)+1;
         
         Td=Ta(:,:,p1idx)-Ta(:,:,p2idx);
         Hd=Ha(:,:,p1idx)-Ha(:,:,p2idx);
      else
         p1idx=pidx+1;
         Td=Ta(:,:,p1idx);
         Hd=Ha(:,:,p1idx); % -Ha(:,:,1);
      end
      
      subplot(rc2,4,pidx*4-3);
      imagesc(Td,[-1 1].*(max(abs(Td(:))) + (max(abs(Td(:)))==0) ) );
      axis xy
      title(sprintf('targ %d - mean',pidx));
      
      subplot(rc2,4,pidx*4-2);
      imagesc(Hd,[-1 1].*(max(abs(Hd(:))) + (max(abs(Hd(:)))==0) ) );
      axis xy
      title(sprintf('srf %d - base',pidx));
      
      
      if 1,
         Hsf0=mean(Ha(:,:,1),2)./std(mean(Ha(:,:,1),2));
         Hor0=mean(Ha(:,:,1),1)./std(mean(Ha(:,:,1),1));
         
         Tsf=mean(Td,2)./(std(mean(Td,2)) + (std(mean(Td,2))==0) );
         Tor=mean(Td,1)./(std(mean(Td,1)) + (std(mean(Td,1))==0) );
         Hsf=mean(Hd,2)./std(mean(Ha(:,:,1),2));
         Hor=mean(Hd,1)./std(mean(Ha(:,:,1),1));
      else
         [u,s,v]=svd(Td);
         if sum(Td,2)'*u(:,1)>0,
            Tsf=u(:,1);
            Tor=v(:,1);
         else
            Tsf=-u(:,1);
            Tor=-v(:,1);
         end
         [u,s,v]=svd(Hd);
         if sum(Hd,2)'*u(:,1)>0,
            Hsf=u(:,1);
            Hor=v(:,1);
         else
            Hsf=-u(:,1);
            Hor=-v(:,1);
         end
      end
      
      subplot(rc2,4,pidx*4-1);
      plot(Tsf,'k-'); hold on; 
      plot(Hsf0,'b--'); plot(Hsf,'b-'); hold off
      title('sf');
      
      subplot(rc2,4,pidx*4-0);
      plot(Tor,'k--'); hold on; 
      plot(Hor0,'b--'); plot(Hor,'b-'); hold off
      title('or');
   end
   colormap(redblue);
   
   figure(2);
   
else
   showkern(eigmtx,kernfmt,iconside,{},0);
end

maxr=max(bresp(:,1));
xx=linspace(-maxr,maxr,100);

if 0,
subplot(rowcount,colcount,1);
cla; axis off;
for ii=1:size(eigmtx,3),
   subplot(rowcount,colcount,ii*colcount);
   cla; axis off;
end

end


if 1
   
   subplot(rowcount,colcount,colcount);
   plot(targpreds(:,1).*hzscale,'k-');
   hold on
   plot([0 attcount],mean(globaldc(:)).*[1 1],'k--');
   tt=diag(targpreds(:,2:end)).*hzscale;
   plot(tt,'bx');
   %plot(targpreds(:,2:end));
   hold off
   
   a=axis;
   axis([0 attcount+1 a(3) a(4)]);
   xlabel('attidx');
   %sleg={'a'};
   %for aa=2:attcount,
   %   sleg{aa}=sprintf('%d',aa-1);
   %end
   %legend(sleg);
end

subplot(rowcount,colcount,2);
if 1,
   
   bincount=7;
   ppp=zeros(bincount,attcount-1);
   rrr=zeros(bincount,attcount-1);
   rrre=zeros(bincount,attcount-1);
   fmtstr={'k-','r-','b-','g-'};
   h=zeros(size(fmtstr));
   
   cla
   for ii=1:attcount-1,
      attuse=find(~isnan(bresp(:,ii+1)));
      rrb=round(linspace(1,length(attuse),bincount+1));
      pdata=sortrows([rprec(attuse,1) bresp(attuse,1)]).*1000;
      
      for bb=1:bincount,
         ppp(bb,ii)=mean(pdata(rrb(bb):(rrb(bb+1)-1),1));
         rrr(bb,ii)=mean(pdata(rrb(bb):(rrb(bb+1)-1),2));
         %rrre(bb,ii)=std(diff(pdata(rrb(bb):(rrb(bb+1)-1),:),[],2))./ ...
         %    sqrt(rrb(bb+1)-1-rrb(bb));
         rrre(bb,ii)=std(pdata(rrb(bb):(rrb(bb+1)-1),2)) ./ ...
             sqrt(rrb(bb+1)-1-rrb(bb));
      end
      ht=errorbar(ppp(:,ii),rrr(:,ii),rrre(:,ii),fmtstr{ii});
      h(ii)=ht(1);
      hold on
   end
   hold off
   title('lin vs att resp');
   a=axis;
   axis([min(ppp(:))-max(ppp(:)).*0.1 max(ppp(:)).*1.1 ...
         min(rrr(:)-rrre(:))-max(rrr(:)+rrre(:)).*0.1 max(rrr(:)+rrre(:)).*1.1]);
   legend(h,'1','2','3','4');
   
   subplot(rowcount,colcount,3);
   cla
   axis off
else
   
   %plot(linspace(0,attcount,params.bootcount),globaldc','--');
   plot([0 attcount],mean(globaldc(:)).*[1 1],'k--');
   hold on
   ht=errorbar(mdc,edc,'k+');
   set(ht,'LineWidth',2);
   hold off
   title('att vs dc');
   a=axis;
   axis([0 attcount a(3) a(4)]);
   set(gca,'XTickLabel',[]);
   
   subplot(rowcount,colcount,3);
   
   %plot(linspace(0,attcount,params.bootcount),globalgain','--');
   plot([0 attcount],mean(globalgain(:)).*[1 1],'k--');
   hold on
   ht=errorbar(mgain,egain,'k+');
   set(ht,'LineWidth',2);
   hold off
   title('att vs gain');
   a=axis;
   axis([0 attcount a(3) a(4)]);
   set(gca,'XTickLabel',[]);
end

subplot(rowcount,colcount,1);
title(sprintf('%s no-att kernel',cellid));
for attidx=2:attcount,
   subplot(rowcount,colcount,attidx*colcount-colcount+1);
   ylabel(sprintf('attidx=%d',attidx-1));
end


subplot(rowcount,1,rowcount);
axis off;
axis ij;

tpredxc=predxc;
if tpredxc(2,1)>0,
   tpredxc(1:2,2)=tpredxc(1:2,2)./tpredxc(2,2).*tpredxc(1,1);
   tpredxc(1:2,3)=tpredxc(1:2,3)./tpredxc(2,3).*tpredxc(1,1);
   tpredxc(1:2,4)=tpredxc(1:2,4)./tpredxc(2,4).*tpredxc(1,1);
   tpredxc(1:2,end)=tpredxc(1:2,end)./tpredxc(2,end).*tpredxc(1,end-1);
end

tfracxc=tpredxc(1,[1 2 4 5]);
tfracxc(find(tfracxc<0))=0;
fracxc=zeros(size(tfracxc));
         
fracmax=max(tfracxc);
if fracmax==0,
   tfracxc(1)=1;
else
   if tfracxc(2)>tfracxc(1),
      fracxc(2)=(tfracxc(2).^2-tfracxc(1).^2)./fracmax.^2;
   else
      tfracxc(2)=tfracxc(1);
   end
   if tfracxc(3)>tfracxc(2),
      fracxc(3)=(tfracxc(3).^2-tfracxc(2).^2)./fracmax.^2;
   else
      tfracxc(3)=tfracxc(2);
   end
   if tfracxc(4)==fracmax,
      fracxc(4)=(tfracxc(4).^2-tfracxc(3).^2)./fracmax.^2;
   end
   fracxc(1)=1-sum(fracxc(2:end));
end
fracxc=[fracxc(1:3) 0 fracxc(4)];

details={};
details{1}=sprintf('%s: pred xc: indiv att state MSE (att/rand):',cellid);
details{2}=[sprintf('attidx    :            FR'),...
            sprintf('%11d ',1:attcount)];
details{3}=[sprintf('framecount:              '),...
            sprintf('%11d ',ncount)];

for nlidx=1:nlcount,
   details{3+nlidx}=[...
      sprintf('%-10s:',nlnames{nlidx}),...
      sprintf('%5.2f/%5.2f %.2f',predxc(1:2,nlidx),fracxc(nlidx)),...
      sprintf(' (p<%4.2f):',pxc(nlidx)),...
      sprintf(' %5.2f/%5.2f',predxccross(1:2,nlidx,2:attcount))];
end

ht=text(-0.05,0.1,details,'VerticalAlignment','top','FontName','Courier');
fullpage('portrait');


if strcmp(showextra,'bigpatches'),
   if (params.runclassid==3 | params.runclassid==10) & ~showSFOR,
      if params.runclassid==10,
         targlist=1:size(bigpatches,3);
      end
      
      mmin=min(min(min(bigpatches(:,:,targlist))));
      mmax=max(max(max(bigpatches(:,:,targlist))));
      
      figure(4);
      clf
      for ii=1:length(targlist),
         subplot(length(targlist)+2,1,ii);
         imagesc(bigpatches(:,:,targlist(ii)),[mmin mmax]);
         title(sprintf('%s target %d (id %d)',cellid,ii,targlist(ii)));
         axis off
         axis image
      end
      colormap(gray);
      fullpage('portrait');
   end
elseif strcmp(showextra,'srfdiff'),
   
   figure(4);
   clf
   
   eigmtx=repmat(H2,[1 attcount]);
   for rr=2:attcount,
      eigmtx(:,rr,1)=H2(:,1,rr);
      for cc=2:attcount,
         eigmtx(:,rr,cc)=H2(:,1,rr)-H2(:,1,cc);
      end
   end
   colcount=size(eigmtx,2);
   rowcount=size(eigmtx,3);
   
   
   if showSFOR,
      obins=linspace(0,180,15+1);
      obins=obins(1:end-1);
      sfbins=linspace(1,8,8);
      sformtx=sf2gr(eigmtx,length(obins),length(sfbins),0,0,params.kernfmt);
      sformtx=permute(sformtx,[2 1 3 4]);
      %sformtx=flipdim(sformtx,1);
      
      ttmax=max(abs(sformtx(:)));
      %ttmax=max(max(max((sformtx(:,:,2,:)))));
      for rr=1:rowcount,
         for cc=1:colcount,
            subplot(rowcount,colcount,(rr-1)*colcount+cc);
            ta=sformtx(:,:,cc,rr);
            xx=-2:2;
            gsf=exp(-(xx./0.8).^2/2);
            gsf=gsf./sum(gsf(:))';
            gor=exp(-(xx./1.0).^2/2);
            gor=(gor./sum(gor(:)));
            
            astd=std(ta(:));
            
            smta=conv2(ta,gsf,'same');
            smta=cconv2(smta,gor);
            
            ta=repmat(ta,[1 1 3])./ttmax;
            ta(:,:,1)=1+ta(:,:,1).*(ta(:,:,1)<0);
            ta(:,:,2)=1-abs(ta(:,:,2));
            ta(:,:,3)=1-ta(:,:,3).*(ta(:,:,3)>0);
            
            if sum(ta(:)>1 | ta(:)<0)
               keyboard
            end
            
            ta(find(ta>1))=1;
            ta(find(ta<0))=0;
            
            %keyboard
            imagesc(obins,sfbins,ta);
            %imagesc(obins,sfbins,ta,[-mm mm]);
            hold on
            slev1=1.0;
            slev2=2.0;
            contour(obins,sfbins,smta,[-astd*slev2  -astd*slev2],'k-');
            contour(obins,sfbins,smta,[-astd*slev1  -astd*slev1],'k--');
            contour(obins,sfbins,smta,[ astd*slev1   astd*slev1],'k--');
            contour(obins,sfbins,smta,[ astd*slev2   astd*slev2],'k-');
            hold off
            
            axis xy
            if (rr==1 & cc==1),
               title(cellid);
            end
            if rr<rowcount,
               set(gca,'XTickLabel',[]);
            else
               xlabel('orientation');
            end
            if cc>1,
               set(gca,'YTickLabel',[]);
            else
               ylabel('spatial freq');
            end   
         end
      end
      
      colormap(redblue);
      set(gcf,'PaperOrientation','portrait',...
              'PaperPosition',[1.25 2.5 6 6]);
   else
      showkern(eigmtx,kernfmt,iconside,{},0);
   end
end

return

%
% compare preds by different att states
%
if exist('attxc','var'),
   disp('local att preds');
   for nlidx=1:nlusecount,
      fprintf('nlidx=%d:\n',nluse(nlidx));
      for attidx=2:attcount,
         pin=attxc(attidx,attidx,nlidx);
         pout=nanmean(attxc([2:attidx-1 attidx+1:end],attidx,nlidx));
         fprintf('attidx=%d: %6.3f/%6.3f',attidx,pin,pout);
         
         for att2=2:attcount,
            if 0 & att2==attidx,
               fprintf('  x.xxx/ x.xxx');
            else
               p1=valattpreds(:,attidx,nlidx);
               p2=valattpreds(:,att2,nlidx);
               a1idx=find(~isnan(bresp(anyokidx,1,attidx)));
               a2idx=find(~isnan(bresp(anyokidx,1,att2)));
               
               pin=xcov(bresp(anyokidx([a1idx;a2idx]),:,1),...
                        [p1(a1idx);p2(a2idx)],0,'coeff');
               pout=xcov(bresp(anyokidx([a1idx;a2idx]),:,1),...
                         [p2(a1idx);p1(a2idx)],0,'coeff');
               
               fprintf(' %6.3f/%6.3f',pin,pout);
            end
         end
         fprintf('\n');
      end   
   end
end

tstrf=vstrf(1);
clear tunedata
for attidx=1:attcount,
   tstrf.h=H2(:,:,attidx);
   tunedata(attidx)=kern2tune(tstrf);
end

figure(4);
clf
obins=tunedata(1).obins;
sfbins=tunedata(1).sfbins;

patorsf=sf2gr(mpatches,length(obins),length(sfbins),1,0,kernfmt);
for attidx=1:attcount
   
   if attidx>1,
      subplot(attcount+1,4,attidx*4-3);
      
      ta=patorsf(:,:,attidx);
      ta=ta.*(ones(length(obins),1)*sqrt(sfbins));
      ta=ta.^2;
      
      mrange=max(abs(ta(:)));
      imagesc(obins,sfbins,ta',[-mrange mrange]);
      axis xy
      ylabel(sprintf('attidx=%d or',attidx));
      if attidx==2,
         title(sprintf('%s targets',cellid));
      end
   end
   
   subplot(attcount+1,4,attidx*4-2);
   mrange=max(max(abs(tunedata(attidx).orsf(:,:,1))));
   imagesc(obins,sfbins,tunedata(attidx).orsf(:,:,1)',[-mrange mrange]);
   
   if 1,
      xx=-2:2;
      gsf=exp(-(xx./0.8).^2/2);
      gsf=gsf./sum(gsf(:))';
      gor=exp(-(xx./1.0).^2/2);
      gor=(gor./sum(gor(:)));
      
      ta=tunedata(attidx).orsf(:,:,1)';
      astd=std(ta(:));
      
      ta=conv2(ta,gsf,'same');
      ta=cconv2(ta,gor);
      
      slev1=1.0;
      slev2=2.0;
      
      hold on
      contour(obins,sfbins,ta,[-astd*slev2  -astd*slev2],'k-');
      contour(obins,sfbins,ta,[-astd*slev1  -astd*slev1],'k--');
      contour(obins,sfbins,ta,[ astd*slev1   astd*slev1],'k--');
      contour(obins,sfbins,ta,[ astd*slev2   astd*slev2],'k-');
      hold off
   end
   
   axis xy
   if attidx==1,
      title(sprintf('%s spatial tuning',cellid));
   end
   
   subplot(attcount+1,4,attidx*4-1);
   plot(tunedata(attidx).obins,tunedata(attidx).por(:,1));
   ostep=(tunedata(attidx).obins(2)-tunedata(attidx).obins(1))./2;
   axis([tunedata(attidx).obins(1)-ostep ...
         tunedata(attidx).obins(end)+ostep -1 1]);
   
   if attidx==1,
      title('orientation');
   end
   
   subplot(attcount+1,4,attidx*4);
   plot(tunedata(attidx).sfbins,tunedata(attidx).psf(:,1));
   sfstep=(tunedata(attidx).sfbins(2)-tunedata(attidx).sfbins(1))./2;
   axis([tunedata(attidx).sfbins(1)-sfstep ...
         tunedata(attidx).sfbins(end)+sfstep -1 1]); 
   
   if attidx==1,
      title('spatial freq');
   end
end
colormap(redblue);

subplot(attcount+1,1,attcount+1);
axis off;
axis ij;
ht=text(0,0,details,'VerticalAlignment','top','FontName','Courier');
fullpage('landscape');


return


if OUTNLMODE==4,
   
   resplen=size(bresp,1);
   nluse=[4 5];
   nlusecount=length(nluse);
   
   attpreds=zeros(resplen,attcount,nlusecount,params.bootcount);
   attxc=zeros(attcount,attcount,nlusecount);
   
   for attidx=1:attcount,
      for nlidx=1:nlusecount,
         attpdummy=zeros(resplen,1);
         for bootidx=1:params.bootcount,
            if attidx==1,
               tstrf=vstrf(1,1,1,bootidx);
            else
               tstrf=vstrf(nluse(nlidx),attidx-1,1,bootidx);
            end
            
            tstim=bstim'-repmat(tstrf.mS,[1 resplen]);
            tpred=kernpredict(tstrf.h,tstim,1,0);
            
            tnltype=tstrf.nltype;
            tnlparms=tstrf.nlparms;
            if ~isempty(tnltype) & ~strcmp(tnltype,'none'),
               attpreds(:,attidx,nlidx,bootidx)=feval(tnltype,tnlparms,tpred);
            else
               attpreds(:,attidx,nlidx,bootidx)=tpred;
            end
            
            if nlidx==2,
               attpreds(:,attidx,nlidx,bootidx)=...
                   attpreds(:,attidx,nlidx,bootidx) + ...
                   attpreds(:,attidx,1,bootidx);
            end
            
            a0idx=vidx((round((bootidx-1)/params.bootcount*vcount)+1):...
                       round(bootidx/params.bootcount*vcount));
            attpdummy(a0idx)=attpreds(a0idx,attidx,nlidx,bootidx);
         end
         
         for att2=1:attcount,
            if att2==1,
               tokidx=anyokidx;
            else
               tokidx=find(~isnan(bresp(:,att2)));
            end
            attxc(attidx,att2,nlidx)=...
                xcov(attpdummy(tokidx),bresp(tokidx,1),0,'coeff');
         end
      end
   end
   attpreds(:,2:end,2,:)=attpreds(:,2:end,1,:)+attpreds(:,2:end,2,:);
   attpreds=attpreds(anyokidx,:,:,:);
else
   
   
   resplen=length(anyokidx);
   nluse=[3 6];
   nlusecount=length(nluse);
   
   attpreds=zeros(resplen,attcount,attcount-1,nlusecount,params.bootcount);
   
   for attidx=1:attcount,
      for nlidx=1:nlusecount,
         for bootidx=1:params.bootcount,
            if attidx==1,
               tstrf=vstrf(1,1,1,bootidx);
            else
               tstrf=vstrf(nluse(nlidx),attidx-1,1,bootidx);
            end
            
            tstim=bstim(anyokidx,:)'-repmat(tstrf.mS,[1 resplen]);
            tpred=kernpredict(tstrf.h,tstim,1,0);
            
            tnltype=tstrf.nltype;
            for att2idx=1:attcount-1,
               tnlparms=tstrf.nlparmsalt{att2idx};
               if ~isempty(tnltype) & ~strcmp(tnltype,'none'),
                  attpreds(:,attidx,att2idx,nlidx,bootidx)=feval(tnltype,tnlparms,tpred);
               else
                  attpreds(:,attidx,att2idx,nlidx,bootidx)=tpred;
            end
            end
         end
      end
   end
end


figure(4);
clf

targpredsraw=zeros(size(mpatches,2),params.bootcount,attcount);
for bootidx=1:params.bootcount,
   for attidx=1:attcount
      ts=mpatches-repmat(mS(:,bootidx,attidx),[1 size(mpatches,2)]);
      
      targpredsraw(:,bootidx,attidx)=ts'*H(:,bootidx,attidx);
      
      targpredsraw(:,bootidx,attidx)=...
          dcgain(vstrf(1,1,1,bootidx).nlparms,targpredsraw(:,bootidx,attidx));
   end
end

plotcount=attcount-1;
linstr={'k:','b-','r-','g-','k-'};
mp=mean(attpreds,4);
for attidx=1:plotcount,
   
   aokidx=find(~isnan(bresp(anyokidx,1,attidx+1)));
   %aokidx=1:length(anyokidx);
   
   % 1H, att dcg: mp(aokidx,attidx+1,attidx,1),...
   % 4H, att dcg: mean(mp(aokidx,[2:attidx attidx+2:attcount],attidx,2),2),...
   
   if 1,
      smp=cat(2,mp(aokidx,1),...
              mp(aokidx,attidx+1,1),...
              mp(aokidx,attidx+1,2));
   else
      smp=zeros(length(aokidx),3);
      for attidx2=[1:attidx-1 attidx+1:attcount-1],
         smp(:,1)=smp(:,1) + mp(aokidx,attidx2+1,attidx2,2);
         smp(:,2)=smp(:,2) + mp(aokidx,attidx2+1,attidx,2);
      end
      smp(:,1:2)=smp(:,1:2)./(attcount-2);
      smp(:,3)=mp(aokidx,attidx+1,attidx,2);
   end
   
   tr=bresp(anyokidx(aokidx),attidx+1);
   
   if 1,
      [oo,ss]=sort(smp(:,1));
      
      smp=smp(ss,:);
      smp(:,end)=gsmooth(smp(:,end),5);
      %smp(:,2:end)=gsmooth(smp(:,2:end),3);
      
      tr=gsmooth(tr(ss),5);
   else
      [tr,ss]=sort(tr);
      
      smp=smp(ss,:);
      smp=gsmooth(smp,3);
   end
   
   subplot(2,2,attidx);
   for ii=size(smp,2):-1:1,
      plot(smp(:,ii),linstr{ii});
      hold on
   end
   %plot(tr,'r:');
   
   plot([0,round(size(smp,1)/3)],[1 1].*mean(targpredsraw(attidx+1,:,1),2),'r:');
   
   hold off
   title(sprintf('attidx=%d',attidx));
   
   axis([0 size(smp,1) 0 max(max(max(mean(attpreds,4))))]);
end

legend('noatt','dcg','local',2);



keyboard


return







figure(4);
clf
subplot(2,1,1);
plotstyle={'k:','r-','g--','b-','m--'};
attstate={'all','1','2','3','4'};
xx=linspace(nanmin(rprec(:,1)),nanmax(rprec(:,1)));
for attidx=1:attcount,
    plot(xx,sigmoid(sigparms(:,attidx),xx),plotstyle{attidx});
    hold on
end

subplot(2,1,2);
for attidx=1:attcount,
    errorbar((1:4)',sigparms(:,attidx),esigparms(:,attidx),plotstyle{attidx});
    hold on
end




if 0 & ~exist('asigparms','var'),
   spcount=3;
   asigparms=zeros(spcount,params.bootcount,attcount-1,spcount+1);
   valpred=zeros(size(bresp,1),spcount+1);
   valxc=zeros(1,size(asigparms,1)+1);
   
   % find baseline (all att) sigmoid parms
   fprintf('fitting hinge4, attidx=1');
   attidx=1;
   tokidx=find(sum(~isnan(bresp(:,1,2:end)),3)>0);
   for jj=1:params.bootcount,
      % use same jackknife sets as for fitting 
      useidx=vidx([1:round((jj-1)/params.bootcount*vcount) ...
                   round(jj/params.bootcount*vcount+1):end]);
      useidx=useidx(find(~isnan(bresp(useidx,1,attidx))));
      predidx=vidx(round((jj-1)/params.bootcount*vcount+1): ...
                   round(jj/params.bootcount*vcount));
      predidx=predidx(find(~isnan(bresp(predidx,1,attidx))));
      
      ract=bresp(useidx,1,1);
      rpred=rprec(useidx,1);
      
      tfitparms=fithinge4(rpred,ract,0);
      asigparms(:,jj,:,1)=repmat(tfitparms,[1 1 attcount-1]);
      
      valpred(predidx,1)=hinge4(asigparms(:,jj,attidx),...
                                rprec(predidx,1));
   end
   
   for sigidx=1:spcount,
      
      filler=ones(size(bresp,1),spcount);
      for fitidx=1:attcount-1,
         filler(find(~isnan(bresp(:,fitidx+1))),sigidx)=fitidx;
      end
      
      for jj=1:params.bootcount,
         % use same jackknife sets as for fitting 
         useidx=vidx([1:round((jj-1)/params.bootcount*vcount) ...
                      round(jj/params.bootcount*vcount+1):end]);
         predidx=vidx(round((jj-1)/params.bootcount*vcount+1): ...
                      round(jj/params.bootcount*vcount));
         
         ract=bresp(useidx,1,1);
         rpred=rprec(useidx,1);
         
         tfitparms=fithinge4(rpred,ract,0,filler(useidx,:));
         
         for attidx=1:attcount-1,
            useidx=[1:sigidx-1 sigidx+attidx-1 ...
                    (sigidx+attcount-1):(spcount+attcount-2)]';
            
            asigparms(:,jj,attidx,sigidx+1)=tfitparms(useidx);
            
            apredidx=predidx(find(~isnan(bresp(predidx,1,attidx+1))));
            valpred(apredidx,sigidx+1)=...
                hinge4(asigparms(:,jj,attidx,sigidx+1),rprec(apredidx,1));
         end
      end
   end
   for sigidx=1:spcount+1,
      valxc(sigidx)=xcov(bresp(anyokidx,1),valpred(anyokidx,sigidx),...
                         0,'coeff');
   end
   
   sigparms=squeeze(mean(asigparms,2));
   esigparms=squeeze(std(asigparms,1,2)).*sqrt(params.bootcount-1);
   nsigparms=sigparms./repmat(sigparms(:,1),[1 attcount-1 spcount+1]);
   nesigparms=esigparms./repmat(sigparms(:,1),[1 attcount-1 spcount+1]); 
elseif ~exist('asigparms','var'),
   spcount=4;
   asigparms=zeros(spcount,params.bootcount,attcount);
   valpred=zeros(size(bresp,1),spcount+1);
   valxc=zeros(1,size(asigparms,1)+1);
   
   % find baseline (all att) sigmoid parms
   fprintf('fitting sigmoid, attidx=1');
   attidx=1;
   tokidx=find(sum(~isnan(bresp(:,1,2:end)),3)>0);
   for jj=1:params.bootcount,
      % use same jackknife sets as for fitting 
      useidx=vidx([1:round((jj-1)/params.bootcount*vcount) ...
                   round(jj/params.bootcount*vcount+1):end]);
      useidx=useidx(find(~isnan(bresp(useidx,1,attidx))));
      predidx=vidx(round((jj-1)/params.bootcount*vcount+1): ...
                   round(jj/params.bootcount*vcount));
      predidx=predidx(find(~isnan(bresp(predidx,1,attidx))));
      
      ract=bresp(useidx,1,1);
      rpred=rprec(useidx,1);
      
      asigparms(:,jj,attidx)=fitsigmoid(rpred,ract,0);
      
      valpred(predidx,1)=sigmoid(asigparms(:,jj,attidx),...
                                 rprec(predidx,1));
   end
   
   for attidx=2:attcount,
      
      fprintf(' %d',attidx);
      
      tokidx=find(~isnan(bresp(:,1,attidx)));
      
      for jj=1:params.bootcount,
         for sigidx=1:spcount,
            % use same jackknife sets as for fitting 
            useidx=vidx([1:round((jj-1)/params.bootcount*vcount) ...
                         round(jj/params.bootcount*vcount+1):end]);
            useidx=useidx(find(~isnan(bresp(useidx,1,attidx))));
            
            predidx=vidx(round((jj-1)/params.bootcount*vcount+1): ...
                         round(jj/params.bootcount*vcount));
            predidx=predidx(find(~isnan(bresp(predidx,1,attidx))));
            
            ract=bresp(useidx,1,1);
            rpred=rprec(useidx,1);
            
            forcevalues=asigparms(:,jj,1);
            forcevalues(sigidx)=nan;
            
            tsigparms=fitsigmoid(rpred,ract,0,forcevalues);
            asigparms(sigidx,jj,attidx)=tsigparms(sigidx);
            
            valpred(predidx,sigidx+1)=sigmoid(tsigparms,rprec(predidx,1));
         end
      end
   end
   
   for sigidx=1:spcount+1,
      valxc(sigidx)=xcov(bresp(anyokidx,1),valpred(anyokidx,sigidx),...
                         0,'coeff');
   end
   
   sigparms=squeeze(mean(asigparms,2));
   esigparms=squeeze(std(asigparms,1,2)).*sqrt(params.bootcount-1);
   nsigparms=sigparms./repmat(sigparms(:,1),[1 attcount]);
   nesigparms=esigparms./repmat(sigparms(:,1),[1 attcount]);
else
   spcount=size(sigparms,1);
   
end
fprintf('\n');

figure(4);
clf

xmin=min(rprec(anyokidx,1))-std(rprec(anyokidx,1));
xmax=max(rprec(anyokidx,1))+std(rprec(anyokidx,1));
ymin=nanmin(bresp(:,1));
ymax=nanmax(bresp(:,1));
xrange=linspace(xmin,xmax,100);
plotstyle={'k:','r-','g--','b-','m--'};
attstate={'all','1','2','3','4'};
spstring={'x10','sigma','baseline','amp'};

for sigidx=1:spcount,
   subplot(spcount,1,sigidx);
   for attidx=1:attcount,
      sparms=sigparms(:,1);
      sparms(sigidx)=sigparms(sigidx,attidx);
      plot(xrange,sigmoid(sparms,xrange),plotstyle{attidx});
      hold on
   end
   
   for attidx=1:attcount,
      if attidx==1,
         tokidx=find(sum(~isnan(bresp(:,1,2:end)),3)>0);
      else
         tokidx=find(~isnan(bresp(:,1,attidx)));
      end
      
      plotcnt=min([100 length(tokidx)]);
      plot(rprec(tokidx(1:plotcnt),1),bresp(tokidx(1:plotcnt),1,1),...
           [plotstyle{attidx}(1),'.']);
   end
   hold off
   
   title(sprintf('%s att %s: valxc: %.3f->%.3f',cellid,...
                 spstring{sigidx},valxc(1),valxc(sigidx+1)));
end
legend(attstate{1:attcount});
fullpage('portrait');

if 0
subplot(2,1,2);

for attidx=2:attcount,
   
   errorbar(nsigparms(:,attidx),nesigparms(:,attidx),plotstyle{attidx});
   hold on
   
end
hold off
end


if 0
figure(4);
clf
for ii=1:attcount,
   subplot(attcount,1,ii);
   mm=max(max(abs(lth(:,:,ii))));
   if mm>0,
      imagesc(lth(:,:,ii)',[-mm mm]);
      title(sprintf('att %d local kernels',ii));
      axis image
   end
   axis off
end
colormap(gray);
fullpage('Portrait');
end


