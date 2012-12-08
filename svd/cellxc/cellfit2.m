% cellfit2.m
%
% given :  H, bstim, rstartidx, rstopidx
%
% meant to be a sort of plug-in replacement for cellfit, doing
% FET-inspired selection of the regularization parameter.
%
% results saved in cellfit2data
%


% should be set by cellxc
if ~exist('attidx'),
   attidx=1;
   attcount=1;
end
fprintf('cellfit2.m: for attidx=%d/%d\n',attidx,attcount);

if params.respfmtcode==0,
   tbincount=min([params.maxlag(2)+1 14]);
   latcount=1;
   respuse=1;
else
   tbincount=size(H,2);
   latcount=tbincount;
   respuse=respcount;
end
nlcount=3;
DOTHRESH=0;

%sigrange=exp(linspace(0,1.5,sffiltsigma));
sigrange=exp(linspace(log(0.9),log(2.0),params.sffiltsigma));
sigcount=length(sigrange);

if attidx==1,
   xc=zeros(params.sfscount,sigcount,latcount,nlcount,attcount);
   xcnl=zeros(params.nloutparm,latcount,nlcount,attcount);
   expxc=zeros(latcount,nlcount,attcount);
   sfsfit=zeros(latcount,nlcount,attcount);
   sigfit=zeros(latcount,nlcount,attcount);
   threshfit=zeros(params.nloutparm,latcount,nlcount,attcount);
   hf=zeros(spacecount,tbincount,nlcount,attcount);
end

if VCELLXC==1,
   fitfilecount=sum(fitfile(:,attidx)>0);
   fitlen=0;
   fituseidx=1;
   for fitidx=1:fitfilecount,
      if times(2).stop(fitidx)-times(2).start(fitidx)>fitlen
         fituseidx=fitidx;
      end
   end
   
   tfitstartframe=fitstartframe(fituseidx,attidx);
   
   % if scale_pix not defined, use crf parms
   tstimloadparms=stimloadparms;
   if strcmp(stimloadcmd,'loadimfile') & ...
         length(stimloadparms)>0 & stimloadparms{1}==0,
      tstimloadparms{1}=stimloadparms{3} .* ...
       stimcrfs(fitfile(fituseidx,attidx))./stimwindowcrf;
   end
   stim=feval(stimloadcmd,stimfiles{fitfile(fituseidx,attidx)},...
               tfitstartframe,fitstopframe(fituseidx,attidx),...
               tstimloadparms{:});
   
   % filter stimulus segment if selected
   if ~isempty(stimfiltercmd),
      stim=feval(stimfiltercmd,stim,stimfilterparms{:});
   end
   iconside=size(stim);
   tmovlen=iconside(end);
   iconside=iconside(1:(end-1));
   if length(size(stim))>2,
      % reshape all spatial dims into one
      stim=reshape(stim,prod(iconside),tmovlen);
   end
   stim=stim';
   
   resp=feval(resploadcmd,respfiles{fitfile(fituseidx,attidx)},...
              resploadparms{:});
   resp=resp(tfitstartframe:fitstopframe(fituseidx,attidx),:,:);
   
   % filter response (eg smooth, pick out ICs, PCs? currently not used)
   if ~isempty(respfiltercmd),
      resp=feval(respfiltercmd,resp,respfilterparms{:});
   end
   predlen=fitstopframe(fituseidx,attidx)-tfitstartframe+1;
   
   [rstartidx,rendidx]=resampsegs(resp,params.resampcount);
   if strcmp(boundary,'zero'),
      rstartidx=rstartidx-params.maxlag(2);
      rendidx=rendidx-params.maxlag(1);
      rstartidx(find(rstartidx<1))=1;
      rendidx(find(rendidx>blen))=blen;
   end
end

blen=size(stim,1);
bsmean=mean(stim,1);
mod_psth=zeros(blen,respuse,nlcount,params.sfscount,sigcount);
for resampidx=1:params.resampcount,
   fprintf('resampidx=%d',resampidx);
   userange=rstartidx(resampidx):rendidx(resampidx);
   if params.meansub,
      % subtract exploratory mean from stimulus to get DC level
      % set correctly
      tstim=stim(userange,:)-repmat(bsmean,[length(userange) 1]);
      %tstim=stim(userange,:)-repmat(mSall(:,r,attidx),[length(userange) 1]);
   else
      tstim=stim(userange,:);
   end
   for sfsidx=1:params.sfscount,
      fprintf('.');
      
      for sigidx=1:sigcount,
         
         % shrink across resamples as needed
         thf=H(:,(1:tbincount)-params.maxlag(1),sfsidx,1,resampidx);
         smm=mH(:,(1:tbincount)-params.maxlag(1),sfsidx,1,attidx);
         sms=eH(:,(1:tbincount)-params.maxlag(1),sfsidx,1,attidx);
         smd=abs(smm)./(sms+(sms==0)) ./ sigrange(sigidx);
         smd=(1-smd.^(-2));
         smd=smd.*(smd>0);
         smd(find(isnan(smd)))=0;
         thf=thf.*smd;
         
         % linear prediction, rectified?
         tmod_psth=kernpredict(thf,tstim',spacecount,0,1);
         valididx=find(~isnan(tmod_psth(:,1)));
         tmod_psth=tmod_psth(valididx,:);
         linpred=sum(tmod_psth,2);
         
         mod_psth(userange(valididx),:,1,sfsidx,sigidx)=linpred;
         mod_psth(userange(valididx),:,2,sfsidx,sigidx)=...
             linpred.*(linpred>0);
         mod_psth(userange(valididx),:,3,sfsidx,sigidx)=...
             sum(tmod_psth.*(tmod_psth>0),2);
      end
   end
   fprintf('\n');
   if exist('BATQUEUEID') & BATQUEUEID>0,
      dbsetqueue(BATQUEUEID,2000+resampidx*10);
   end
end

disp('calculating xc...');

for respidx=1:respuse,
   tgoodidx=find(~isnan(resp(:,respidx)));
   for nlidx=1:nlcount,
      for sfsidx=1:params.sfscount,
         for sigidx=1:sigcount,
            
            xc(sfsidx,sigidx,respidx,nlidx,attidx)=...
                xcov(mod_psth(tgoodidx,respidx,nlidx,sfsidx,sigidx),...
                     resp(tgoodidx,respidx),0,'coeff');
         end
      end
   end
   
   for nlidx=1:nlcount,
      % find (sfs,sig) with maximum correlation
      maxidx=min(find(xc(:,:,respidx,nlidx,attidx)==...
                      max(max(xc(:,:,respidx,nlidx,attidx)))));
      [sfsfituse,sigfituse]=ind2sub([params.sfscount,sigcount],maxidx);
      sfsfit(respidx,nlidx,attidx)=sfsfituse;
      sigfit(respidx,nlidx,attidx)=sigfituse;
      
      % generate linear filter
      thf=mH(:,(1:tbincount)-params.maxlag(1),sfsfituse,1,attidx);
      if params.resampcount>1,
         smm=mH(:,(1:tbincount)-params.maxlag(1),sfsfituse,1,attidx);
         sms=eH(:,(1:tbincount)-params.maxlag(1),sfsfituse,1,attidx) .* ...
             sigrange(sigfituse);
         smd=abs(smm)./(sms+(sms==0));
         smd=1-smd.^(-2);
         smd=smd.*(smd>0);
         smd(find(isnan(smd)))=0;
         hf(:,:,nlidx,attidx)=thf.*smd;
      else
         hf(:,:,nlidx,attidx)=thf;
      end

      % fit output nonlinearities for each hf
      for nloutidx=1:params.nloutparm,
         if nloutidx==1,
            % no rectification
            seplinpred=kernpredict(hf(:,:,nlidx,attidx),stim',spacecount,0);
            linpred=sum(seplinpred,2);
            tmod_psth=linpred;
            t0=min(tmod_psth(:));
         elseif nloutidx==2,
            % rectify summed output
            tmod_psth=linpred;
            [t0,res]=findthresh(tmod_psth(tgoodidx),resp(tgoodidx,1));
         elseif nloutidx==3,
            % rectify each spatial channel
            tmod_psth=seplinpred;
            [t0,res]=findthresh(tmod_psth(tgoodidx,:),resp(tgoodidx,1));
         end
         
         tmod_psth=thresh(t0,tmod_psth);
         threshfit(nloutidx,1,nlidx,attidx)=t0;
         xcnl(nloutidx,1,nlidx,attidx)=...
             xcov(resp(tgoodidx,1),tmod_psth(tgoodidx),0,'coeff');
         if nloutidx==nlidx,
            expxc(1,nlidx,attidx)=xcnl(nloutidx,1,nlidx,attidx);
         end
      end
   end
end
%for nlidx=2:nlcount,
%   xc(:,:,:,nlidx,attidx)=xc(:,:,:,1,attidx);
%   sfsfit(:,nlidx,attidx)=sfsfit(:,1,attidx);
%   sigfit(:,nlidx,attidx)=sigfit(:,1,attidx);
%   hf(:,:,nlidx,attidx)=hf(:,:,1,attidx);
%   threshfit(:,1,nlidx,attidx)=threshfit(:,1,1,attidx);
%   xcnl(:,1,nlidx,attidx)=xcnl(:,1,1,attidx);
%end





