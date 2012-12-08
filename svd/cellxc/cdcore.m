% xccore.m
%
% script to compute unbiased kernel from preprocessed stimulus and
% response matrices (ie, stim X x T and resp T x M supplied)
% 
% returns:
%
% mean and stderr on kernel estimates: (X x tbincount x M x sfscount)
% mH=mean(H,5);
% eH=std(H,1,5) .* sqrt((resampcount-1)/resampcount);
% 
% mean and stderr on raw (biased) XC: (X x tbincount x M x sfscount)
% mSR=mean(SR,4);
% eSR=std(SR,1,4).* sqrt((resampcount-1)/resampcount);
% 
% % mean stim/response across resamples
% mSall=mean(mS,3);
% mRall=mean(mR,2);
% ntot=n;
%
% created SVD 1/23/03 - ripped off of cellxc.m
%

disp([mfilename,':']);

%
% define resampling regimes here
%
% option to branch according to resampfmt
if params.resampcount==1 | length(resp)==0,
   % bootstrap resampling
   rstartidx=1;
   rendidx=size(resp,1);
elseif params.resampfmt==1,
   [rstartidx,rendidx]=resampsegs(resp,params.resampcount);
   
elseif params.resampfmt==2,
   % shuffled response resampling
   rstartidx=starttimes(fidx,attidx);
   rendidx=stoptimes(fidx,attidx);
   resp=resampnoise(resp,params.resampcount);
   rsize=[rsize,size(resp,3)];
   resp=resp(:,:); % reshape to 2D matrix
   
elseif params.resampfmt==3,
   % gamma/poisson spiking model resampling   
   rstartidx=1;
   rendidx=size(resp,1);
   respcount=1;
   rsize(2)=1;
   [mu,alpha,beta]=reversepoisson(resp);
end

mS=nanmean(stim);
sS=nanstd(stim);
sS(mS==1)=1;

if 0,
   % don't subtract mean from stim.  does this allow for fairer
   % comparison between stimulus classes (ie, torc vs speech)?
   mS=mS.*0;
   sS=nanmean(stim);
end

nanstim=find(sum(isnan(stim),2)>0);
if length(nanstim)>0,
   disp('NaN values found in stim, NaN-ing corresponding r bins');
   resp(nanstim)=nan;
end

mR=nanmean(resp);
mSall=repmat(mS',[1 respcount]);
mRall=mR;
stim0=stim-repmat(mS,[size(stim,1),1]);
resp=resp-repmat(mR,[size(resp,1) 1]);

% KLUDGE: zero out channels that are just way too weak??
if isfield(params,'runclassid') & ismember(params.runclassid,[1 4 42 58]),
   ss=std(stim0);
   ss=ss./median(ss);
   if params.runclassid==4,
       ff=find(ss<0.4 & mS~=1);
   else
       ff=find(ss<0.3 & mS~=1);
   end
   if ~isempty(ff),
      disp('zeroing out weak channels!!!');
      ff
   else
      disp('no weak channels to zero out');
   end
   stim0(:,ff)=0;
end

H=zeros(spacecount,diff(params.maxlag)+1,params.sfscount,...
        respcount,params.resampcount);

% force old algorithm
%expstim=stim0;
%expspacecount=spacecount;

% larger stim
params.allowneg=getparm(params,'allowneg',1);
params.basisfilter=getparm(params,'basisfilter',{[1]});
params.costfun=getparm(params,'costfun','mse');
if isfield(params,'basisfilter'),
   expstim=[];
   if iscell(params.basisfilter),
      basisfilter=params.basisfilter;
   else
      basisfilter={param.basisfilter};
   end
   % transform stimulus into basis set constraint for the STRF
   for gg=1:length(basisfilter),
      expstim=[expstim rconv2(stim0,(basisfilter{gg})')];
   end
   expspacecount=size(expstim,2);
else
   basisfilter={};
   gset={{[0.5 0.01],0},{[1.5 0.1],0}};
   expstim=[];
   for gg=1:length(gset),
      expstim=[expstim gsmooth(stim0,gset{gg}{:},0)];
   end
   expspacecount=size(expstim,2);
end

%
% new format. each successive channel chosen across jackknifes
% ADDED 2006-04-23
%
if isfield(params,'sfsstep'),
    iterpersfscount=params.sfsstep;
else
    iterpersfscount=6;
end

params.activeoutnl=getparm(params,'activeoutnl',0);
params.stepfrac=getparm(params,'stepfrac',30);
fprintf('stepfrac=%d\n',params.stepfrac);

for respidx=1:respcount,
   coeffsize=nanstd(resp(:,respidx))./sS ./ params.stepfrac;
   normexpstim=expstim.*repmat(coeffsize,[length(expstim) 1]);
   
   tresp=repmat(resp(:,respidx),[1 params.resampcount]);
   predsofar=zeros(size(tresp));
   fprintf('\nrespidx=%d\n',respidx);
   
   if respidx>1,
      disp('restricting fit to subspace fit by attidx==1');
      %keyboard
      ttH=sum(sum(sum(H(:,:,:,1,:),3),5),2);
      ff=find(ttH==0);
      ff
      expstim(:,ff)=0;
   end
   
   lasterror=repmat(nanstd(abs(tresp(:))).^2,[params.resampcount,1]).*100;
   for sfsidx=1:params.sfscount,
      % in order to get denser sampling without wasting a lot of time,
      % do multiple iterations per strf that will actually be tested
      % at the next stage
      
      for subidx=1:iterpersfscount,
         
         % can't use quick & dirty calculations because preds for each
         % resample create different residuals for each resample on
         % subsequent iterations!
         disp('computing separate cross correlation for all jackknifes');
         
         ttmse=zeros(expspacecount,diff(params.maxlag)+1,2,...
                      params.resampcount);
         USECC=1;
         for rridx=1:params.resampcount,
            tr=tresp(:,rridx);
            tr(rstartidx(rridx):rendidx(rridx),:)=nan;
            rgoodidx=find(~isnan(tr));
            
            for tt=params.maxlag(1):params.maxlag(2),
               trg=rgoodidx;
               trg=trg(find(trg-tt>0 & trg-tt<size(normexpstim,1)));
               if strcmp(params.costfun, 'mse'),
                  for xx=1:expspacecount,
                     % test increment of each channel
                     ttmse(xx,tt-params.maxlag(1)+1,1,rridx)=...
                         mean(abs(tr(trg)-normexpstim(trg-tt,xx)).^2);
                     % test decrement of each channel
                     ttmse(xx,tt-params.maxlag(1)+1,2,rridx)=...
                         mean(abs(tr(trg)+normexpstim(trg-tt,xx)).^2);
                     %ttmse(xx,tt-params.maxlag(1)+1,2,rridx)=...
                     %    mean((repmat(tr(rgoodidx),[1 expspacecount])+...
                     %          normexpstim(trg-tt,:)).^2);
                  end
               elseif strcmp(params.costfun, 'mi'),
                  % mutual info cost function, a la Tatyana
                  ttr=sqrt(tresp(trg)+mR(respidx));
                  
                  for xx=1:expspacecount,
                     ttpred=predsofar(trg,rridx)+normexpstim(trg-tt,xx);
                     tpstd=std(ttpred).*3;
                     if tpstd>0,
                        ttpred(find(ttpred-mean(ttpred)<-tpstd))=...
                            mean(ttpred)-tpstd;
                        ttpred(find(ttpred-mean(ttpred)>tpstd))=...
                            mean(ttpred)+tpstd;
                        
                        ttmse(xx,tt-params.maxlag(1)+1,1,rridx)=...
                            quick_mi(ttpred,ttr);
                     end
                     
                     ttpred=predsofar(trg,rridx)-normexpstim(trg-tt,xx);
                     tpstd=std(ttpred).*3;
                     if tpstd>0 & (sfsidx>1 | subidx>1),
                        ttpred(find(ttpred-mean(ttpred)<-tpstd))=...
                            mean(ttpred)-tpstd;
                        ttpred(find(ttpred-mean(ttpred)>tpstd))=...
                            mean(ttpred)+tpstd;
                        
                        ttmse(xx,tt-params.maxlag(1)+1,2,rridx)=...
                            quick_mi(ttpred,ttr);
                     end
                  end
                  
               elseif strcmp(params.costfun, 'cc')
                  ttmse(:,tt-params.maxlag(1)+1,1,rridx)=...
                      normexpstim(trg-tt,:)'*tr(trg);
               end
            end
            
            % and compute prediction for this new channel
            rpred=ones(size(tresp)).*nan;
            if strcmp(params.costfun, 'mse'),
               maxh=min(find(ttmse(:,:,:,rridx)==...
                             min(min(min(ttmse(:,:,:,rridx))))));
               
               [xx,tt0,chsign]=...
                   ind2sub([spacecount,diff(params.maxlag)+1,2],maxh);
            elseif strcmp(params.costfun, 'mi'),
               % find max mutual info
               maxh=min(find(ttmse(:,:,:,rridx)==...
                             max(max(max(ttmse(:,:,:,rridx))))));
               
               [xx,tt0,chsign]=...
                   ind2sub([spacecount,diff(params.maxlag)+1,2],maxh);
               
               % since this loop can take a while, record a tick.
               dbsetqueue;
               
            elseif strcmp(params.costfun, 'cc'),
               maxh=min(find(abs(ttmse(:,:,:,rridx))==...
                             max(max(abs(ttmse(:,:,1,rridx))))));
               [xx,tt0]=ind2sub([spacecount,diff(params.maxlag)+1],maxh);
               if (ttmse(xx,tt0,1,rridx)>0),
                  chsign=1;
               else
                  chsign=2;
               end
            end
            tt=tt0+params.maxlag(1)-1;
            
            %fprintf('[xx,tt,rridx]=[%d,%d,%d]\n',xx,tt,rridx);
            
            % take a baby step in the appropriate direction
            if chsign==1,
               coeff=coeffsize(xx);
            else
               coeff=-coeffsize(xx);
            end
            
            % if that's bigger than regression result, take regression
            % result???
            
            % predict response for subtraction, skip cross-val segments
            tridx=[1:(rstartidx(rridx)-1) (rendidx(rridx)+1):length(tresp)]';
            %tridx=rstartidx(rridx):rendidx(rridx);
            tridx=tridx(find(tridx-tt>0 & tridx-tt<=length(tresp) & ...
                             ~isnan(tresp(tridx,rridx))));
            
            rpred=ones(size(tresp,1),1).*nan;
            rpred(tridx)=coeff.*expstim(tridx-tt,xx);
            
            newerror=std(tr(tridx)-rpred(tridx)).^2;
            
            gg=ceil(xx./spacecount);
            tSR=zeros(spacecount,diff(params.maxlag)+1);
            tSR(mod(xx-1,spacecount)+1,tt0)=coeff;
            
            if ~isempty(basisfilter),
               H(:,:,sfsidx,respidx,rridx)=H(:,:,sfsidx,respidx,rridx)+...
                   rconv2(tSR,basisfilter{gg});
               
            elseif sum(gset{gg}{1})==0,
               H(xx,tt0,sfsidx,respidx,rridx)=...
                   H(xx,tt0,sfsidx,respidx,rridx)+coeff;
            else
               H(:,:,sfsidx,respidx,rridx)=...
                   H(:,:,sfsidx,respidx,rridx) + ...
                   gsmooth(tSR,gset{gg}{:},0);
            end
            
            ff=find(~isnan(rpred) & ~isnan(tresp(:,rridx)));
            predsofar(ff,rridx)=predsofar(ff,rridx)+rpred(ff);
            
            if strcmp(params.costfun, 'mi'),
               % a couple little tweaks.
               % first, find optimal gain
               tcc=predsofar(ff,rridx)'*resp(ff);
               if tcc<0,
                  predsofar(ff,rridx)=-predsofar(ff,rridx);
                  H(:,:,sfsidx,respidx,rridx)=-H(:,:,sfsidx,respidx,rridx);
                  fprintf('flipped sign - ');
               end
               
               [dcgparms]=fitdcgain(predsofar(ff,rridx),resp(ff,respidx));
               H(:,:,sfsidx,respidx,rridx)=...
                   dcgparms(2).*H(:,:,sfsidx,respidx,rridx);
               predsofar(ff,rridx)=dcgparms(2).*predsofar(ff,rridx);
               fprintf('scaling by %.3f:\n',dcgparms(2));
               
               newerror=std(resp(tridx,respidx)-predsofar(tridx,rridx)).^2;
               
            elseif params.activeoutnl,
               tstrf.h=sum(H(:,:,:,respidx,rridx),3);
               tstrf.zerobin=-params.maxlag(1)+1;
               tstrf.architecture='linear';
               
               linpred=strfpredict(tstrf,stim0);
               
               %traw=find_raw_nl(linpred(ff),resp(ff),0);
               %tmod_psth=raw_nl(t0,linpred(ff));
               [t0,res]=findthreshsat(linpred(ff),resp(ff),0,10,1);
               tmod_psth=threshsat(t0,linpred(ff));
               
               predsofar(ff,rridx)=tmod_psth;
               tresp(ff,rridx)=resp(ff)-tmod_psth;
               if ~onseil && rridx==1,
                  sfigure(1);
                  clf
                  subplot(2,1,1);
                  imagesc(tstrf.h);
                  axis xy;
                  if length(t0)>1,
                     title(sprintf('lminmax=[%.2f %.2f] thresh=[%.2f %.2f]',...
                                   min(linpred(ff)), max(linpred(ff)),t0));
                  else
                     title(sprintf('lminmax=[%.2f %.2f] thresh=[%.2f]',...
                                   min(linpred(ff)), max(linpred(ff)),t0));
                  end
                  
                  traw=find_raw_nl(linpred(ff),resp(ff),0);
                  subplot(2,1,2);
                  plot(traw{1},traw{2});
                  drawnow;
                  %keyboard
                end
                
                %newerror=std(resp(tridx,respidx)-predsofar(tridx,rridx)).^2;
               
            elseif 1
               
               % subtract off prediction to produce residual for next round of
               % fitting
               tresp(ff,rridx)=tresp(ff,rridx)-rpred(ff);
               
            else
               [dcgparms]=fitdcgain(predsofar(ff,rridx),resp(ff,respidx));
               H(:,:,sfsidx,respidx,rridx)=...
                   dcgparms(2).*H(:,:,sfsidx,respidx,rridx);
               predsofar(ff,rridx)=dcgparms(2).*predsofar(ff,rridx);
               fprintf('scaling by %.3f:\n',dcgparms(2));
               
               % subtract off prediction to produce residual for next round of
               % fitting
               tresp(ff,rridx)=resp(ff,respidx)-predsofar(ff,rridx);
            end
            
            
            fprintf('sfs=%d: xx=%d tt=%d step=%0.3f mse=%.4f -> %.4f\n',...
                    sfsidx,xx,tt,coeff,lasterror(rridx),newerror);
            if newerror>lasterror(rridx),
               disp('backwards??');
               %keyboard;
            end
            
            lasterror(rridx)=newerror;
         end
      end
   end
end

H=cumsum(H,3);

% don't do this since now reping for all respcounts
%H=permute(H,[1 2 3 5 4]);

%figure(gcf);
%imagesc(sum(H(:,:,:),3));

mH=mean(H,5);
eH=std(H,1,5) .* sqrt(params.resampcount-1);
lambda=1:params.sfscount;

% take mean stim/response across resamples
ntot=rsize(1);

% create to make compatible with other programs.
sSA2=zeros(spacecount,spacecount);

clear predsofar normexpstim ff expstim stim0 tresp validx trg tridx




return

%
% old format, each jackknife picks channels to include independently
% REPLACED 2006-04-23
%

for rridx=1:params.resampcount,
   tr=resp;
   %tr=stim*[0 0 0 0.5 0 0 0 0 0 2 zeros(1,20)]';
   tr(rstartidx(rridx):rendidx(rridx),:)=nan;
   
   for sfsidx=1:params.sfscount,
      
      rgoodidx=find(~isnan(tr));
      tSR=zeros(spacecount,diff(params.maxlag)+1);
      for tt=params.maxlag(1):params.maxlag(2),
         %fprintf('.');
         trg=rgoodidx;
         trg=trg(find(trg-tt>0 & trg-tt<size(stim,1)));
         
         d1=(sum(stim(trg-tt,:).^2,1));
         tSR(:,tt-params.maxlag(1)+1)=(tr(trg,:)'*stim(trg-tt,:)./d1)';
      end
      
      ttSR=abs(tSR).*repmat(sS',[1 diff(params.maxlag)+1]);
      maxh=find(ttSR(:)==max(ttSR(:)));
      
      [xx,tt0]=ind2sub(size(tSR),maxh);
      tt=tt0+params.maxlag(1)-1;
      
      rpred=tSR(maxh).*stim(:,xx);
      tr1=tr;
      if tt>=0,
         tr((1+tt):end)=tr((1+tt):end)-rpred(1:(end-tt));
      else
         tr(1:(end+tt))=tr(1:(end+tt))-rpred((-tt+1):end);
      end
      fprintf('rr=%d sfs=%d: xx=%d tt=%d\n',rridx,sfsidx,xx,tt);
      fprintf('start: %.4f pred: %.4f res: %.4f\n',...
              sqrt([nanmean(tr1.^2) nanmean(rpred.^2) nanmean(tr.^2)]));
      H(xx,tt0,sfsidx,rridx)=tSR(xx,tt0);
   end
end
H=cumsum(H,3);

H=permute(H,[1 2 3 5 4]);
%figure(gcf);
%imagesc(sum(H(:,:,:),3));

mH=mean(H,5);
eH=std(H,1,5) .* sqrt(params.resampcount-1);
lambda=1:params.sfscount;

% take mean stim/response across resamples
ntot=rsize(1);


return

irS=squeeze(std(SR,1,1));
irH=squeeze(std(H,1,1));

allt=std(reshape(permute(H,[1 2 4 5 3]),spacecount*(diff(params.maxlag)+1)* ...
                 respcount*params.resampcount,params.sfscount));
alls=std(reshape(permute(H,[1 3 4 5 2]),spacecount*params.sfscount* ...
                 respcount*params.resampcount,(diff(params.maxlag)+1)));
tH=squeeze(std(reshape(permute(H,[1 4 5 2 3]),spacecount* ...
                       respcount*params.resampcount,...
                       (diff(params.maxlag)+1),params.sfscount)));


mSR=mean(SR,4);
eSR=std(SR,1,4).* sqrt(params.resampcount-1);










