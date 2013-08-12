% mcdcore.m
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
%test run:
% cellxcmaster('min020a-b1',179);


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

mS=nanmean(stim(~isnan(resp),:));
sS=nanstd(stim);
sS(mS==1)=1;

if 0,
   % don't subtract mean from stim.  does this allow for fairer
   % comparison between stimulus classes (ie, torc vs speech)?
   mS=mS.*0;
   sS=nanmean(stim);
end

mR=nanmean(resp);
mSall=repmat(mS',[1 respcount]);
mRall=mR;
stim0=stim-repmat(mS,[size(stim,1),1]);
resp=resp-repmat(mR,[size(resp,1) 1]);

% KLUDGE: zero out channels that are just way too weak??
if isfield(params,'runclassid') & ismember(params.runclassid,[1 42 58]),
   ss=std(stim0);
   ss=ss./median(ss);
   ff=find(ss<0.3 & mS~=1);
   if ~isempty(ff),
      disp('zeroing out weak channels!!!');
      ff
   else
      disp('no weak channels to zero out');
   end
   stim0(:,ff)=0;
end

dimcount=getparm(params,'dimcount',1);
if dimcount<=2,
   H=zeros(spacecount,diff(params.maxlag)+1,params.sfscount,...
           dimcount,params.resampcount);
   params.architecture='linear';
else
   H=zeros(spacecount+diff(params.maxlag)+1,diff(params.maxlag)+1,...
           params.sfscount,1,params.resampcount);
   params.architecture='order2';
end

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

params.stepfrac=getparm(params,'stepfrac',30);
fprintf('stepfrac=%d\n',params.stepfrac);

respidx=1;
coeffsize=nanstd(resp(:,respidx))./sS ./ params.stepfrac;
normexpstim=expstim.*repmat(coeffsize,[length(expstim) 1]);
%absnormexpstim=abs(normexpstim);
absnormexpstim=normexpstim;
absnormexpstim(absnormexpstim<0)=0;

tresp=repmat(resp(:,respidx),[1 params.resampcount]);
lin_predsofar=zeros(size(tresp));
fprintf('\nrespidx=%d\n',respidx);
stimlen=size(stim,1);

if dimcount>2,
   % compute 2nd order stim for different time lags.
   shiftstim=zeros(length(normexpstim),diff(params.maxlag)+1);
   normnorm=std(normexpstim).*2;
   for os=0:diff(params.maxlag),
      % t2 is always >= t1
      shiftstim((1+os):end,os+1)=normexpstim(1:(end-os),:).*...
          normexpstim((1+os):end,:)./normnorm;
   end
end

lasterror=repmat(nanstd(abs(tresp(:))).^2,[params.resampcount,1]).*100;
for sfsidx=1:params.sfscount,

    %keyboard
    
    % in order to get denser sampling without wasting a lot of time,
    % do multiple iterations per strf that will actually be tested
    % at the next stage
    
    for subidx=1:iterpersfscount,
       
        % can't use quick & dirty calculations because preds for each
        % resample create different residuals for each resample on
        % subsequent iterations!
        fprintf('iteration %d/%d:\n',sfsidx,subidx);
        
        if dimcount==2,
           ttmse=zeros(expspacecount,diff(params.maxlag)+1,2*dimcount,...
                       params.resampcount);
        elseif dimcount>2,
           ttmse=zeros(expspacecount+diff(params.maxlag)+1,...
                       diff(params.maxlag)+1,2,...
                       params.resampcount);
        end
        
        for rridx=1:params.resampcount,
            tr=tresp(:,rridx);
            tr(rstartidx(rridx):rendidx(rridx),:)=nan;
            rgoodidx=find(~isnan(tr));
            

            % go through each possible time lag
            for tt=params.maxlag(1):params.maxlag(2),
                trg=rgoodidx;
                trg=trg(find(trg-tt>0 & trg-tt<size(normexpstim,1)));
                for xx=1:expspacecount,
                    % test increment of each channel
                    ttmse(xx,tt-params.maxlag(1)+1,1,rridx)=...
                        var(tr(rgoodidx)-normexpstim(trg-tt,xx));
                    % test decrement of each channel
                    ttmse(xx,tt-params.maxlag(1)+1,2,rridx)=...
                        var(tr(rgoodidx)+normexpstim(trg-tt,xx));
                    
                    if dimcount==2,
                       % test increment of each channel
                       ttmse(xx,tt-params.maxlag(1)+1,3,rridx)=...
                         mean(abs(tr(rgoodidx)-absnormexpstim(trg-tt,xx)).^2);
                       % test decrement of each channel
                       ttmse(xx,tt-params.maxlag(1)+1,4,rridx)=...
                         mean(abs(tr(rgoodidx)+absnormexpstim(trg-tt,xx)).^2);
                    end
                end
                if dimcount>2,
                   for tt2=tt:params.maxlag(2),
                      os=tt2-tt;
                      
                      trg2=trg(find(trg-tt2>0));
                      tstim=shiftstim(trg2-tt,os+1);
                      
                      xxidx=2+tt2-params.maxlag(1);
                      % test increment of each channel
                      ttmse(xxidx,tt-params.maxlag(1)+1,1,rridx)=...
                          var(tr(trg2)-tstim);
                      % test decrement of each channel
                      ttmse(xxidx,tt-params.maxlag(1)+1,2,rridx)=...
                          var(tr(trg2)+tstim);
                   end
                end
                
            end
            
            % and compute prediction for this new channel
            rpred=ones(size(tresp)).*nan;
            ttmse(ttmse==0)=max(max(max(ttmse(:,:,:,rridx))));
            maxh=min(find(ttmse(:,:,:,rridx)==...
                min(min(min(ttmse(:,:,:,rridx))))));
            
            [xx,tt0,chsign,dim]=...
                ind2sub([size(ttmse,1),diff(params.maxlag)+1,2,2],maxh);
            tt=tt0+params.maxlag(1)-1;
            
            %fprintf('[xx,tt,rridx]=[%d,%d,%d]\n',xx,tt,rridx);
            
            % take a baby step in the appropriate direction
            coeff=-(chsign*2-3) * coeffsize;
            
            % predict response for subtraction, skip cross-val segments
            tridx=[1:(rstartidx(rridx)-1) (rendidx(rridx)+1):length(tresp)]';
            %tridx=rstartidx(rridx):rendidx(rridx);
            tridx=tridx(find(tridx-tt>0 & tridx-tt<=length(tresp) & ...
                ~isnan(tresp(tridx,rridx))));
            
            rpred=ones(size(tresp,1),1).*nan;
            if xx>spacecount,
               tt2=xx-spacecount+params.maxlag(1)-1;
               tridx=tridx(find(tridx-tt2>0 & tridx-tt2<=length(tresp)));
               coeff=coeff./normnorm(1).*coeffsize;
               rpred(tridx)=coeff.*expstim(tridx-tt,1).*expstim(tridx-tt2,1);
               %keyboard
            elseif dim==1,
               rpred(tridx)=coeff.*expstim(tridx-tt,xx);
            else
               rpred(tridx)=expstim(tridx-tt,xx);
               rpred(tridx(rpred(tridx)<0))=0;
               rpred(tridx)=coeff.*rpred(tridx);
            end
            newerror=std(tr(tridx)-rpred(tridx)).^2;
            
            gg=ceil(xx./spacecount);
            
            tSR=zeros(spacecount,diff(params.maxlag)+1);
            tSR(mod(xx-1,spacecount)+1,tt0)=coeff;
            if dimcount>2
               if xx<=spacecount,
                  H(xx,tt0,sfsidx,respidx,rridx)=...
                      H(xx,tt0,sfsidx,respidx,rridx)+coeff;
               else
                  % split second-order
                  H(xx,tt0,sfsidx,respidx,rridx)=...
                      H(xx,tt0,sfsidx,respidx,rridx)+coeff./2;
                  H(tt0+spacecount,xx-spacecount,sfsidx,respidx,rridx)=...
                      H(tt0+spacecount,xx-spacecount,sfsidx, ...
                        respidx,rridx)+coeff./2;
               end
            elseif ~isempty(basisfilter),
                H(:,:,sfsidx,dim,rridx)=H(:,:,sfsidx,dim,rridx)+...
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
            lin_predsofar(ff,rridx)=lin_predsofar(ff,rridx)+rpred(ff);

            
            % subtract off prediction to produce residual for next round of
            % fitting
            tresp(ff,rridx)=tresp(ff,rridx)-rpred(ff);
            fprintf('sfs=%d: xx=%d tt=%d stp=%0.3f dim=%d mse=%.4f -> %.4f\n',...
                    sfsidx,xx,tt,coeff,(xx>spacecount)+1,...
                    lasterror(rridx),newerror);
            if newerror>lasterror(rridx),
                disp('backwards??');
                %keyboard;
            end
            %if rridx==1,
            %   keyboard
            %end

            lasterror(rridx)=newerror;
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

clear shiftstim stim0 tresp expstim normexpstim absnormexpstim 
clear lin_predsofar ttmse rpred tridx trg trg2 


return


% code below shoul be transferred/morphed into xcfit2

keyboard

figure;
imagesc(mean(H(1:end,:,:,:,12),3));
colorbar

% debug 2nd order prediction code
shiftstim=zeros(size(stim0,1),diff(params.maxlag)+1);
for os=0:diff(params.maxlag),
   % t2 is always >= t1
   shiftstim((1+os):end,os+1)=stim0(1:(end-os),:).*...
       stim0((1+os):end,:);
end
pp=zeros(size(resp));
for t1=params.maxlag(1):params.maxlag(2),
   gidx=max(t1+1,1):min(t1+stimlen,stimlen);
   pp(gidx)=pp(gidx)+stim(gidx-t1,:)*H(1,t1-params.maxlag(1)+1,end);
   
   for t2=t1:params.maxlag(2),
      os=t2-t1;
      if os>0
         pp(gidx)=pp(gidx)+shiftstim(gidx-t1,os+1)* 2 *...
            H(t2-params.maxlag(1)+2,t1-params.maxlag(1)+1,end);
      else
         pp(gidx)=pp(gidx)+shiftstim(gidx-t1,os+1)*...
            H(t2-params.maxlag(1)+2,t1-params.maxlag(1)+1,end);
      end
   end
end

rgoodidx=find(~isnan(resp));
xcov(pp(rgoodidx),resp(rgoodidx),0,'coeff')

if 1,
    % find optimal output NL
    pred1=kernpredict(mH(:,:,12,1),stim0',1,0,1,-params.maxlag(1)+1);
    pred2=kernpredict(mH(:,:,12,2),stim0',1,0,1,-params.maxlag(1)+1);
    ff=find(~isnan(resp));
    
    bincount=10;
    pp1=zeros(bincount,1);
    pp2=zeros(bincount,1);
    rr1=zeros(bincount,1);
    rr2=zeros(bincount,1);
    rr=zeros(bincount,bincount);
    [ss,si1]=sort(pred1(ff));
    edges1=round(linspace(1,length(si1)+1,bincount+1));
    [ss,si2]=sort(pred2(ff));
    edges2=round(linspace(1,length(si2)+1,bincount+1));
    for bb1=1:bincount,
        pp1(bb1)=mean(pred1(ff(si1(edges1(bb1):(edges1(bb1+1)-1)))));
        rr1(bb1)=mean(resp(ff(si1(edges1(bb1):(edges1(bb1+1)-1)))));
        for bb2=1:bincount,
            if bb1==1,
                pp2(bb2)=mean(pred2(ff(si2(edges2(bb2):(edges2(bb2+1)-1)))));
                rr2(bb2)=mean(resp(ff(si2(edges2(bb2):(edges2(bb2+1)-1)))));
            end
            tt=intersect(si1(edges1(bb1):(edges1(bb1+1)-1)),si2(edges2(bb2):(edges2(bb2+1)-1)));
            if length(tt)>0,
                rr(bb1,bb2)=mean(resp(ff(tt)));
            end
        end
    end
    figure;
    subplot(2,2,1);
    plot(squeeze(mH(:,:,end,:)));
    legend('lin','squ');
    subplot(2,2,2);
    plot(pp1,rr1);
    subplot(2,2,3);
    plot(rr2,pp2);
    subplot(2,2,4);
    imagesc(pp1,pp2,rr');
    xlabel('pred1');
    ylabel('pred2');
    axis xy
end



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










