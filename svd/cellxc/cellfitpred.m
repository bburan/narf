% cellfitpred.m
%
% now it's time to fit the data!!!  ack!   ????

% we have two things to fit: the svd cutoff (sfsidx) and the
% shrinkage filter fudge factor (sfiltsigma)

MAXPREDLEN=10000;
spacecount=size(mH,1);
tbincount=size(mH,2);
sfscount=size(mH,3);
respcount=size(mH,4);
attcount=size(mH,5);
lagrange=(1-maxlag(1)):(tbincount-maxlag(2)+10);

sfsidxtouse=(1:sfscount)';  % test all sfs cutoffs

% define a range of sigmas for combining resampled estimates (via
% shrinkage filter)
%sigrange=[1.0 1.2 1.4 1.7 2.0 2.4 2.8 3.3 3.8 4.6 5.2 6.0 7.0 8.0
%9.0];
if respfmtcode==0,
   sigrange=exp(linspace(0,1.5,7));
else
   sigrange=exp(linspace(0,1.5,10));
   %sigrange=exp(linspace(-1,log(2.5),10));
end  
if 0, % run fast  for debugging
   sigrange=exp([0 1]);
   sfsidxtouse=0:2:10;
   sfscount=length(sfsidxtouse);
end
sigcount=length(sigrange);

% pred scores for various fit values
nlcount=2;

% for PFTH preds, separate score for each latency bin
if respfmtcode==0,
   latcount=1;
else
   latcount=tbincount;
end
xc=zeros(sfscount,sigcount,latcount,nlcount,attcount);
sigmoidparms=zeros(4,sfscount,sigcount,latcount,attcount);

for attidx=1:attcount,
   fprintf('Fitting for attentional state %d.\n',attidx);
   
   % load stim and response segments to fit. although this will
   % generally be the same file for all attentional states, it
   % theoretically could be different. thus reload the stimulus and
   % response for each attentional state
   
   % adjust fitstart frame with extra frames to beginning so that
   % first valid prediction bin lines up with the beginning of fit
   % data.
   if respfmtcode==0,
      leadbincount=min([fitstartframe(attidx)-1 tbincount-1]);
   else
      leadbincount=0;
   end
   fitstartframe(attidx)=fitstartframe(attidx)-leadbincount;
   
   stim=feval(stimloadcmd,stimfiles{fitfile(attidx)},fitstartframe(attidx),...
              fitstopframe(attidx),stimloadparms{:});
   
   % filter stimulus segment if selected
   if ~isempty(stimfiltercmd),
      stim=feval(stimfiltercmd,stim,stimfilterparms{:});
   end
   iconside=size(stim);
   movlen=iconside(end);
   iconside=iconside(1:(end-1));
   if length(size(stim))>2,
      % reshape all spatial dims into one
      stim=reshape(stim,prod(iconside),movlen);
   end
   
   resp=feval(resploadcmd,respfiles{fitfile(attidx)},resploadparms{:});
   resp=resp(fitstartframe(attidx):fitstopframe(attidx),:,:);
   
   % filter response (eg smooth, pick out ICs, PCs? currently not used)
   if ~isempty(respfiltercmd),
      resp=feval(respfiltercmd,resp,respfilterparms{:});
   end
   predlen=fitstopframe(attidx)-fitstartframe(attidx)+1;

   % smooth actual psth if movformat.predsmoothsigma > 0
   %if movformat.predsmoothsigma > 0 & ~ismember(movformat.predtype,[2,3]),
   %   tt=(-10:10)';
   %   pfilt=exp(-tt.^2/(2*movformat.predsmoothsigma.^2))./...
   %         (sqrt(2*pi)*movformat.predsmoothsigma);
   %   pfilt=pfilt(find(pfilt > 0.1*max(pfilt)));
   %   pflen=(length(pfilt)-1)/2;
   %   fresp_psth=conv(tresp_psth.*(tresp_psth>=0),pfilt);
   %   fresp_psth=fresp_psth(pflen+1:length(fresp_psth)-pflen);
   %   sresp_psth=(tresp_psth>=0).*fresp_psth+(tresp_psth<0)*(-1);
   %else
   %   sresp_psth=tresp_psth;
   %end
   
   % figure out actual valid range of response being predicted
   tr=resp(leadbincount+1:predlen,:,attidx);
   tgoodidx={};
   % normalize valid bins of actual response
   for r=1:latcount,
      tgoodidx{r}=find(~isnan(tr(:,r)));
      % tr=tr(tgoodidx,:);
      if sum(abs(tr(tgoodidx{r},r)))>0,
         tr(tgoodidx{r},r)=(tr(tgoodidx{r},r)-mean(tr(tgoodidx{r},r))) ...
             ./std(tr(tgoodidx{r},r),1);
      end
   end
   %tr=tr(:);   % reshape as a vector for quick xc DON'T DO IT FOR PFTH
   
   % loop through each attentional state, sfs cutoff and sigma fudge
   % factor and predict the response to the fit data
   for sfsidx=1:sfscount,
      for sigidx=1:sigcount,
         for r=1:respcount,
            fprintf('Testing fit (att=%d sfsidx=%d sig=%.2f r=%d)...\n',...
                    attidx,sfsidxtouse(sfsidx),sigrange(sigidx),r);
            
            % current kernel is shrinkage-filtered from current (sfs,r,attidx)
            % only use causal time lags for pred!
            smm=mH(:,(1-maxlag(1)):end,sfsidx,r,attidx);
            sms=eH(:,(1-maxlag(1)):end,sfsidx,r,attidx) .* sigrange(sigidx);
            smd=abs(smm)./(sms+(sms==0));
            smd=1-smd.^(-2);
            smd=smd.*(smd>0);
            smd(find(isnan(smd)))=0;
            hf=smm.*smd;
            
            if respfmtcode==1,
               
               % don't sum over time dimension... this returns a
               % time X 1 X latency matrix
	       tmod_psth=kernpredict(hf,stim,1,0,0);
	       tmod_psth=squeeze(tmod_psth);
               
	       % "resp" space comes back in columns. time still in
               % rows. need to adjust normalization routines to
               % deal with more than one column in tmod_psth
	    else
	       tmod_psth=kernpredict(hf,stim,spacecount,1);
	       %tmod_psth=kernpredict(hf,stim,1,0);
            end
            
            % remove leading bins that don't matter
	    tmod_psth=tmod_psth(leadbincount+1:predlen,:);
	    for ridx=1:latcount,
	       % normalize for appropriate correlation coefficient calculation
	       tmod_psth(:,ridx)=(tmod_psth(:,ridx) - ...
				  mean(tmod_psth(tgoodidx{ridx},ridx)));
	       if sum(abs(tmod_psth(tgoodidx{ridx},ridx)))>0,
		  tmod_psth(:,ridx)=tmod_psth(:,ridx) ./ ...
		      std(tmod_psth(tgoodidx{ridx},ridx),1);
		  
                  xc(sfsidx,sigidx,ridx,1,attidx)=...
                      tr(tgoodidx{ridx},ridx)' * ...
                      tmod_psth(tgoodidx{ridx},ridx)./length(tgoodidx{ridx});
		  
		  % fit output nl (fit sigmoid), renormalize, calc score
		  [outnl,xx,sigmoidparms(:,sfsidx,sigidx,ridx,attidx)]=...
		      findoutnl(tmod_psth(tgoodidx{ridx},ridx),...
                                tr(tgoodidx{ridx},ridx),0);
		  tmod_psth_rec=sigmoid(sigmoidparms(:,sfsidx,sigidx,...
                              ridx,attidx),tmod_psth(tgoodidx{ridx},ridx));
		  tmod_psth_rec=tmod_psth_rec-mean(tmod_psth_rec);
		  if std(tmod_psth_rec)>0,
		     tmod_psth_rec=tmod_psth_rec ./ std(tmod_psth_rec,1);
                     
                     % calculate corr coeff for NL prediction
                     xc(sfsidx,sigidx,ridx,2,attidx)=...
                         tr(tgoodidx{ridx},ridx)' * ...
                         tmod_psth_rec./length(tmod_psth_rec);
		  end
               else
                  disp('zero pred cond');
	       end
            end
         end
      end
      
      if not(isempty(BATQUEUEID)),
         % record latest (sfsidx,attidx) predicted
         dbsetqueue(BATQUEUEID,sfsidx+(attidx-1)*100+1000);
      end
   end % for sfsidx
end % for attidx

% display fit results for different sfs and sigsmooth values in
% each attentional state with lin and NL outputs
%
% also compute and extract the corresponding kernel for each
% attentional/linearity condition
figure(1);
clf
latuse=min([latcount 7]);
maxc=max(max(max(max(xc(:,:,latuse,:,:)))));
minc=min(min(min(min(xc(:,:,latuse,:,:)))));
sfsfit=zeros(latcount,nlcount,attcount);
sigfit=zeros(latcount,nlcount,attcount);
hf=zeros(spacecount,size(mH,2)+maxlag(1),attcount);
for attidx=1:attcount,
   for nlidx=1:nlcount,
      for latidx=1:latcount,
         txc=xc(:,:,latidx,nlidx,attidx);
         
         %smooth txc?
         if 0,
            [XX,YY]=meshgrid(-2:2,-2:2);
            sx=1.5; sy=2;
            g=exp(-XX.^2./(2*sx.^2) - YY.^2./(2*sy.^2));
            g=g./sum(g(:));
            txc=rconv2(txc,g);
         end
         maxidx=min(find(txc==max(txc(:))));
         [sfsfit(latidx,nlidx,attidx),sigfit(latidx,nlidx,attidx)]=...
             ind2sub([sfscount,sigcount],maxidx);
      
         if latidx==latuse,
            subplot(2,attcount,(nlidx-1)*attcount+attidx);
            imagesc(xc(:,:,latidx,nlidx,attidx),[minc maxc]);
            
            hold on
            plot(sigfit(latidx,nlidx,attidx),sfsfit(latidx,nlidx,attidx),'x');
            hold off
            title(sprintf('FIT nl=%d attidx=%d',nlidx,attidx));
         end
         if nlidx==1,
            % use best fit for non-attentional state
            if latcount>1,
               smm=mH(:,latidx,sfsfit(latidx,nlidx,attidx),1,attidx);
               sms=eH(:,latidx,sfsfit(latidx,nlidx,attidx),1,attidx) .* ...
                   sigrange(sigfit(latidx,nlidx,attidx));
               smd=abs(smm)./(sms+(sms==0));
               smd=1-smd.^(-2);
               smd=smd.*(smd>0);
               smd(find(isnan(smd)))=0;
               hf(:,latidx,attidx)=smm.*smd;
               
            else
               smm=mH(:,-maxlag(1)+1:end,sfsfit(latidx,nlidx,1),1,attidx);
               sms=eH(:,-maxlag(1)+1:end,sfsfit(latidx,nlidx,1),1,attidx) .* ...
                   sigrange(sigfit(latidx,nlidx,1));
               %smm=mH(:,-maxlag(1)+1:end,sfsfit(latidx,nlidx,attidx),1,attidx);
               %sms=eH(:,-maxlag(1)+1:end,sfsfit(latidx,nlidx,attidx),1,attidx) .* ...
               %    sigrange(sigfit(latidx,nlidx,attidx));
               smd=abs(smm)./(sms+(sms==0));
               smd=1-smd.^(-2);
               smd=smd.*(smd>0);
               smd(find(isnan(smd)))=0;
               hf(:,:,attidx)=smm.*smd;
            end
         end
      end 
   end   
end
colorbar
colormap(hot);

figure(2)
showkern(hf(:,:,:,1),kernfmt,iconside);
%showkern(squeeze(mSR),kernfmt,iconside)
%showkern(squeeze(mH(:,:,end,:,:)),kernfmt,iconside)

%now... try predicting other attentional states!

% mod_psth: output of predictor for each kernel
mod_psth={};
mod_psth_rec={};
predxc=zeros(attcount,attcount,latcount,nlcount);

for attidx=1:attcount,
   fprintf('Predicting responses to attentional state %d.\n',attidx);
   
   % load stim and response segments to pred. although this will
   % generally be the same file for all attentional states, it
   % theoretically could be different. thus reload the stimulus and
   % response for each attentional state
   
   % adjust predstart frame with extra frames to beginning so that
   % first valid prediction bin lines up with the beginning of pred data.
   leadbincount=min([predstartframe(attidx)-1 tbincount-1]);
   predstartframe(attidx)=predstartframe(attidx)-leadbincount;
   
   stim=feval(stimloadcmd,stimfiles{predfile(attidx)},...
              predstartframe(attidx),predstopframe(attidx),stimloadparms{:});
   
   % filter stimulus segment if selected
   if ~isempty(stimfiltercmd),
      stim=feval(stimfiltercmd,stim,stimfilterparms{:});
   end
   iconside=size(stim);
   movlen=iconside(end);
   iconside=iconside(1:(end-1));
   if length(size(stim))>2,
      % reshape all spatial dims into one
      stim=reshape(stim,prod(iconside),movlen);
   end
   
   resp=feval(resploadcmd,respfiles{predfile(attidx)},resploadparms{:});
   resp=resp(predstartframe(attidx):predstopframe(attidx),:,:);
   
   % filter response (eg smooth, pick out ICs, PCs? currently not used)
   if ~isempty(respfiltercmd),
      resp=feval(respfiltercmd,resp,respfilterparms{:});
   end
   predlen=predstopframe(attidx)-predstartframe(attidx)+1;
   
   % smooth actual psth if movformat.predsmoothsigma > 0
   %if movformat.predsmoothsigma > 0 & ~ismember(movformat.predtype,[2,3]),
   %   tt=(-10:10)';
   %   pfilt=exp(-tt.^2/(2*movformat.predsmoothsigma.^2))./...
   %         (sqrt(2*pi)*movformat.predsmoothsigma);
   %   pfilt=pfilt(find(pfilt > 0.1*max(pfilt)));
   %   pflen=(length(pfilt)-1)/2;
   %   fresp_psth=conv(tresp_psth.*(tresp_psth>=0),pfilt);
   %   fresp_psth=fresp_psth(pflen+1:length(fresp_psth)-pflen);
   %   sresp_psth=(tresp_psth>=0).*fresp_psth+(tresp_psth<0)*(-1);
   %else
   %   sresp_psth=tresp_psth;
   %end
   
   %disp('stopping to figure out how to do PFTH preds');
   %keyboard
   
   % mod_psth: output of predictor for each kernel
   mod_psth{attidx}=zeros(predlen-leadbincount,latcount,nlcount,attcount);

   % figure out actual valid range of response being predicted
   tr=resp(leadbincount+1:predlen,:,attidx);
   tgoodidx={};
   % normalize valid bins of actual response
   for r=1:latcount,
      tgoodidx{r}=find(~isnan(tr(:,r)));
      % tr=tr(tgoodidx,:);
      if sum(abs(tr(tgoodidx{r},r)))>0,
         tr(tgoodidx{r},r)=(tr(tgoodidx{r},r)-mean(tr(tgoodidx{r},r))) ...
             ./std(tr(tgoodidx{r},r),1);
      end
   end
   %tr=tr(:);   % reshape as a vector for quick xc
   
   % loop through each attentional state, sfs cutoff and sigma fudge
   % factor and predict the response to the fit data
   for kidx=1:attcount,
      for r=1:respcount,
         fprintf('Predicting (kidx=%d att=%d r=%d)...\n',kidx,attidx,r);
         
         if respfmtcode==1,
            % don't sum over time dimension... this returns a
            % time X 1 X latency matrix
            tmod_psth=kernpredict(hf(:,:,kidx,1),stim,1,0,0);
            tmod_psth=squeeze(tmod_psth);
         else
            tmod_psth=kernpredict(hf(:,:,kidx,1),stim,spacecount,1);
            %tmod_psth=kernpredict(hf(:,:,kidx,1),stim,1,0);
         end
         
         tmod_psth=tmod_psth(leadbincount+1:predlen,:);
         
         % cycle through each latency ... for each one, compare
         % predicted response to actual response
         for ridx=1:latcount,
            % normalize for appropriate correlation coefficient calculation
            tmod_psth(:,ridx)=(tmod_psth(:,ridx) - ...
                               mean(tmod_psth(tgoodidx{ridx},ridx)));
            if sum(abs(tmod_psth(:,ridx)))>0,
               tmod_psth(:,ridx)=tmod_psth(:,ridx) ./ ...
                   std(tmod_psth(tgoodidx{ridx},ridx),1);
               % save results
               mod_psth{attidx}(:,ridx,1,kidx)=tmod_psth(:,ridx);
               predxc(attidx,kidx,ridx,1)=tr(tgoodidx{ridx},ridx)' * ...
                   tmod_psth(tgoodidx{ridx},ridx)./length(tgoodidx{ridx});
               
               % fit output nl (fit sigmoid), renormalize, calc score
               tmod_psth_rec=sigmoid(sigmoidparms(:,sfsfit(ridx,nlidx,kidx),...
                    sigfit(ridx,nlidx,kidx),ridx,kidx),tmod_psth(tgoodidx{ridx},ridx));
               tmod_psth_rec=tmod_psth_rec-mean(tmod_psth_rec);
               if sum(abs(tmod_psth_rec))>0,
                  tmod_psth_rec=tmod_psth_rec ./ std(tmod_psth_rec,1);
                  
                  mod_psth{attidx}(tgoodidx{ridx},ridx,2,kidx)=...
                      tmod_psth_rec;
                  % calculate corr coeff for NL prediction
                  predxc(attidx,kidx,ridx,2)=tr(tgoodidx{ridx},ridx)' * ...
                      tmod_psth_rec./length(tmod_psth_rec);
               end
            end
         end
      end % for kidx
      
      if not(isempty(BATQUEUEID)),
         % record latest (sfsidx,kidx) predicted
         dbsetqueue(BATQUEUEID,(attidx-1)*100+2000);
      end
   end % for sfsidx
end % for attidx

if respfmtcode==1 & attcount>1,
   %predxc
   figure(3);
   clf
   minx=min(min(min(predxc(:,:,:,1))));
   maxx=max(max(max(predxc(:,:,:,1))));
   
   subplot(ceil(latcount/3)+1,2,1);
   plot(nanmean(resp(:,:,1)));
   title('mean resp');
   subplot(ceil(latcount/3)+1,2,2);
   plot(squeeze(predxc(1,1,:,1)));
   title('full xc');
   
   for ii=1:latcount,
      subplot(ceil(latcount/3)+1,3,ii+3);
      imagesc(predxc(:,:,ii),[minx maxx]);
      axis image
      axis off
      if ii==1,
         title(sprintf('cell %s',cellid));
      else
         title(sprintf('lat=%d',ii));
      end
      
      if ii==latcount,
         colorbar
      end
   end
   colormap(hot);
   
   keyboard
   
   print -f3
else
   cellid
   predxc
   keyboard
end
