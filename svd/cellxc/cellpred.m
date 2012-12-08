% cellpred.m
%
% now... try predicting other attentional states!
%
% this might be functionalizable!
%
% inputs:
% hf - matrix containing linear filter of kernels
% threshparms - threshold info used for output nonlinearity
% mSall - mean stimulus used to calculate kernel, to be subtracted
%         off stimuli used for prediction
% params - global params with info about data files and other
%          important stuff (see cellxcnodb)
%
% returns:
% predxc - correletion coeff for each kernel and data set
% prederr - stderr of pred cc
% predp - p value of predictions, ie, whether they're significant
%         or not
% mod_psth - the predicted psth
%

spacecount=size(strf(1,1).h,1);
tbincount=size(strf(1,1).h,2);
attcount=size(strf,2);

% somewhat spaghetti-ish code allows doing either multiple batch
% preds or multiple attentional state preds
if attcount==1 & length(params.predbatch)>1,
   predcount=length(params.predbatch);
   multbatch=1;
   predokbatch=zeros(1,predcount);
   for ii=1:predcount,
      sql=['SELECT * FROM sRunData where cellid="',params.cellid,'" AND ',...
           'batch=',num2str(params.predbatch{ii})];
      rr=mysql(sql);
      if length(rr)>0,
         predokbatch(ii)=1;
      end
   end
else
   predcount=attcount;
   predokbatch=ones(1,attcount);
   multbatch=0;
end

% mod_psth: output of predictor for each kernel
mod_psth={};
act_resp={};
predxc=ones(predcount,params.nloutparm,attcount,latcount) * nan;
predp=ones(predcount,params.nloutparm,attcount,latcount);
prederr=ones(predcount,params.nloutparm,attcount,latcount) * nan;
predmse=ones(predcount,params.nloutparm,attcount,latcount) * nan;
predxc2=zeros(predcount,params.nloutparm,attcount,latcount);

for attidx=find(predokbatch),
   
   if multbatch,
      % figure out pred files for the current batch and load them
      
      fprintf('Predicting responses to batch %d.\n',params.predbatch{attidx});
      [pcellfiledata,ptimes,pbatchdata]=...
          cellfiletimes(params.cellid,params.predbatch{attidx});
      if length(pcellfiledata)==0,
         disp('error no pred file data ');
         keyboard;
      end
      
      predparams=params;
      predparams.stimfiles={};
      predparams.respfiles={};
      predparams.stimcrfs=[];
      for ii=1:length(pcellfiledata),
         predparams.stimfiles{ii}=[pcellfiledata(ii).stimpath,...
                    pcellfiledata(ii).stimfile];
         predparams.respfiles{ii}=[pcellfiledata(ii).path,...
                    pcellfiledata(ii).respfile];
         predparams.stimcrfs(ii)=pcellfiledata(ii).stimfilecrf;
      end
      
      tpredstartframe=ptimes(3).start;
      tpredstopframe=ptimes(3).stop;
      tpredfile=ptimes(3).fileidx;
      
      predparams.stimloadcmd=pbatchdata.stimloadcmd;
      predparams.stimloadparms=strsep(pbatchdata.stimloadparms,',');
      predparams.stimfiltercmd=pbatchdata.stimfiltercmd;
      predparams.stimfilterparms=strsep(pbatchdata.stimfilterparms,',');
      
      [stim,resp]=xcloadstimresp(tpredfile,tpredstartframe,...
                                 tpredstopframe,predparams);
   else
      % load stim and response segments to pred. although this will
      % generally be the same file for all attentional states, it
      % theoretically could be different. thus reload the stimulus and
      % response for each attentional state
      
      fprintf('Predicting responses to attentional state %d.\n',attidx);
      
      tpredstartframe=predstartframe(attidx);
      tpredstopframe=predstopframe(attidx);
      tpredfile=predfile(attidx);
      [stim,resp]=xcloadstimresp(tpredfile,tpredstartframe,...
                                 tpredstopframe,params);
   end
   
   predlen=tpredstopframe-tpredstartframe+1;
   
   % figure out actual valid range of response being predicted
   if multbatch,
      tr=resp(:,:,1);
   else
      tr=resp(:,:,attidx);
   end
   act_resp{attidx}=tr;
   
   % smooth actual psth if predsmoothsigma > 0
   if params.predsmoothsigma > 0 & params.respfmtcode==0,
      %fprintf('Smoothing actual response with optimal filter...\n');
      %pfilt=[1/9 2/9 1/3 2/9 1/9]';
      
      fprintf('Smoothing pred response sigma=%.1f bins...\n',...
              params.predsmoothsigma);
      tt=(-10:10)';
      pfilt=exp(-tt.^2/(2*params.predsmoothsigma.^2))./...
            (sqrt(2*pi)*params.predsmoothsigma);
      pfilt=pfilt(find(pfilt > 0.1*max(pfilt)));
      fr=tr;
      fr(isnan(fr))=0;
      fr=conv2(fr,pfilt,'same');
      tr(~isnan(tr))=fr(~isnan(tr));
   end
   
   % normalize valid bins of actual response
   tgoodidx={};
   for r=1:latcount,
      tgoodidx{r}=find(~isnan(tr(:,r)));
      %tgoodidx{r}=tgoodidx{r}(find(tgoodidx{r}>leadbincount));
      tgoodlen(r,attidx)=length(tgoodidx{r});
   end
   
   % mod_psth: output of predictor for each kernel
   mod_psth{attidx}=zeros(length(tgoodidx{r}),attcount,...
                          latcount,nlcount*params.nloutparm);
   
   % loop through each attentional state, sfs cutoff and sigma fudge
   % factor and predict the response to the fit data
   for kidx=1:attcount,
      for nlidx=1:params.nloutparm,
         fprintf('Predicting (kidx=%d att=%d nl=%d)...\n',...
                 kidx,attidx,nlidx);
         
         % subtract appropriate mean
         
         if params.meansub,
            tstim=stim'-repmat(strf(nlidx,kidx).mS,[1 size(stim,1)]);
         elseif 0,
            %mlocal=mean(stim,2);
            %tstim=stim-repmat(mlocal,[1 size(stim,2)]);
         else
            tstim=stim';
         end
         
         %nlstr={'none','post-sum thresh','pre-sum thresh'};
         if params.respfmtcode==1,
            % don't sum over time dimension... this returns a
            % time X 1 X latency matrix
            tmod_psth=kernpredict(strf(nlidx,kidx).h,tstim,1,0,0);
            if nlidx==2,
               tmod_psth(tmod_psth<0)=0;
            end
            % option to rectify... doesn't help at the moment (4/28/02)
            %fprintf(' R');
            %tmod_psth(find(tmod_psth<0))=0;
            tmod_psth=squeeze(tmod_psth);
            
         elseif strcmp(strf(nlidx,kidx).nltype,'none'),
            seplinpred=kernpredict(strf(nlidx,kidx).h,tstim,spacecount,0);
            linpred=sum(seplinpred,2);
            tmod_psth=linpred;
         elseif strcmp(strf(nlidx,kidx).nltype,'post-sum thresh'),
            linpred=kernpredict(strf(nlidx,kidx).h,tstim,1,0);
            tmod_psth=thresh(strf(nlidx,kidx).nlparms,linpred);

         elseif strcmp(strf(nlidx,kidx).nltype,'pre-sum thresh'),
            tmod_psth=kernpredict(strf(nlidx,kidx).h,tstim,spacecount,0);
            tmod_psth=thresh(strf(nlidx,kidx).nlparms,tmod_psth);
            
            
            
         elseif nlidx==1,    % no rectification
            seplinpred=kernpredict(strf(nlidx,kidx).h,tstim,spacecount,0);
            linpred=sum(seplinpred,2);
            tmod_psth=linpred;
         %elseif nlidx==2,    % rectify summed output
         %   tmod_psth=kernpredict(strf(nlidx,kidx).h,tstim,1,1);
         
         elseif nlidx==3,    % rectify each spatial channel individually
            tmod_psth=kernpredict(strf(nlidx,kidx).h,tstim,spacecount,1);
            
         elseif ceil(nlidx/nlcount)==2,
            nloutidx=ceil(nlidx/nlcount);
            nlinidx=mod(nlidx-1,nlcount)+1;
            linpred=kernpredict(hf(:,:,nlinidx,kidx),tstim,1,0);
            tmod_psth=thresh(threshfit(nloutidx,1,nlinidx,kidx),linpred);
            
         elseif ceil(nlidx/nlcount)==3,
            nloutidx=ceil(nlidx/nlcount);
            nlinidx=mod(nlidx-1,nlcount)+1;
            tmod_psth=kernpredict(hf(:,:,nlinidx,kidx),tstim,spacecount,0);
            tmod_psth=thresh(threshfit(nloutidx,1,nlinidx,kidx),tmod_psth);
            
         elseif nlidx==4,
            % don't just rectify at 0... fit a threshold
            % instead. how does this compare to nlidx=2?
            tmod_psth=linpred;
            tmod_psth=thresh(threshfit(ridx,nlidx,kidx),tmod_psth);
         elseif nlidx==5,
            tmod_psth=seplinpred;
            tmod_psth=thresh(threshfit(ridx,nlidx,kidx),tmod_psth);
         end
         
         % cycle through each latency ... for each one, compare
         % predicted response to actual response
         for ridx=1:latcount,
            
            pp=tmod_psth(tgoodidx{ridx},ridx);
            rr=tr(tgoodidx{ridx},ridx);
            mod_psth{attidx}(:,kidx,ridx,nlidx)=pp;
            
            rok=find(~isnan(pp) & ~isnan(rr));
            pp=pp(rok);
            rr=rr(rok);
            
            % compute correlation between pred and actual
            % response. only look at max lag of 0
            [cxy,exy,tt,p]=randxcov(rr,pp,0,200);
            predxc(attidx,nlidx,kidx,ridx)=cxy(1);
            prederr(attidx,nlidx,kidx,ridx)=exy(1);
            predp(attidx,nlidx,kidx,ridx)=p(1);
            
            predmse(attidx,nlidx,kidx,ridx)=var(pp-rr)./var(rr);
            
            %predxc(attidx,nlidx,kidx,ridx)=...
            %       xcov(tr(tgoodidx{ridx},ridx),...
            %            tmod_psth(tgoodidx{ridx},ridx),0,'coeff');
            
            % new crap to do hinge fit
            if multbatch,
               htr=resp(:,ridx,1);
            else
               htr=resp(:,ridx,attidx);
            end
            if sum(abs(tmod_psth(tgoodidx{ridx},ridx)))>0 & ...
                  mean(htr(tgoodidx{ridx}))>0 & exist('hingeparms','var'),
               tm=hinge(hingeparms(:,ridx,kidx),tmod_psth(:,ridx));
               
               predxc2(attidx,nlidx,kidx,ridx)=1-mean( ...
                  abs(tm(tgoodidx{ridx})-htr(tgoodidx{ridx})))...
                   ./ mean(htr(tgoodidx{ridx}));
               
               %keyboard
            else
               %disp('skip');
               predxc2(attidx,nlidx,kidx,ridx)=-1;
            end
            
            if 0 & nlcount>1,
               % fit output nl (fit sigmoid), renormalize, calc score
               tmod_psth_rec=...
                   sigmoid(sigmoidparms(:,sfsfit(ridx,nlidx,kidx),...
                                        sigfit(ridx,nlidx,kidx),ridx,kidx),...
                           tmod_psth(tgoodidx{ridx},ridx));
               tmod_psth_rec=tmod_psth_rec-mean(tmod_psth_rec);
               if sum(abs(tmod_psth_rec))>0,
                  tmod_psth_rec=tmod_psth_rec ./ std(tmod_psth_rec,1);
                  
                  mod_psth{attidx}(tgoodidx{ridx},ridx,2,kidx)=...
                      tmod_psth_rec;
                  % calculate corr coeff for NL prediction
                  predxc(attidx,2,kidx,ridx)=tr(tgoodidx{ridx},ridx)' * ...
                      tmod_psth_rec./length(tmod_psth_rec);
               end
            end
         end
      end
      fprintf('\n');
   end % for kidx
   
   if not(isempty(BATQUEUEID)),
      % record latest (sfsidx,kidx) predicted
      dbsetqueue(BATQUEUEID,(attidx-1)*100+2000);
   end
end % for attidx

clear r kidx tmod_psth attidx nlidx tpredstartframe