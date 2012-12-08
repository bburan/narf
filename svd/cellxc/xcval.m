% function res=xcval(strf,params,data);
%
% strf - strf structure. important parts of strf:
%        strf(nn).h    space X timelag kernel
%                .mS   space X 1 mean subtracted from each stim
%                      frame before convolving with h
%                .nltype  name of output NL ('none' means nothing)
%                .nlparms parameters for output NL
%                         called rout=feval(nltype,nlparms,rin)
%                         where rin = h * (data.stim - mS)
%
% data - data.stim    time X space stimulus
%        data.resp    time X 1 response recorded simulutaneously
%        OR
%        parameters explaining how to load and preprocess stimdata
%
function res=xcval(strf,params,data);

global ESTIMATIONPHASE VALIDATIONPHASE

if params.nrandxcov>1,
   % ie, this is gonna take a while, report that xcval is executing
   verbosity=1;
else
   verbosity=0;
end

if verbosity,
   disp('xcval.m:');
end

global BATQUEUEID

strfcount=size(strf,1);
attcount=size(strf,2);
spacecount=size(strf(1,1).h,1);
tbincount=size(strf(1,1).h,2);

if params.respfmtcode==0,
   latcount=1;
else
   latcount=tbincount;
end

% mod_psth: output of predictor for each kernel
mod_psth={};
act_resp={};
predxc=ones(strfcount,attcount,latcount) * nan;
predinf=ones(strfcount,attcount,latcount) * nan;
predone=ones(strfcount,attcount,latcount) * nan;
predfix=ones(strfcount,attcount,latcount) * nan;
predfixerr=ones(strfcount,attcount,latcount) * nan;
predp=ones(strfcount,attcount,latcount);
prederr=ones(strfcount,attcount,latcount) * nan;
predsignal=ones(1,attcount,latcount);
predmse=ones(strfcount,attcount,latcount) * nan;
predxc2=zeros(strfcount,attcount,latcount);

if exist('data','var') & isfield(data,'stim')
   
   if size(data.stim,1)==0,
      res.predxc=predxc;
      res.predone=predone;
      res.predinf=predinf;
      res.predp=predp;
      res.prederr=prederr;
      res.predsignal=predsignal;
      res.predmse=predmse;
      res.predxc2=predxc2;
      res.mod_psth=mod_psth;
      res.act_resp=act_resp;
      res.predfix=predfix;
      res.predfixerr=predfixerr;
      return
   end
   
   % data.stim and data.resp already loaded
   tpredstartframe=1;
   tpredstopframe=size(data.stim,1);
   tpredfile=1;
else
   % figure out what stim/resp data to load
   if exist('data','var'),
      times=data;
   else
      times=params.times(3);
   end
   
   tpredstartframe=times(1).start;
   tpredstopframe=times(1).stop;
   tpredfile=times(1).fileidx;
   
   % load the stim/resp data
   [data.stim,data.resp]=xcloadstimresp(tpredfile,tpredstartframe,...
                                        tpredstopframe,params);
end

tresp=[];
if VALIDATIONPHASE && isfield(params,'times'),
   times=params.times(3);
   
   tpredstartframe=times.start;
   tpredstopframe=times.stop;
   tpredfile=times.fileidx;
   if strcmp(params.resploadcmd,'respload'),
      tresp=respload(params.respfiles{tpredfile(1)},'r',1,1,0);
   elseif strcmp(params.resploadcmd,'loadspikeraster'),
      if isstruct(params.resploadparms{1})
         params.resploadparms{1}.psthonly=0;
         tresp=feval(params.resploadcmd,params.respfiles{tpredfile(1)},...
                     params.resploadparms{:});
      else
         tresp=loadspikeraster(params.respfiles{tpredfile(1)},...
                      params.resploadparms{1:4},{'Reference'},0);
      end
   else
      tresp=feval(params.resploadcmd,params.respfiles{tpredfile(1)},...
                  params.resploadparms{:});
      %tresp=feval(params.resploadcmd,params.respfiles{tpredfile(1)},...
      %            10,0);
   end
   if size(tresp,1)==1,
      tresp=tresp';
   end
   if size(tresp,2)>1,
      if ~strcmp(params.resploadcmd,'loadspikeraster'),
         tresp=tresp(tpredstartframe:tpredstopframe,2:end);
         tresp=compact_raster_matrix3(tresp);
         firstwithnans=min([find(sum(isnan(tresp))>0) size(tresp,2)+1]);
         
         if firstwithnans>2,
            tresp=tresp(:,1:firstwithnans-1);
         end
      else
         tresp=tresp(tpredstartframe:tpredstopframe,:);
      end
   else
      tresp=tresp(tpredstartframe:tpredstopframe);
   end
   
   repcount=size(tresp,2);
   fprintf('(%d reps): ceiling... ',repcount);
   if 0 && repcount>10,
       fprintf('(%d reps): ceiling direct measure... ',repcount);
       
   elseif ~strcmp(params.resploadcmd,'load_ecog') && ...
         ~strcmp(params.resploadcmd,'loadgammaraster') && ...
         ~strcmp(params.resploadcmd,'loadevpraster') && ...
         length(find(tresp(:,1)<0))==0,
       fprintf('(%d reps): ceiling gamma model... ',repcount);
      gg=find(~isnan(tresp(:,1)));
      [mu,alpha,beta]=reversepoisson(tresp(gg,1));
      rmax=singletrialceiling(tresp(gg,1),alpha,beta);
   else
      rmax=0;
   end
end

predlen=size(data.stim,1);

% each attentional condition might give a different response to
% measure against. er... diff stim too, eh?
for attidx=1:attcount,
   tr=data.resp(:,attidx);
   %traster=data.raster;
   % smooth actual psth if predsmoothsigma > 0
   if VALIDATIONPHASE & params.predsmoothsigma > 0 & params.respfmtcode==0,
      %fprintf('Smoothing actual response with optimal filter...\n');
      %pfilt=[1/9 2/9 1/3 2/9 1/9]';
      
      if verbosity,
         fprintf('(smsig=%.1f)',params.predsmoothsigma);
      end
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
   %mod_psth{attidx}=zeros(length(tgoodidx{r}),latcount,strfcount);
   mod_psth{attidx}=zeros(size(tr,1),latcount,strfcount);
   
   % loop through each attentional state, sfs cutoff and sigma fudge
   % factor and predict the response to the fit data
   %fprintf('Predicting (strfidx=%d attidx=%d)...\n',...
   %        strfidx,attidx);
   
   for strfidx=1:strfcount,      
      % subtract appropriate mean
      if params.meansub,
         if strfidx==1 | sum(abs(strf(strfidx,attidx).mS- ...
                                 strf(strfidx-1,attidx).mS))>0,
            tstim=data.stim'-repmat(strf(strfidx,attidx).mS,...
                                    [1 size(data.stim,1)]);
         end
      else
         tstim=data.stim';
      end

      %nlstr={'none','post-sum thresh','exp thresh','pre-sum thresh'};
      if ~isfield(strf(strfidx,attidx),'zerobin'),
         strf(strfidx,attidx).zerobin=1;
      end
      
      if params.respfmtcode==1,
         
         % don't sum over time dimension... this returns a
         % time X 1 X latency matrix
         tmod_psth=strfpredict(strf(strfidx,attidx),tstim');
         
         %if nlidx==2,
         %   tmod_psth(tmod_psth<0)=0;
         %end
         % option to rectify... doesn't help at the moment (4/28/02)
         %fprintf(' R');
         %tmod_psth(find(tmod_psth<0))=0;
         tmod_psth=squeeze(tmod_psth);
         
      else
         [tmod_psth,linpred]=strfpredict(strf(strfidx,attidx),tstim');
      end
      
      act_resp{attidx}=tr;
      if size(tmod_psth,2)>size(tr,2) && isfield(strf.parms,'maxlag'),
         % add delay line to response
         plagcount=strf.parms.maxlag(2)-strf.parms.maxlag(1)+1;
         rlen=length(tr);
         new_tr=zeros(rlen,plagcount);
         for ridx=1:plagcount,
            os=ridx-strf.parms.maxlag(2)-1;
            if os<0
               new_tr(1-os:rlen,ridx)=tr(1:rlen+os);
            else
               new_tr(1:rlen-os,ridx)=tr(os+1:rlen);
            end
            tgoodidx{ridx}=tgoodidx{1};
         end
         tr=new_tr;
         latcount=plagcount;
      end
      
      % cycle through each latency ... for each one, compare
      % predicted response to actual response
      for ridx=1:latcount,
         pp=tmod_psth(tgoodidx{ridx},ridx);
         rr=tr(tgoodidx{ridx},ridx);
         
         mod_psth{attidx}(:,ridx,strfidx)=tmod_psth(:,ridx);
         
         rok=find(~isnan(pp) & ~isnan(rr));
         pp=pp(rok);
         rr=rr(rok);
         
         if isfield(data,'raster'),
            tresp=data.raster(tgoodidx{ridx},:);
            tresp=tresp(rok,:);
         end
        
         % compute correlation between pred and actual
         % response. only look at max lag of 0
         if params.nrandxcov>1,
            [cxy,exy,tt,p]=randxcov(rr,pp,0,params.nrandxcov);
            predxc(strfidx,attidx,ridx)=cxy(1);
            prederr(strfidx,attidx,ridx)=exy(1);
            predp(strfidx,attidx,ridx)=p(1);
         else
            if var(pp)==0 | var(rr)==0,
               cxy=0;
            else
               cxy=xcov(rr,pp,0,'coeff');
            end
            predxc(strfidx,attidx,ridx)=cxy(1);
            prederr(strfidx,attidx,ridx)=nan;
            predp(strfidx,attidx,ridx)=nan;
         end
         
         if 0 && VALIDATIONPHASE,
            %[epochxc,ract,rpred]=epochxc(rpredpsth,resp,stim);
            %pp=tmod_psth(tgoodidx{ridx},ridx);
            %rr=tr(tgoodidx{ridx},ridx);
            
            [predfix(strfidx,attidx,ridx),...
             predfixerr(strfidx,attidx,ridx)]=...
               epochxc(tmod_psth(:,ridx),tr(:,ridx),tstim);
         elseif VALIDATIONPHASE
            %disp('epochxc disabled');
            predfix(strfidx,attidx,ridx)=0;
            predfixerr(strfidx,attidx,ridx)=0;
         end
         
         if VALIDATIONPHASE && length(tresp)==length(pp),
            
            % norm pred power from sahani and linden 2003
            tpsth=nanmean(tresp,2);
            tpred=pp;
            %tpred=tpred-nanmean(tpred);
            %tpred=tpred./nanstd(tpred).*nanstd(tpsth);
            %tpred=tpred+nanmean(tpsth);
            %tresidual=(tpred-nanmean(tpred))./nanstd(tpred).*nanstd(tpsth)-...
            %          (tpsth-nanmean(tpsth));
            tresidual=tpred-tpsth;
            psignal=(repcount*var(tpsth)-nanmean(nanvar(tresp)))./...
                    (repcount-1);
            pnoise=nanmean(nanvar(tresp))-psignal;
            if (repcount*var(tpsth))/nanmean(nanvar(tresp))<1.02,
               disp('less than 2% signal power, setting normpredpower to nan');
               normpredpower=nan;
               normprederr=nan;
            else
               normpredpower=(nanvar(tpsth)-nanvar(tresidual))./ psignal;
               normprederr=prederr(strfidx,attidx,ridx).*normpredpower./...
                   predxc(strfidx,attidx,ridx);
               if normpredpower<0,
                   normpredpower=nan;
                   normprederr=nan;
               end
            end
            predfix(strfidx,attidx,ridx)=normpredpower;
            predfixerr(strfidx,attidx,ridx)=normprederr;
            
            % hsu & theunissen-inspired single trial prediction
            % extrapolation
            xct=zeros(repcount,1);
            for xcidx=1:repcount,
               gg=find(~isnan(tpred) & ~isnan(tresp(:,xcidx)));
               if length(gg)>0 & std(tpred(gg))>0 & std(tresp(gg,xcidx))>0,
                  xct(xcidx)=xcov(tpred(gg),tresp(gg,xcidx),0,'coeff');
               end
            end
            
            projr=sqrt(mean(xct.^2)./rmax.^2);
            fprintf('%.3f/%.3f ---> %.3f (vs. %.3f)\n',...
                    sqrt(mean(xct.^2)),rmax,projr,...
                    predxc(strfidx,attidx,ridx));
            predone(strfidx,attidx,ridx)=sqrt(mean(xct.^2));
            predinf(strfidx,attidx,ridx)=projr;
         end
         
         vrr=var(rr);
         if vrr>0,
            predmse(strfidx,attidx,ridx)=var(pp-rr)./vrr;
         end
         %keyboard
      end
      
   end % for strfidx
   if VALIDATIONPHASE & not(isempty(BATQUEUEID)),
      % record latest (sfsidx,kidx) predicted
      dbsetqueue(BATQUEUEID,1);
   end
end % for attidx

res.predxc=predxc;
res.predone=predone;
res.predinf=predinf;
res.predp=predp;
res.prederr=prederr;
res.predsignal=predsignal;
res.predmse=predmse;
res.predxc2=predxc2;
res.mod_psth=mod_psth;
res.act_resp=act_resp;
res.predfix=predfix;
res.predfixerr=predfixerr;

% svd tweak for demo
if strcmp(strf(1).architecture,'ML'),
    res.stim=tstim;
end
