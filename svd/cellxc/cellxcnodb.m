% function cellxcnodb(params);
%
% General flow:
% 0. figure out parameter values (defaults) and determine which
%    files to use in which stage of the analysis (xcfilefracs.m)
% 1. load and preprocess stimulus and response files as specified
%    in stimfiles,respfiles,times,stimloadcmd,stimloadparms,etc.
%    (xcloadfiles.m)
% 2. calculate stimulus-response XC and normalize. (xccore.m)
% 3. fit best svd regularization parameter and shrinkage factor (xcfit.m)
% 4. evaluate kernel by predicting validation data (xcval.m)
%
% INPUTS: params structure (described below)
%
% OUTPUTS: (or at least the important ones)
%  strf(nlidx): structure array where each member represents the
%               STRF with with a different output nonlinearity
%               nlidx=1: none
%                     2: single threshold after filter (best?)
%                     3: sigmoid
%                     4: threshold each spatial channel before sum
%               fields:
%               h: linear(ized) filter (space X time matrix)
%               mS: mean subtracted from each stimulus frame
%                   before applying h (space X 1 vector)
%               nltype: type of output nl
%               nlparms: parameters for output nl (eg, threshold)
%               so:      rpred= NL( h * (stim-mS) )
%               parms.kernfmt: code for linearizing transform indicating
%                              appriate display format for showkern
% predres:  results of predictions on validatation data
% 
% Fields for params (default value in parens):
% params.stimfiles : {list} of stimulus files *REQUIRED* 
% params.stimloadcmd : matlab command to load stim ('loadimfile')
%                stimloadcmd(file,startframe,stopframe,p1,p2,...)
% params.stimloadparms : {p1,p2,...} parameters passed to stimloadcmd ([])
% params.stimwindowcrf : special parameter for stimloadparms (0)
% params.stimcrfs : special parameter for stimloadparms ([0 ... ])
% params.stimfiltercmd : filter applied to stim after load ('')
% params.stimfilterparms : parameters for stimfiltercmd ([]);
% params.respfiles : {list} of response files *REQUIRED*
% params.resploadcmd : command to load responses ('respload')
% params.resploadparms : parameters to resploadcmd ({'',1,1,1})
% params.respfiltercmd : command to filter response ('')
% params.respfilterparms : parameters to respfiltercmd ([])
% params.kernfmt : string describing spatial format ('space')
%                  eg, 'pfft','sfft','wav'
% params.minlag : min temporal latency of filter (-6)
% params.maxlag : maxiumum temporal latency of filter (13)
% params.resampcount : number of bootstrap (20);
% params.resampfmt : 0 no bootstrap, 1 bootstrap (1)
% params.expfrac : fraction for regression (1-fitfrac-predfrac)
% params.fitfrac : fraction for fitting regularization parms (0)
%                  (zero means use exp data)
% params.predfrac : fraction to reserve for validation (0.1)
% params.decorrspace : order of spatial decorrelation (2)
% params.decorrtime : order of temporal decorrelation (1)
% params.sffiltsigma : number of shrinkage values to test (7)
% params.sfscount : number of regularization parms to test (30)
% params.sfsstep : order of magnitude (range) of reg parms to test (4)
% params.smoothtime : 0 don't smooth STRF in time
%                     1 smooth STA in time
%                     2 smooth decorrelated STRFs in time   (2)
% params.predsmoothsigma : std of gaussian for smoothing responses (0)
% params.predtype : uh... dunno what this is for (0)
% params.nloutparm : number of output NLs to test (3)
% params.sffiltthresh : nothing, i think (0)
% params.sffiltsmooth : ditto (0)
% params.boundary : how to deal with boundaries in XC ('zero')
% params.repexclude : 0 single iteration of regression 
%      1 remove spatial channels that aren't significant, ARD-ish (0)
% params.DOCELLFIT2 : 0 do shrinkage then fit regularization parms
%                     1 =getparm(params,'docellfit2',0);
% params.stimfmtcode=getparm(params,'stimfmtcode',0);
% params.respfmtcode=getparm(params,'respfmtcode',0);
% params.nlidxsave : output nl to save 1=none, 2=single thresh,
%                    3=thresh each spatial channel
%                    4=rect each channel (2)
% params.outfile : save output to this file ('/tmp/cellxcout.mat')
% params.showres : 1=run xcresult(params.outfile) at end (default=1)
%                : 0=skip results display 
% params.cellid : name of cell for labeling ('CELL')
% params.shrinkage : 0 threshold 1 shrinkage (1)
% params.predbatch : ({}) used by cellxcqueue for multiple
%                    validation sets
% params.fitboot : if 1, bootstrap validation (0)
%
% TODOs:
% 1. make sure fitboot parameter is working up to snuff
% 2. cope with big stim/response sets by loading partial segments
% into memory. this would require a loop around xcloadfiles and
% xccore. probably would make sense to pull out the matrix
% initialization and normalization from xccore and put them in
% routines on the outside of this loop.
%
% created SVD 3/10/03 - ripped off of cellxc3.m
%
function cellxcnodb(params);

VCELLXC=3;

clear ESTIMATIONPHASE VALIDATIONPHASE
global ESTIMATIONPHASE VALIDATIONPHASE
ESTIMATIONPHASE=1;
VALIDATIONPHASE=0;

disp('cellxcnodb.m: INITIALIZING');
if ~exist('params','var'),
   error('param argument required.');
   return
end

% set everything in params structure to default values if they
% haven't been specified
params.stimfiles=getparm(params,'stimfiles','ERROR');
params.stimcrfs=getparm(params,'stimcrfs',zeros(1,length(params.stimfiles)));
params.stimloadcmd=getparm(params,'stimloadcmd','loadimfile');
params.stimloadparms=getparm(params,'stimloadparms',[]);
params.stimfiltercmd=getparm(params,'stimfiltercmd','');
params.stimfilterparms=getparm(params,'stimfilterparms',[]);
params.respfiles=getparm(params,'respfiles','ERROR');
params.resploadcmd=getparm(params,'resploadcmd','respload');
params.resploadparms=getparm(params,'resploadparms',{'',1,1,1});
params.respfiltercmd=getparm(params,'respfiltercmd','');
params.respfilterparms=getparm(params,'respfilterparms',[]);
params.kernfmt=getparm(params,'kernfmt','space');
if isfield(params,'maxlag') & length(params.maxlag)>=2,
   % do nothing, this is a good format for running cellxc
else
   params.maxlag=[getparm(params,'minlag',-6) getparm(params,'maxlag',13)];
end
params.resampcount=getparm(params,'resampcount',20);
params.resampfmt=getparm(params,'resampfmt',1);
params.expfrac=getparm(params,'expfrac',0);
params.fitfrac=getparm(params,'fitfrac',0);
params.predfrac=getparm(params,'predfrac',0.1);
params.decorrspace=getparm(params,'decorrspace',2);
params.decorrtime=getparm(params,'decorrtime',1);
params.sffiltsigma=getparm(params,'sffiltsigma',7);
params.sfscount=getparm(params,'sfscount',30);
params.sfsstep=getparm(params,'sfsstep',4);
params.smoothtime=getparm(params,'smoothtime',2);
params.predsmoothsigma=getparm(params,'predsmoothsigma',0);
params.predtype=getparm(params,'predtype',0);
params.stimwindowcrf=getparm(params,'stimwindowcrf',2);
params.nloutparm=getparm(params,'nloutparm',3);
params.sffiltthresh=getparm(params,'sffiltthresh',0);
params.sffiltsmooth=getparm(params,'sffiltsmooth',0);
params.boundary=getparm(params,'boundary','zero');
params.repexclude=getparm(params,'repexclude',0);
params.docellfit2=getparm(params,'docellfit2',0);
params.stimfmtcode=getparm(params,'stimfmtcode',0);
params.respfmtcode=getparm(params,'respfmtcode',0);
params.nlidxsave=getparm(params,'nlidxsave',2);
params.outfile=getparm(params,'outfile','/tmp/cellxcout.mat');
params.showres=getparm(params,'showres',0);
params.zipoutfile=getparm(params,'zipoutfile',0);
params.cellid=getparm(params,'cellid','CELL');
params.shrinkage=getparm(params,'shrinkage',1);
params.predbatch=getparm(params,'predbatch',{});
params.batch=getparm(params,'batch',[]);
params.fitboot=getparm(params,'fitboot',0);
params.meansub=getparm(params,'meansub',1);
params.tbinms=getparm(params,'tbinms',16);
params.nrandxcov=getparm(params,'nrandxcov',200);
params.sharpspacenorm=getparm(params,'sharpspacenorm',0);
params.keepneg=getparm(params,'keepneg',0);
params.cutextendedsilence=getparm(params,'cutextendedsilence',0);
%=getparm(params,'','');

% BATQUEUEID is the marker to an entry in tQueue, set by queuerun.m
% and cellxcmaster.m so that cellxcnodb can occasionally tell the
% db that the job is still alive and making progress
global BATQUEUEID
if BATQUEUEID<=0,
   BATQUEUEID=[];
end

% given the list of response files and exp-/fit-/predfracs, figure
% out the relevant start and stop times of each rc segment

if isfield(params,'times'),
   times=params.times;
   
elseif params.fitboot,
   tparams=params;
   tparams.fitfrac=0;
   %tparams.predfrac=0;
   times=xcfilefracs(tparams);
   %times(1).stop(end)=times(1).stop(end)+1;
   params.times=times;
else
   
   % call xcfilefracs to get times
   times=xcfilefracs(params);
   params.times=times;
end

starttimes=times(1).start;
stoptimes=times(1).stop;
fitfile=times(2).fileidx;
fitstartframe=times(2).start;
fitstopframe=times(2).stop;
predfile=times(3).fileidx;
predstartframe=times(3).start;
predstopframe=times(3).stop;

attcount=size(times(1).start,2);  % ignore attcount>1 for the time being
attidx=1;

if params.fitboot==0,
   %%
   %% NON-BOOTSTRAPPED KERNEL ESTIMATION SECTION
   %% (this usually runs. i haven't debugged the bootstrapped code
   %% recently) 
   %% 
   
   % for big memory jobs, stick all this stuff in a loop.
   xcloadfiles;
   % keyboard
   
   if isfield(params,'spacetimesep') & params.spacetimesep,
      
      % space time separable model. xcsepcore does xccore and xcfit
      % in one fell swoop
      nlidxsavefinal=params.nlidxsave;
      params.nlidxsave=2;
      xcsepcore;
      params.nlidxsave=nlidxsavefinal;
   else

      if params.repexclude,
         % do repeated exclusions
         xcitercore;
      elseif isfield(params,'usecd') & params.usecd,
         cdcore;
      elseif isfield(params,'altcore'),
         eval(params.altcore);
      else
         % simply do the XC. standard
         xccore;
      end
      
      %%
      %% FIT SECTION - FIND OPTIMAL REGULARIZATION AND SHRINKAGE PARAMETERS
      %%
      if isfield(params,'altfit'),
         eval(params.altfit);
      else
         xcfit;
      end
      
      %if isfield(params,'usecd') & params.usecd,
      %   keyboard
      %end
   end
   
   
else
   
   xcloadfiles;
   times(1).start=stimstart;
   times(1).stop=stimstop;
 
   resp0=resp;
   
   filelens=stimstop-stimstart+1;
   filecount=max(find(filelens>1));
   filelens=filelens(1:filecount);
   
   totlen=sum(filelens);
   clen=[0; cumsum(filelens)];
   
   bootcount=params.fitboot;
   bootstep=filelens./bootcount;
   spacer=diff(params.maxlag);
   vexpxc=[];
   clear bootresults
   
   for bootidx=1:bootcount,
      
      [resp,cdata,fdata]=xcgetbootset(params,bootidx,times,resp0,stim,raster);
      %keyboard
      if 0 & strcmp(params.respfiltercmd,'resptakefracs'),
         respfracs=params.respfilterparms{1};
         goodidx=find(~isnan(resp(:,1)));
         goodlen=length(goodidx);
         shift=round(goodlen./bootcount.* ...
                     mod(-bootidx+round(bootcount/4),bootcount));
         resp(:,2:end)=nan;
         for ff=2:length(respfracs),
            keepidx=1:round(length(goodidx).*respfracs(ff));
            keepidx=mod(keepidx+shift-1,goodlen)+1;
            resp(goodidx(keepidx),ff)=resp(goodidx(keepidx),1);
         end
      end
      
      firstseg=1;
      if isfield(params,'spacetimesep') & params.spacetimesep,
         nlidxsavefinal=params.nlidxsave;
         params.nlidxsave=2;
         xcsepcore;
         params.nlidxsave=nlidxsavefinal;
      else
         if params.repexclude,
            % do repeated exclusions
            xcitercore;
         elseif isfield(params,'altcore'),
            eval(params.altcore);
         else
            % simply do the XC. standard
            xccore;
         end
         
         if isfield(params,'altfit'),
            eval(params.altfit);
         else
            xcfit;
         end
      end
      
      % save imporant stuff from this bootidx
      %keyboard
      
      vexpxc=cat(4,vexpxc,expxc);
      vstrf(:,bootidx,:)=strf;
      
      % only do validation here if cdata changes with each boot iteration
      if params.predfrac==0,
         predres(bootidx)=xcval(strf,params,cdata);
         
         if bootidx==1,
            rrfull=zeros(size(resp)).*nan;
            ppfull=zeros(length(resp),...
                         size(predres(bootidx).mod_psth{1},3)).*nan;
         end
         
         rrfull(cdata.fromidx)=predres(bootidx).act_resp{1}(cdata.toidx);
         ppfull(cdata.fromidx,:)=predres(bootidx).mod_psth{1}(cdata.toidx,:);
      end
      
   end
   resp=resp0;
   strf=vstrf;
   clear vstrf resp0
end

clear sH H ttSA tsSA2 tsSAfull tSR tmR tn tSA1
clear tSA sSA1 sSA2 sSAfull sSA0 SR
clear fdata

%%
%% VALIDATION SECTION - PREDICT NOVEL RESPONSES
%%

ESTIMATIONPHASE=0;
VALIDATIONPHASE=1;

batchcount=length(params.predbatch);

if batchcount==0 && params.fitboot>0 && params.predfrac==0,
   % do nothing, preds already done.
   % now skip down and piece into standard batchcount>0 area
   
elseif batchcount==0 && params.predfrac>0,
   
   %
   % this code usually runs for validation unless you're doing
   % something funky with multiple stimulus types
   %
   
   if isempty(BATQUEUEID),
      %keyboard
   end
   
   % just assume current batch
   if ~isempty(params.batch),
      params.predbatch={params.batch};
   else
      params.predbatch={};  % "not a batch"
   end
   
   nlcount=size(strf,1);
   attcount=size(strf,2);
   respcount=size(strf,3);
   
   if params.fitboot,
      predres=xcval(strf(:),params,times(3));
   else
      predres=xcval(reshape(permute(strf,[1 3 2]),nlcount*respcount,...
                            attcount),params,times(3));
   end
elseif batchcount>0,
   
   if params.fitboot,
      rrfull={};
      ppfull={};
      clear predres;
   end
   
   % multiple validation run classes
   for predidx=1:batchcount,
      
      fprintf('Predicting responses to batch %d.\n',params.predbatch{predidx});
      
      % figure out pred files for the current batch
      [pcellfiledata,ptimes,pbatchdata]=...
          cellfiletimes(params.cellid,params.predbatch{predidx});
      
      % does this cell have data for batchid=predidx?
      if length(pcellfiledata)>0,
         
         % predict whole response for datasets with zero predfrac
         if pbatchdata.predfrac==0,
            ptimes(3)=ptimes(1);
         end
         
         predparams=params;
         if params.batch~=params.predbatch{predidx},
            predparams.stimfiles=pbatchdata.stimfiles;
            predparams.respfiles=pbatchdata.respfiles;
            predparams.stimcrfs=pbatchdata.stimcrfs;
            
            predparams.stimloadcmd=pbatchdata.stimloadcmd;
            predparams.stimloadparms=pbatchdata.stimloadparms;
            if strcmp(predparams.stimloadcmd,'loadsiteraster'),
               predparams.stimloadparms{1}.channel=params.stimloadparms{1}.channel;
               predparams.stimloadparms{1}.unit=params.stimloadparms{1}.unit;
            end
            predparams.stimfiltercmd=pbatchdata.stimfiltercmd;
            if ~isempty(pbatchdata.stimfilterparms),
               pbatchdata.stimfilterparms=...
                   strrep(pbatchdata.stimfilterparms,'cellfiledata','pbatchdata');
            end
            predparams.stimfilterparms=strsep(pbatchdata.stimfilterparms,',');
            predparams.resploadcmd=pbatchdata.resploadcmd;
            predparams.resploadparms=pbatchdata.resploadparms;
            predparams.respfiltercmd=pbatchdata.respfiltercmd;
            predparams.respfilterparms=pbatchdata.respfilterparms;
            predparams.times=ptimes;
            
% $$$             % allow for exection of cellid-specific identifiers in strsep
% $$$             if ~isempty(pbatchdata.stimloadparms),
% $$$                pbatchdata.stimloadparms=...
% $$$                    strrep(pbatchdata.stimloadparms,'cellfiledata','pbatchdata');
% $$$             end
% $$$             predparams.stimloadparms=strsep(pbatchdata.stimloadparms,',');
% $$$             predparams.stimfiltercmd=pbatchdata.stimfiltercmd;
% $$$             if ~isempty(pbatchdata.stimfilterparms),
% $$$                pbatchdata.stimfilterparms=...
% $$$                    strrep(pbatchdata.stimfilterparms,'cellfiledata','pbatchdata');
% $$$             end
% $$$             predparams.stimfilterparms=strsep(pbatchdata.stimfilterparms,',');
% $$$             predparams.resploadcmd=pbatchdata.resploadcmd;
% $$$             predparams.resploadparms=strsep(pbatchdata.resploadparms,',');
% $$$             predparams.respfiltercmd=pbatchdata.respfiltercmd;
% $$$             predparams.times=ptimes;
% $$$             predparams.respfilterparms=strsep(pbatchdata.respfilterparms,',');
         end
         
         tpredstartframe=ptimes(3).start;
         tpredstopframe=ptimes(3).stop;
         tpredfile=ptimes(3).fileidx;
         cdata=[];         
         [cdata.stim,cdata.resp,extras]=...
             xcloadstimresp(tpredfile,tpredstartframe,...
                            tpredstopframe,predparams);
         if isfield(extras,'raster'),
            cdata.raster=extras.raster;
         end
      else
         cdata.stim=[];
         cdata.resp=[];
         predparams=params;
         %if predidx>1,
         %   predres(predidx).predxc=[];
         %end
      end
      nlcount=size(strf,1);
      attcount=size(strf,2);
      respcount=size(strf,3);
      
      if params.fitboot,
         for bootidx=1:params.fitboot,
            if length(cdata.resp)>0
               [tresp,tcdata]=xcgetbootset(predparams,bootidx,ptimes,...
                                           cdata.resp,cdata.stim,raster);
               
               predres(predidx,bootidx)=...
                   xcval(strf(:,bootidx),predparams,tcdata);
               
               if bootidx==1,
                  %keyboard
                  rrfull{predidx}=zeros(size(tresp)).*nan;
                  ppfull{predidx}=zeros(length(tresp),...
                          size(predres(predidx,bootidx).mod_psth{1},3)).*nan;
               end
               
               rrfull{predidx}(tcdata.fromidx)=...
                   predres(predidx,bootidx).act_resp{1}(tcdata.toidx);
               ppfull{predidx}(tcdata.fromidx,:)=...
                   predres(predidx,bootidx).mod_psth{1}(tcdata.toidx,:)-...
                   repmat(nanmean(predres(predidx,bootidx).mod_psth{1}(tcdata.toidx,:)),[length(tcdata.toidx),1]);
            else
               predres(predidx,bootidx)=xcval(strf(:,bootidx),predparams,cdata);
            end
         end
         attcount=size(predres(1).predxc,2); 
      else
          %keyboard
         predres(predidx)=xcval(reshape(permute(strf,[1 3 2]),...
                                        nlcount*respcount,attcount),...
                                predparams,cdata);
      end
   end
else
   predres.predxc=[];
   predres.predone=[];
   predres.predinf=[];
   predres.predp=[];
   predres.prederr=[];
   predres.predfix=[];
   predres.predfixerr=[];
   predres.predmse=[];
end

% consoldiate predres into compact matrices for saving to sResults
% in the database.
predxc=cat(4,predres.predxc);
predinf=cat(4,predres.predinf);
predone=cat(4,predres.predone);
predp=cat(4,predres.predp);
prederr=cat(4,predres.prederr);
predfix=cat(4,predres.predfix);
predfixerr=cat(4,predres.predfixerr);
predmse=cat(4,predres.predmse);

bcount=size(predxc,4);
if params.predfrac>0 & ~params.fitboot,
   predxc=permute(reshape(predxc,nlcount,respcount,size(predxc,3),...
                          bcount),[4 1 3 2]);
   predinf=permute(reshape(predinf,nlcount,respcount,size(predinf,3),...
                          bcount),[4 1 3 2]);
   predone=permute(reshape(predone,nlcount,respcount,size(predone,3),...
                          bcount),[4 1 3 2]);
   predp=permute(reshape(predp,nlcount,respcount,size(predp,3),...
                         bcount),[4 1 3 2]);
   prederr=permute(reshape(prederr,nlcount,respcount,size(prederr,3),...
                           bcount),[4 1 3 2]);
   predmse=permute(reshape(predmse,nlcount,respcount,...
                           size(predmse,3),bcount),[4 1 3 2]);
   predfix=permute(reshape(predfix,nlcount,respcount,...
                           size(predfix,3),bcount),[4 1 3 2]);
   predfixerr=permute(reshape(predfixerr,nlcount,respcount,...
                              size(predfixerr,3),bcount),[4 1 3 2]);
   
elseif params.predfrac>0,
   
   % predres ( batch x 1) . predxc ( nlcount*bootcount*respcount x 1)
   
   % vpredxc should be batch X nl X boot X resp
   % predxc should be batch X nl X 1 X resp   (averaged over boots!)
   vpredxc=permute(reshape(predxc,nlcount,attcount,respcount,...
                           bcount),[4 1 2 3]);
   vpredinf=permute(reshape(predinf,nlcount,attcount,respcount,...
                           bcount),[4 1 2 3]);
   vpredone=permute(reshape(predone,nlcount,attcount,respcount,...
                           bcount),[4 1 2 3]);
   vpredp=permute(reshape(predp,nlcount,attcount,respcount,...
                          bcount),[4 1 2 3]);
   vprederr=permute(reshape(prederr,nlcount,attcount,...
                            respcount,bcount),[4 1 2 3]);
   vpredmse=permute(reshape(predmse,nlcount,attcount,...
                            respcount,bcount),[4 1 2 3]);
   vpredfix=permute(reshape(predfix,nlcount,attcount,...
                            respcount,bcount),[4 1 2 3]);
   vpredfixerr=permute(reshape(predfixerr,nlcount,attcount,...
                               respcount,bcount),[4 1 2 3]);
   
   predxc=mean(vpredxc,3);
   predinf=mean(vpredinf,3);
   predone=mean(vpredone,3);
   predp=mean(vpredp,3);
   prederr=mean(vprederr,3);
   predmse=mean(vpredmse,3);
   predfix=mean(vpredfix,3);
   predfixerr=mean(vpredfixerr,3);
   
   expxc=mean(vexpxc,4);
   
elseif params.fitboot,
   
   % predres ( batch x 1) . predxc ( nlcount*bootcount*respcount x 1)
   
   % vpredxc should be batch X nl X boot X resp
   % predxc should be batch X nl X 1 X resp   (averaged over boots!)
   
   vpredxc=permute(reshape(predxc,nlcount,attcount,respcount,...
                           bcount),[4 1 2 3]);
   vpredinf=permute(reshape(predinf,nlcount,attcount,respcount,...
                           bcount),[4 1 2 3]);
   vpredone=permute(reshape(predone,nlcount,attcount,respcount,...
                           bcount),[4 1 2 3]);
   vpredp=permute(reshape(predp,nlcount,attcount,respcount,...
                          bcount),[4 1 2 3]);
   vprederr=permute(reshape(prederr,nlcount,attcount,...
                            respcount,bcount),[4 1 2 3]);
   vpredmse=permute(reshape(predmse,nlcount,attcount,...
                            respcount,bcount),[4 1 2 3]);
   vpredfix=permute(reshape(predfix,nlcount,attcount,...
                            respcount,bcount),[4 1 2 3]);
   vpredfixerr=permute(reshape(predfixerr,nlcount,attcount,...
                               respcount,bcount),[4 1 2 3]);
   
   if batchcount==0,
      batchcount=1;
      rr=rrfull;
      pp=ppfull;
   else
      rr=rrfull{1};
      pp=ppfull{1};
   end
   predxc=permute(mean(reshape(vpredxc,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   predinf=permute(mean(reshape(vpredinf,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   predone=permute(mean(reshape(vpredone,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   predp=permute(mean(reshape(vpredp,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   prederr=permute(mean(reshape(vprederr,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   predmse=permute(mean(reshape(vpredmse,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   predfix=permute(mean(reshape(vpredfix,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   predfixerr=permute(mean(reshape(vpredfixerr,batchcount,params.fitboot,nlcount),2),[1 3 2]);
   
   expxc=mean(vexpxc,4);

   %disp('untangling prediction data');
   %for predidx=1:length(rrfull),
   %   rr=rrfull{predidx};
   %   pp=ppfull{predidx};
   %   
   %   for ii=1:size(pp,2),
   %      kk=find(~isnan(rr)&~isnan(pp(:,ii)));
   %      [predxc(predidx,ii),prederr(predidx,ii),tt,predp(predidx,ii)]=randxcov(rr(kk),pp(kk,ii),0,100);
   %      vrr=var(rr(kk));
   %      if vrr>0,
   %         predmse(predidx,ii)=var(pp(kk,ii)-rr(kk))./vrr;
   %      end
   %   end
   %end
else
   predxc=ones([1 size(strf)]).*nan;
   predinf=ones([1 size(strf)]).*nan;
   predone=ones([1 size(strf)]).*nan;
   predp=zeros([1 size(strf)]);
   prederr=ones([1 size(strf)]).*nan;
   predmse=ones([1 size(strf)]).*nan;
   predfix=ones([1 size(strf)]).*nan;
   predfixerr=ones([1 size(strf)]).*nan;
end

clear tstim tr th tH seplinpred rgoodidx respsave linpred tstrf
clear tmod_psth tgoodidx rpred rbadidx psthval lengthVec 
clear fstim fstim_spike fitidx CSR_ns CS_ns cdata tresp valpred

global SUBXC
if ~isempty(SUBXC),
   subxc=SUBXC;
   subxc(end,:,:)=squeeze(predxc);
end

% save results, zipped if necessary
daterun=date;
if params.zipoutfile,
   skernfile=basename(params.outfile);
   fprintf('cellxcnodb.m: SAVING to %s\n',params.outfile);
   
   % then save everything (inclusive to avoid future screw-ups
   % involving accidentally not saving new important stuff)
   save([tempdir,skernfile]);
   
   if params.showres,
      xcresult([tempdir,skernfile],[],1);
   end
   
   % compress and copy over network if necessary
   unix(['gzip -cf ',tempdir,skernfile,' > ',...
         params.outfile,'.gz']);
   delete([tempdir,skernfile]);
else
   save(params.outfile);
   if ~params.showres
      disp('skipping figure output');
   elseif params.batch>0,
      xcresult(params.outfile,[],1);
   else
      xcresult(params.outfile); % no batch--don't print
   end
end

disp('cellxcnodb.m: DONE');



