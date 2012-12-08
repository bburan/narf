% function [mH,eH,pctvar,aH]=cdwrapper(stim,resp,params)
%
% standalone wrapper for cdcore.m
%
% created SVD 9/3/07 - ripped off cellxcnodb.m
%
function [mH,eH,pctvar,H]=cdwrapper(stim,resp,params);

clear ESTIMATIONPHASE VALIDATIONPHASE
global ESTIMATIONPHASE VALIDATIONPHASE
ESTIMATIONPHASE=1;
VALIDATIONPHASE=0;

disp('cdwrapper.m: INITIALIZING');
if ~exist('params','var'),
   error('stim, resp and params arguments required.');
   return
end

% set everything in params structure to default values if they
% haven't been specified
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
params.smoothtime=getparm(params,'smoothtime',0);
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
params.iconside=getparm(params,'iconside',[size(stim,2) 1]);
%=getparm(params,'','');

% BATQUEUEID is the marker to an entry in tQueue, set by queuerun.m
% and cellxcmaster.m so that cellxcnodb can occasionally tell the
% db that the job is still alive and making progress
global BATQUEUEID
if BATQUEUEID<=0,
   BATQUEUEID=[];
end


rsize=size(resp);
times(1).start=1;
times(1).stop=rsize(1);
starttimes=times(1).start;
stoptimes=times(1).stop;
stimstart=starttimes;
stimstop=stoptimes;
spacecount=size(stim,2);
respcount=rsize(2);

if params.fitboot==0,
   %%
   %% NON-BOOTSTRAPPED KERNEL ESTIMATION SECTION
   %% (this usually runs. i haven't debugged the bootstrapped code
   %% recently) 
   %% 
   
   firstseg=1;
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
   
   
   

else
   error('fitboot>0 not yet supported');
   
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
      
      [resp,cdata,fdata]=xcgetbootset(params,bootidx,times,resp0,stim);
      
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
   
   strf=vstrf;
   clear vstrf
end

if ~exist('pctvar','var'),
   pctvar=zeros(params.sfscount,1);
end
