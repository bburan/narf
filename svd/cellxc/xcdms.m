% xcdms.m: remake of kerncompN.m to make sure it works and to
% make the code/analysis a bit cleaner.
%
% created SVD 7/1/04 (hacked coarsely from kernvsatt2.m)
%
% BASIC FLOW:
% 1. load relevant stimuli and resposnes (these can be jamie's fvvs
%    data or ben's dms data)
%
% 2. reshuffle stimuli so as to get mean and difference responses.
% 3. fit various models:
% 4. predict validation data for each jackknife and randomized
%    attention
% 5. test for significant improvements with each model component.
%
% created 8/04 - adapted from kernvsatt2.m for dms data
%
disp('xcdms.m : dms attention tester');

if isfield(params,'maxlag') & length(params.maxlag)>=2,
   % do nothing, this is a good format for running cellxc
else
   params.maxlag=[getparm(params,'minlag',0) getparm(params,'maxlag',0)];
end
params.meansub=getparm(params,'meansub',1);
params.tbinms=getparm(params,'tbinms',16);
params.nrandxcov=getparm(params,'nrandxcov',200);
params.randexp=getparm(params,'randexp',0);
params.randcnf=getparm(params,'randcnf',1);
params.runclassid=getparm(params,'runclassid',10);
params.noisecount=getparm(params,'noisecount',200);
params.dodamping=getparm(params,'dodamping',0);
params.bootcount=getparm(params,'bootcount',20);
params.spacelims=getparm(params,'spacelims',128);
params.sfscount=getparm(params,'sfscount',25);
params.sfsstep=getparm(params,'sfsstep',5);
params.sffiltsigma=getparm(params,'sffiltsigma',5);
params.decorrspace=getparm(params,'decorrspace',2);
params.decorrtime=getparm(params,'decorrtime',0);
params.resampcount=getparm(params,'resampcount',20);
params.resampfmt=getparm(params,'resampfmt',1);
params.smoothtime=getparm(params,'smoothtime',0);
params.sharpspacenorm=getparm(params,'sharpspacenorm',0);
params.shrinkage=getparm(params,'shrinkage',1);

VCELLXC=3;
clear ESTIMATIONPHASE VALIDATIONPHASE
global ESTIMATIONPHASE VALIDATIONPHASE
ESTIMATIONPHASE=1;
VALIDATIONPHASE=0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. LOAD ALL STIM AND RESPONSE FILES AND APPLY THE APPROPRIATE
%    INPUT NL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if params.runclassid==10, % ie, dms
   z=load(params.respfiles{1});
   bigpatches=z.target;
   [fmask,crop]=movfmask(size(z.target,1),0,params.stimloadparms{3}*...
                         size(z.target,1)./params.stimloadparms{1});
   targpatches=movresize(floor(z.target),params.stimloadparms{3},fmask,crop);
   ftargs=feval(params.stimfiltercmd,targpatches,...
                params.stimfilterparms{:});
   
   mpatches=[mean(ftargs,2),ftargs];
   fpatches=ftargs;
   
else   %runclass==4 -- fvvs
   loadfreetargs;
end

btarglist=-1;
bresp=[];
bstim=[];
filecount=length(params.respfiles);
% for each file:
for fidx=1:filecount,
   
   % check to see if params.stim and response exist
   if ~exist(params.respfiles{fidx},'file'),
      disp([params.respfiles{fidx},' missing!!']);
      return
   end
   if ~exist(params.stimfiles{fidx},'file'),
      disp([params.stimfiles{fidx},' missing!!']);
      return
   end
   
   % load stimulus
   tstim=feval(params.stimloadcmd,params.stimfiles{fidx},1,0,...
               params.stimloadparms{:});
   
   % filter stimulus segment if selected
   if ~isempty(params.stimfiltercmd),
      tstim=feval(params.stimfiltercmd,tstim,params.stimfilterparms{:});
   end
   % reshape to space X time if necessary
   iconside=size(tstim);
   iconside=iconside(1:(end-1));
   if length(iconside)>=2,
      % reshape all spatial dims into one
      tstim=reshape(tstim,prod(iconside),size(tstim,length(iconside)+1));
   end
   tstim=tstim'; % take transpose to put time in rows
   
   bstim=cat(1,bstim,tstim);
   
   %
   % figure out atten states:
   % 
   z=load(params.respfiles{fidx});
   if isfield(z,'phase'),
      % FVVS
      %validtargs=unique(z.targets(find(z.phase>-1 & ~z.result)));
      validtargs=unique(z.targets);
      otarglist=btarglist;
      btarglist=union(btarglist,validtargs);
      
      if length(btarglist) > length(otarglist) & size(bresp,1)>0 & ...
            params.resploadparms{2}==1,
         
         oresp=bresp;
         bresp=ones(size(bresp,1),size(bresp,2),length(btarglist)).*nan;
         for ii=1:length(otarglist),
            newidx=find(otarglist(ii)==btarglist);
            bresp(:,:,newidx)=oresp(:,:,ii);
         end
      end
   else
      % DMS
      btarglist=[12 1 2];
   end
   
   fprintf('Current targlist:');
   for ii=1:length(btarglist),
      fprintf(' %d',btarglist(ii));
   end
   fprintf('\n');
   
   %
   % load the response data
   %
   disp(sprintf('Loading response: %s...',params.respfiles{fidx}));
   % resp=resploadatt(respfile,respvarname);
   % resp is time X space(phase) X attcode
   tresp=feval(params.resploadcmd,params.respfiles{fidx},...
               params.resploadparms{:});
   
   % filter response (eg resample, pick attentional state, etc)
   if ~isempty(params.respfiltercmd),
      tresp=feval(params.respfiltercmd,tresp,params.respfilterparms{:});
   end
   
   if size(bresp,1)>0 & params.resploadparms{2}==1 & isfield(z,'phase'),
      
      oresp=tresp;
      tresp=tresp.*nan;
      otarglist=[-1; z.targlist(:)];
      for ii=1:length(otarglist),
         newidx=find(otarglist(ii)==btarglist);
         tresp(:,:,newidx)=oresp(:,:,ii);
      end
      bresp=cat(1,bresp,tresp);
   elseif size(bresp,1)>0,
      
      bresp=cat(1,bresp,tresp);
   else
      bresp=tresp;
   end
end

accounts=sum(~isnan(bresp(:,1,:)),1);

%disp('including stim not shown in all att conds!');
acountgood=find(accounts>0);
bresp=bresp(:,:,acountgood);

clear z

% figure out size of response
blen=size(bresp,1);
respcount=size(bresp,2);
spacecount=size(bstim,2);

% rearrange responses for different attention models
%
% models to test:
%   1: nanmean(ABab)
%   2: nanmean(Aa) - nanmean(Bb)
%   3: nanmean(AB) - nanmean(ab)
%   4: A - B

OUTNLMODE=5;
nlnames={'ABab','Aa-Bb','AB-ab','A-B'};

bresp0=bresp;
attcount=4;
bresp=zeros(size(bresp));

% subtract mean from each attention condition and manipulate means
% independently of indiv stim-resp samples to minimize bias from
% data that was collected from less than four att conds

disp('normalizing mean separately in recombined responses');
gidx=find(sum(isnan(bresp0),3)==0);
mm=nanmean(bresp0(gidx,:,:));
fprintf('m_A=%.2f m_B=%.2f m_a=%.2f m_b=%.2f\n',mm);
bresp0=bresp0-repmat(mm,[size(bresp0,1) 1]);

% include data not shown in all trials
bresp(:,:,1)=permute(nanmean(permute(bresp0,[3 1 2])),[2 3 1]) + ...
    mean(mm);
bresp(:,:,2)=permute(nanmean(permute(bresp0(:,:,[1 3]),[3 1 2]))-...
                     nanmean(permute(bresp0(:,:,[2 4]),[3 1 2])), ...
                     [2 3 1]) + mean(mm([1 3])) - mean(mm([2 4]));
bresp(:,:,3)=permute(nanmean(permute(bresp0(:,:,[1 2]),[3 1 2]))-...
                     nanmean(permute(bresp0(:,:,[3 4]),[3 1 2])),...
                     [2 3 1]) + mean(mm([1 2])) - mean(mm([3 4]));
bresp(:,:,4)=bresp0(:,:,1)-bresp0(:,:,2) + mm(1)-mm(2);

% divide by 2 to make model cleaner (see detrending stuff below)
bresp(:,:,2:end)=bresp(:,:,2:end)./2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. BREAK INTO JACKKNIFE SETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find all time bins with valid responses in an attention state > 1
% (the first set is ALL attention ... throw away data not in the
% other categories)
%

% figure out number of valid fixations in each attentional condition
ncount=squeeze(sum(~isnan(bresp(:,end,:)),1));
anyokidx=find(sum(~isnan(bresp(:,end,:)),3));

% choose a random set of fixations for each bootstrap. this should
% avoid bias in the non-randomized noiseidx=1
vcount=length(anyokidx);
vidx=anyokidx;

origresp=bresp;

if 1,
   %
   % detrend difference data so as not to confound dc/gain with
   % local shifts
   %
   disp('detrending difference data!');
   
   %model:  given mr=(rAB+rab)/2 and dr=(rAa-rBb)/2
   %
   % fit dr=mr * globalgain + globaldc
   % then fit diff srf to dr'=dr-(mr*globalgain+globaldc)
   %
   % to reconstruct:
   % rAa = mr + (mr*globalgain+globaldc) + dr'
   % rBb = mr - (mr*globalgain+globaldc) - dr'
   
   globalgain=zeros(attcount,params.bootcount);
   globaldc=zeros(attcount,params.bootcount);
   globaldconly=zeros(attcount,params.bootcount);
   
   for ii=2:attcount,
      for bootidx=1:params.bootcount;
         
         % generate exploratory response matrix that is resp
         % with cnf data nanned out
         ppidx=vidx([1:(round((bootidx-1)/params.bootcount*vcount)) ...
                     round(bootidx/params.bootcount*vcount+1):end]);
         
         x=bresp(ppidx,1,1);
         y=bresp(ppidx,1,ii);
         
         gidx=find(~isnan(x+y));
         x=x(gidx);
         y=y(gidx);
         
         mpg=polyfit(x,y,1);
         globaldc(ii,bootidx)=mpg(2);
         globalgain(ii,bootidx)=mpg(1);
         globaldconly(ii,bootidx)=mean(y);
      end
      
      % remove mean/gain from difference responses
      
      m=mean(globaldc(ii,:));
      g=mean(globalgain(ii,:));
      
      y=bresp(:,1,ii)-g.*bresp(:,1,1)-m;
      
      bresp(:,1,ii)=y;
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. CROSS-CORRELATE FOR VARIOUS MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% output will be saved to vstrf structure... mimicking strategy of
% cellxc standard analysis

% rendexp is not written in this code, so kerncount should always be 1
if params.randexp,
   kerncount=params.noisecount+1;
else
   kerncount=1;
end

% use for non-cheating local att preds
resplen=size(bresp,1);
valattpreds=zeros(resplen,attcount).*nan;

% cross correlate and apply bias correction, leaving in eig domain
for noiseidx=1:kerncount,
   fprintf('\nnoiseidx=%d: ',noiseidx);
   
   % resp is temporary matrix with attention assigned according to
   % noise (for randomized attention kernels)
   
   eresp=bresp;
   
   % assign resp values from brespnoglob (ie, with dc and gain removed
   % from specific attention state responses)
   if noiseidx==1,
      eresp(:,:,2:end)=bresp(:,:,2:end);
   else
      disp('rand kernel atts not implemented yet!');
   end
   
   vexpxc=[];
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % DOUBLE JACKKNIFED RC
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   rbinidx=1;
   for bootidx=1:params.bootcount;
      fprintf('\nbootidx=%d/%d:\n',bootidx,params.bootcount);
      
      % generate exploratory response matrix that is resp
      % with cnf data nanned out
      ppidx=vidx((round((bootidx-1)/params.bootcount*vcount)+1):...
                 round(bootidx/params.bootcount*vcount));
      
      resp=eresp;
      resp(ppidx,:,:)=nan;
      
      tokidx=find(~isnan(resp(:,1)));
      if diff(params.maxlag)==0,
         % trim nans out of the stim/resp pair. this should
         % speed things up a bit
         stim=bstim(tokidx,:);
         resp=squeeze(resp(tokidx,rbinidx,:));
      else
         % don't trim anything off of the stim/resp pair. this
         % allows for valid temporal lags.
         stim=bstim;
         resp=squeeze(resp(:,rbinidx,:));
      end
      
      fprintf('stim is %d x %d. resp is %d x %d\n',size(stim),size(resp));
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % OLD SCHOOL RC, HACKED FROM CELLXCNODB
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % don't include baseline (all att) resp
      rsize=size(resp);
      respcount=rsize(2);
      attcount=size(resp,2);
      tbincount=1;
      
      firstseg=1;
      
      xccore;
      
      % moosh so that attention dimension is the right one to match
      % xcfit's expectations
      mH=reshape(mH,spacecount,tbincount,params.sfscount,...
                 1,attcount);
      eH=reshape(eH,spacecount,tbincount,params.sfscount,...
                 1,attcount);
      mSall=reshape(mSall,spacecount,tbincount,attcount);
      clear fdata strf tstrf
      
      params.iconside=iconside;
      if isfield(params,'altfit'),
         eval(params.altfit);
      else
         xcfit;
      end
      
      % output nl for each att kernel. float dc/g fit for each attidx
      linpred=zeros(size(resp,1),attcount);
      nlinpred=zeros(size(resp,1),attcount);
      expxc=zeros(attcount,1);
      fprintf('fit xc:');
      for attidx=1:attcount,
         
         % predict with local att kernel
         tstim=stim'-repmat(strf(attidx).mS,[1 size(stim,1)]);
         linpred(:,attidx)=kernpredict(strf(attidx).h,tstim,1,0);
         
         % find only response bins that are in appropriate att cond
         atokidx=find(~isnan(resp(:,attidx)));
         
         afitresp=resp(atokidx,attidx);
         epred=linpred(atokidx,attidx);
         dcgparms=fitdcgain(epred,afitresp);
         
         nlinpred(atokidx,attidx)=dcgain(dcgparms,epred);
         expxc(attidx)=xcov(nlinpred(atokidx,attidx),...
                            resp(atokidx,attidx),0,'coeff');
         fprintf(' %d: %.3f',attidx,expxc(attidx));
         
         strf(attidx).nltype='dcgain';
         strf(attidx).nlparms=dcgparms;
         strf(attidx).parms.attname=[cellid,' ',nlnames{attidx}];
      end
      fprintf('\n');
      
      % save important stuff from this bootidx
      vexpxc=cat(4,vexpxc,expxc);
      vstrf(:,:,:,bootidx)=strf;
      
      %
      % predict validatation data with each kernel
      %
      % save in:
      %valattpreds=zeros(resplen,attcount).*nan;
      
      disp('saving val preds for each att state');
      % ppidx contains ids of val data that were excluded on this boot
      
      % go through each att state and predict response to each stim
      for attidx=1:attcount,
         
         % get kernel
         tstrf=vstrf(attidx,1,1,bootidx);
         
         % do the linear prediction
         estim=bstim(ppidx,:)-repmat(tstrf.mS',length(ppidx),1);
         tpred=estim * tstrf.h;
         
         % output nl if specfied
         tnltype=tstrf.nltype;
         tnlparms=tstrf.nlparms;
         if ~isempty(tnltype) & ~strcmp(tnltype,'none'),
            valattpreds(ppidx,attidx)=feval(tnltype,tnlparms,tpred);
         else
            valattpreds(ppidx,attidx)=tpred;
         end
      end
      
   end % for bootidx
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. MEASURE PREDICTION ACCURACY FOR EACH COND + RANDOMIZED ATT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nStarting predictions...\n\n');

ESTIMATIONPHASE=0;
VALIDATIONPHASE=1;

% set up variables for saving xcs
xc=zeros(params.noisecount+1,attcount);
xcboot=zeros(params.noisecount+1,attcount,params.bootcount);

rprec=ones(blen,attcount).*nan;
rangle=ones(blen,attcount).*nan;

for noiseidx=1:params.noisecount+1,
   fprintf('noiseidx=%3d: ',noiseidx);
   
   resp=bresp;
   
   % assume there's only ONE response dimension!
   rbinidx=1;
   
   if noiseidx>1 & params.randcnf,
      % flip sign (around mean) for a random half of responses
      mm=nanmean(resp);
      for attidx=2:attcount,
         mm=nanmean(resp(:,rbinidx,attidx));
         resp(:,rbinidx,attidx)=resp(:,rbinidx,attidx) - mm;
         rr=round(rand(blen,1))*2-1;
         resp(:,rbinidx,attidx)=resp(:,rbinidx,attidx).*rr + mm;
      end
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % JACKKNIFE PREDICTIONS
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for attidx=1:attcount
      
      % set up vectors to hold preds from different jackknifes
      rpred=zeros(blen,1);
      rpred0=zeros(blen,1);
      
      for bootidx=1:params.bootcount;
         
         % figure out range to predict
         a0idx=vidx(round((bootidx-1)/params.bootcount*vcount+1):...
                    round(bootidx/params.bootcount*vcount));
         aidx=a0idx(find(~isnan(resp(a0idx,rbinidx,attidx))));
         
         % get kernel
         tstrf=vstrf(attidx,1,rbinidx,bootidx);
         tH=tstrf.h;
         tmS=tstrf.mS';
         tnltype=tstrf.nltype;
         tnlparms=tstrf.nlparms;
         attname=tstrf.parms.attname;
         
         % do the linear prediction
         estim=bstim(aidx,:)-repmat(tmS,length(aidx),1);
         tpred=estim * tH;
         
         if ~isempty(tnltype) & ~strcmp(tnltype,'none'),
            tpred=feval(tnltype,tnlparms,tpred);
         end
         
         % measure vector angle between each stimulus and the kernel
         estim0=sqrt(sum(estim.^2,2));
         estim0=estim./repmat(estim0,[1 spacecount]);
         if norm(tH)>0,
            tH0=tH./norm(tH);
         else
            tH0=zeros(size(tH));
         end
         rpred0(aidx)=acos(estim0 * tH0)./pi*180;
         
         rpred(aidx)=tpred;
         
         % maybe this helps? not really... rarely dip below 0
         %rpred(find(rpred<0))=0;
         
         % measure cc for this jackknife segment
         ract=resp(aidx,rbinidx,attidx);
         
         if length(ract)>0 & std(rpred(aidx))>0 & std(ract)>0,
            xcboot(noiseidx,attidx,bootidx)=xcov(rpred(aidx),ract,0,'coeff');
         end
      end
      
      % measure cc for whole attention condition
      aidx=find(~isnan(resp(:,rbinidx,attidx)));
      ract=resp(aidx,rbinidx,attidx);
      
      if noiseidx==1,
         % get expected xc for randomly shuffled pred... to
         % determine whether this is a reasonable pred or not
         [xc(noiseidx,attidx),randxc,tt,p]=...
             randxcov(rpred(aidx),ract,0,100);
      elseif length(ract)>0 & std(ract)>0 & std(rpred(anyokidx))>0,
         xc(noiseidx,attidx)=xcov(rpred(aidx),ract,0,'coeff');
      end
      
      % not doing full xc.  average over booted xc
      xc(noiseidx,attidx)=mean(xcboot(noiseidx,attidx,:));
      
      fprintf('%7.3f',xc(noiseidx,attidx));
      
      if noiseidx==1,
         % save true predictions
         rprec(:,attidx)=rpred;
      end
   end
   
   fprintf('\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. TEST EACH MODEL FOR SIGNIFICANTLY IMPROVED PREDICTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% place to record p values
pxc=zeros(size(xc,2),size(xc,3));
pxcboot=zeros(size(xcboot,2),size(xcboot,3));

% figure out p value for each model
for xx=1:size(xc,2),
   nn=length(find(xc(2:end,xx)>=xc(1,xx)));
   pxc(xx)=(nn+1)./(params.noisecount+1);
end

for xx=1:size(xcboot,2).*size(xcboot,3),
   nn=length(find(xcboot(2:end,xx)>=xcboot(1,xx)));
   pxcboot(xx)=(nn+1)./(params.noisecount+1);
end

fprintf('       p<     ');
fprintf(' %.3f',pxc);
fprintf('\n');

% assemble matrices for saving to sResults
predxc=cat(1,xc(1,:,:),mean(xc(2:end,:,:),1),std(xc(2:end,:,:),0,1));


% predict WHOLE response with kernels, rather than just diffs
%
% origresp(:,:,1)=(rAB+rab)/2
% origresp(:,:,2)=(rAB-rab)/2
%
% true rAB= origresp(:,:,1) + origresp(:,:,2)
%
% model: 
% given mr=(rAB+rab)/2 and dr=(rAa-rBb)/2:
% fit dr=mr * globalgain + globaldc
% then fit diff srf to dr'=dr-(mr*globalgain+globaldc)
%
% predict (no att):
% rAa = rprec(:,1)
% rBb = rprec(:,1)
%
% predict (dcg)
% rAa = rprec(:,1) + globalgain(2,bootidx)*rprec(:,1) + globaldc(2,bootidx)
% rBb = rprec(:,1) - globalgain(2,bootidx)*rprec(:,1) - globaldc(2,bootidx)
%
% predict (local)
% rAa= rprec(:,1) + globalgain(2,bootidx)*rprec(:,1)+globaldc(2,bootidx)+
%        rprec(:,2)
% rBb= rprec(:,1) - globalgain(2,bootidx)*rprec(:,1)-globaldc(2,bootidx)-
%        rprec(:,2)

% rAB/rab same deal for attidx=3 rather than attidx=2
% rA/rB  ?? sort of for 4, but not really.  hmm.
% ra/rb  should be same as for 4 but for 5.

respatt=zeros(blen*2,attcount);
prednoatt=zeros(blen*2,attcount);
preddc=zeros(blen*2,attcount);
preddcg=zeros(blen*2,attcount);
predlocal=zeros(blen*2,attcount);

fullxc=zeros(4,attcount-1);
randfxc=zeros(params.noisecount,4,attcount-1);

disp('calculating pred fracs');
keyboard
for attidx=1:attcount,
   
   prednoatt(1:blen,attidx)=rprec(:,1);
   prednoatt(blen+(1:blen),attidx)=rprec(:,1);
   
   if attidx==1,
      respatt(1:blen,attidx)=origresp(:,1);
      respatt(blen+(1:blen),attidx)=origresp(:,1);
      
      preddc(:,attidx)=prednoatt(:,attidx);
      preddcg(:,attidx)=prednoatt(:,attidx);
      predlocal(:,attidx)=prednoatt(:,attidx);
      
   elseif attidx<=3, % ie, Aa-Bb or AB-ab
      respatt(1:blen,attidx)=origresp(:,1)+origresp(:,attidx);
      respatt(blen+(1:blen),attidx)=origresp(:,1)-origresp(:,attidx);
      
      for bootidx=1:params.bootcount;
         ppidx=vidx((round((bootidx-1)/params.bootcount*vcount)+1):...
                    round(bootidx/params.bootcount*vcount));
         preddc(ppidx,attidx)=prednoatt(ppidx,attidx) + ...
             globaldconly(attidx,bootidx);
         preddc(blen+ppidx,attidx)=prednoatt(ppidx,attidx) - ...
             globaldconly(attidx,bootidx);
         preddcg(ppidx,attidx)=prednoatt(ppidx,attidx) + ...
             globaldc(attidx,bootidx) + ...
             rprec(ppidx,1)*globalgain(attidx,bootidx);
         preddcg(blen+ppidx,attidx)=prednoatt(ppidx,attidx) - ...
             globaldc(attidx,bootidx) - ...
             rprec(ppidx,1)*globalgain(attidx,bootidx);
      end
      
      predlocal(1:blen,attidx)=preddcg(1:blen,attidx) + rprec(:,attidx);
      predlocal(blen+(1:blen),attidx)=preddcg(blen+(1:blen),attidx) - ...
          rprec(:,attidx);
   else
      % ie, A-B or a-b
      if attidx==4,
         % ie, local attend in predictions
         prednoatt(1:blen,attidx)=predlocal(1:blen,3);
         prednoatt(blen+(1:blen),attidx)=predlocal(1:blen,3);
         
         % retrieve [rA; rB]
         % ie, rABab + rAB-ab + rA-B
         respatt(1:blen,attidx)=origresp(:,1)+origresp(:,3)+origresp(:,4);
         % ie, rABab + rAB-ab - rA-B
         respatt(blen+(1:blen),attidx)=origresp(:,1)+origresp(:,3)-origresp(:,4);
      else
         % ie, local attend out predictions
         prednoatt(1:blen,attidx)=predlocal(blen+(1:blen),3);
         prednoatt(blen+(1:blen),attidx)=predlocal(blen+(1:blen),3);
         % retrieve [ra; rb]
         respatt(1:blen,attidx)=origresp(:,1)-origresp(:,3)+origresp(:,4);
         respatt(blen+(1:blen),attidx)=origresp(:,1)-origresp(:,3)-origresp(:,4);
      end
      
      for bootidx=1:params.bootcount;
         ppidx=vidx((round((bootidx-1)/params.bootcount*vcount)+1):...
                    round(bootidx/params.bootcount*vcount));
         preddc(ppidx,attidx)=prednoatt(ppidx,attidx) + ...
             globaldconly(attidx,bootidx);
         preddc(blen+ppidx,attidx)=prednoatt(ppidx,attidx) - ...
             globaldconly(attidx,bootidx);
         preddcg(ppidx,attidx)=prednoatt(ppidx,attidx) + ...
             globaldc(attidx,bootidx) + ...
             rprec(ppidx,1)*globalgain(attidx,bootidx);
         preddcg(blen+ppidx,attidx)=prednoatt(ppidx,attidx) - ...
             globaldc(attidx,bootidx) - ...
             rprec(ppidx,1)*globalgain(attidx,bootidx);
      end
      
      predlocal(1:blen,attidx)=preddcg(1:blen,attidx) + rprec(:,attidx);
      predlocal(blen+(1:blen),attidx)=preddcg(blen+(1:blen),attidx) - ...
          rprec(:,attidx);
      
   end
   
   gidx=find(~isnan(predlocal(:,attidx)) & ~isnan(respatt(:,attidx)));
   fullxc(1,attidx)=xcov(prednoatt(gidx,attidx),respatt(gidx,attidx),0,'coeff');
   fullxc(2,attidx)=xcov(preddc(gidx,attidx),...
                         respatt(gidx,attidx),0,'coeff');
   fullxc(3,attidx)=xcov(preddcg(gidx,attidx),...
                         respatt(gidx,attidx),0,'coeff');
   fullxc(4,attidx)=xcov(predlocal(gidx,attidx),respatt(gidx,attidx),0,'coeff');
   
   if attidx>1,
      for noiseidx=1:params.noisecount,
         shuffidx=find(rand(blen,1)>0.5);
         rpreddc=predlocal(:,attidx);
         rpreddc(shuffidx)=predlocal(shuffidx+blen,attidx);
         rpreddc(shuffidx+blen)=predlocal(shuffidx,attidx);
         
         rpreddcg=predlocal(:,attidx);
         rpreddcg(shuffidx)=preddc(shuffidx,attidx) - ...
             preddc(shuffidx+blen,attidx)+predlocal(shuffidx+blen,attidx);
         rpreddcg(shuffidx+blen)=preddc(shuffidx+blen,attidx) - ...
             preddc(shuffidx,attidx)+predlocal(shuffidx,attidx);
         
         rpredlocal=predlocal(:,attidx);
         rpredlocal(shuffidx)=preddcg(shuffidx,attidx) - rprec(shuffidx,attidx);
         rpredlocal(shuffidx+blen)=preddcg(shuffidx+blen,attidx) + ...
             rprec(shuffidx,attidx);
         
         randfxc(noiseidx,2,attidx)=...
             xcov(rpreddc(gidx),respatt(gidx,attidx),0,'coeff');
         randfxc(noiseidx,3,attidx)=...
             xcov(rpreddcg(gidx),respatt(gidx,attidx),0,'coeff');
         randfxc(noiseidx,4,attidx)=...
             xcov(rpredlocal(gidx),respatt(gidx,attidx),0,'coeff');
      end
   end
end

mfxc=squeeze(mean(randfxc));
disp('testing dcg against no att (rather than dc-only)');
randfxc(:,3,:)=randfxc(:,2,:);
fullxc=[mfxc(2:end,:); fullxc(end,:)];
tfxc=reshape(fullxc,[1 size(fullxc)]);
fullp=squeeze(sum(repmat(tfxc,[params.noisecount+1 1]) <= ...
                  cat(1,tfxc,randfxc))./(params.noisecount+1));
clear randfxc tfxc mfxc



% clean up big matrices
mSA2=mean(sSA2,4);
clear sSA2 timage tstrf tstrf1 tsSA2 tfitres tr tp tokidx tmS threshparm




