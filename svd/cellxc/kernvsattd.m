% kernvsattd.m: remake of kerncompN.m to work for differential
% kernels... controls for problem of different residual noise.
%
% created SVD 9/18/04 (hacked coarsely from kernvsatt2.m)
%
% BASIC FLOW:
% 1. load relevant stimuli and resposnes (these can be jamie's fvvs
%    data or ben's dms data)
%
% 2. transform into eigenvector domain to speed kernel estimation.
% 3. break into jackknife sets
% 4. fit various models:
%     a. linear filters for no attention & tuning shift
%        (use double-jack analysis to get nice clean kernels
%        without cheating and using validation data)
%     b. apply same bias correction to each strf?
%        arbitrary? how to choose regularization parameter???
%     c. fit output nls, with and without attention floating:
%        (dc, gain, threshold...or some subset??)
% 5. predict validation data for each jackknife and randomized
%    attention
% 6. test for significant improvements with each model component.
%
%
disp('kernvsattd.m : differential attention tester');

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
params.oddout=getparm(params,'oddout',0);
params.altfit=getparm(params,'altfit','xcfit');

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
acountgood=find(accounts>0);
bresp=bresp(:,:,acountgood);

mpatches=mpatches(:,acountgood);
btarglist=btarglist(acountgood);

clear z

% figure out size of response
blen=size(bresp,1);
respcount=size(bresp,2);
attcount=size(bresp,3);
spacecount=size(bstim,2);

if ~params.meansub,
   disp('adding on dc channel');
   bstim(:,spacecount+1)=1;
   spacecount=spacecount+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. BREAK INTO JACKKNIFE SETS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find all time bins with valid responses in an attention state > 1
% (the first set is ALL attention ... throw away data not in the
% other categories)
%
aokidx=find(sum(~isnan(bresp(:,end,2:end)),3));

% figure out number of valid fixations in each attentional condition
ncount=squeeze(sum(~isnan(bresp(:,end,:)),1));
cumncount=[0;cumsum(ncount(2:end))];
anyokidx=find(sum(~isnan(bresp(:,end,2:end)),3));

% choose a random set of fixations for each bootstrap. this should
% avoid bias in the non-randomized noiseidx=1
vcount=length(anyokidx);

ss={};
for attidx=2:attcount,
   aokidx=find(~isnan(bresp(:,end,attidx)));
   [vv,ss{attidx-1}]=sort(rand(size(aokidx)));
   ss{attidx-1}=aokidx(ss{attidx-1});
end

vidx=[];
for bootidx=1:params.bootcount,
   for attidx=2:attcount,
      acount=length(ss{attidx-1});
      ppidx=(round((bootidx-1)/params.bootcount*acount)+1):...
            round(bootidx/params.bootcount*acount);
      
      vidx=[vidx; ss{attidx-1}(ppidx)];
   end
end

% generate random attention state sets
sn=zeros(length(anyokidx),params.noisecount+1);
for noiseidx=1:params.noisecount+1;
   [tt,sn(:,noiseidx)]=sort(rand(size(anyokidx)));
end

% pred xc for pairs of att states
paircount=(attcount-1)*(attcount-2)/2;
pairidx=zeros(paircount,2);
pidx=0;
pairnames={};
pairresp=zeros(size(bresp,1),size(bresp,2),paircount);
pairdresp=zeros(size(bresp,1),size(bresp,2),paircount);
for a1idx=1:attcount-1,
   for a2idx=a1idx+1:attcount-1,
      pidx=pidx+1;
      pairidx(pidx,1)=a1idx;
      pairidx(pidx,2)=a2idx;
      pairnames{pidx}=sprintf('%d v %d',a1idx,a2idx);
   end
end

if 0
ll1=sum(~isnan(bresp(:,1,a1idx)));
ll2=sum(~isnan(bresp(:,1,a2idx)));
ll1=ll1./(ll1+ll2);
ll2=1-ll1;
tpair=permute(cat(3,bresp(:,:,a1idx),bresp(:,:,a2idx)),[3 1 2]);
pairresp(:,:,pidx)=nanmean(tpair)';

tpair=permute(cat(3,bresp(:,:,a1idx),-bresp(:,:,a2idx)),[3 1 2]);
pairdresp(:,:,pidx)=nanmean(tpair)';
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

OUTNLMODE=6;
nlnames={'mean','diff'};
nlcount=length(nlnames);

% use for non-cheating local att preds
globalgain=zeros(attcount-1,params.bootcount);
globaldc=zeros(attcount-1,params.bootcount);
resplen=size(bresp,1);
valattpreds=zeros(resplen,paircount,nlcount);
valattdc=zeros(resplen,paircount,nlcount);

% cross correlate and apply bias correction, leaving in eig domain
for noiseidx=1:kerncount,
   fprintf('\nnoiseidx=%d: ',noiseidx);
   
   % resp is temporary matrix with attention assigned according to
   % noise (for randomized attention kernels)
   
   eresp=ones(size(bresp))*nan;
   eresp(anyokidx,:,1)=bresp(anyokidx,:,1);
   % assign resp values from brespnoglob (ie, with dc and gain removed
   % from specific attention state responses)
   if noiseidx==1,
      eresp(:,:,2:end)=bresp(:,:,2:end);
   else
      % pick random atts
      for attidx=2:attcount,
         nn=anyokidx(sn(cumncount(attidx-1)+1:cumncount(attidx),noiseidx));
         eresp(nn,:,attidx)=bresp(nn,:,1);
      end
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
      
      resp0=eresp(:,rbinidx,:);
      resp0(ppidx,:,:)=nan;
      
      params.minsfsidx=1;
      attcount0=attcount;
      
      for pidx=1:paircount,
         fprintf('PAIR : %s\n',pairnames{pidx});
         
         tpair=cat(2,resp0(:,pairidx(pidx,1)+1),resp0(:,pairidx(pidx,2)+1));
         dpair=cat(2,resp0(:,pairidx(pidx,1)+1),-resp0(:,pairidx(pidx,2)+1));
         resp=[nanmean(tpair')',nanmean(dpair')'];
         
         tokidx=find(~isnan(resp(:,1)));
         if diff(params.maxlag)==0,
            % trim nans out of the stim/resp pair. this should
            % speed things up a bit
            stim=bstim(tokidx,:);
            resp=squeeze(resp(tokidx,:));
            
            p1idx=find(~isnan(resp0(tokidx,pairidx(pidx,1)+1)));
            p2idx=find(~isnan(resp0(tokidx,pairidx(pidx,2)+1)));
         else
            % don't trim anything off of the stim/resp pair. this
            % allows for valid temporal lags.
            stim=bstim;
            
            p1idx=find(~isnan(resp0(:,pairidx(pidx,1)+1)));
            p2idx=find(~isnan(resp0(:,pairidx(pidx,2)+1)));
         end
         
         fprintf('stim is %d x %d. resp is %d x %d\n',size(stim),size(resp));
         
         % scale responses from each attention condition so that
         % they contribute equally to kernel estimates. ie, if
         % there are fewer samples, weigh them more.  this should
         % produce appropriate baseline h_0 for mean, and h_d for
         % the difference between pairs.
         %
         ll1=(length(p1idx)+length(p2idx))./(length(p1idx).*2);
         ll2=(length(p1idx)+length(p2idx))./(length(p2idx).*2);
         
         % (by scaling both response and stimulus, actual gains do
         % not change, only the amount that samples from each
         % condition are weighed in the final estimate)
         %
         resp(p1idx,:)=resp(p1idx,:).*ll1;
         resp(p2idx,:)=resp(p2idx,:).*ll2;
         stim(p1idx,:)=stim(p1idx,:).*ll1;
         stim(p2idx,:)=stim(p2idx,:).*ll2;
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % OLD SCHOOL RC, HACKED FROM CELLXCNODB
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         attcount=2;
         firstseg=1;
         rsize=size(resp);
         respcount=rsize(2);
         params.iconside=size(stim,2);
         
         if params.repexclude,
            xcitercore;
         else
            xccore;
         end
         
         % moosh so that attention dimension is the right one to match
         % xcfit's expectations
         tbincount=size(mH,2);
         mH=reshape(mH,spacecount,tbincount,params.sfscount,...
                    1,attcount);
         eH=reshape(eH,spacecount,tbincount,params.sfscount,...
                    1,attcount);
         mSall=reshape(mSall,spacecount,tbincount,attcount);
         
         % clear temp variables so that xcfit doesn't get confused
         clear fdata strf tstrf
         
         % by default, altfit='xcfit';
         eval(params.altfit);
         
         % fit dc and gain to all each pair mean and diff
         
         tstim=stim'-repmat(strf(1).mS,[1 size(stim,1)]);
         linpred=zeros(size(resp));
         nlinpred=zeros(size(resp));
         expxc=zeros(length(nlnames),1);
         for nlidx=1:length(nlnames),
            
            % predic with linear filter
            linpred(:,nlidx)=kernpredict(strf(nlidx).h,tstim,1,0);
            
            fitresp=resp(:,nlidx);
            
            tokidx=find(~isnan(linpred(:,nlidx)) & ~isnan(resp(:,nlidx)));
            epred=linpred(tokidx,nlidx);
            fitresp=resp(tokidx,nlidx);
            
            % actually find dc and gain parameters. should gain
            % always be 1?
            [dcgparms,beta0]=fitdcgain(epred,fitresp);
            nlinpred(tokidx,nlidx)=dcgain(dcgparms,epred);
            strf(nlidx).nltype='dcgain';
            strf(nlidx).nlparms=dcgparms;
            strf(nlidx).parms.attname=[nlnames{nlidx},' - ',pairnames{pidx}];
            
            % save strf for posterity
            savestrf(nlidx,pidx)=strf(nlidx);
            
            expxc(nlidx)=xcov(nlinpred(tokidx,nlidx),fitresp,0,'coeff');
            
            fprintf('nl=%d p=%d dc=%.3f g=%.3f\n',nlidx,pidx,dcgparms);
         end
         
         % save important stuff from this pairidx,bootidx
         vexpxc(:,pidx,:,bootidx)=expxc;
         vstrf(:,pidx,:,bootidx)=savestrf(:,pidx);
      end
      
      %
      % predict validatation data with each kernel
      %
      % save in:
      %valattpreds=zeros(resplen,attcount,nlusecount);
      
      disp('saving val preds for each att state');
      a0idx=vidx((round((bootidx-1)/params.bootcount*vcount)+1):...
                 round(bootidx/params.bootcount*vcount));
      
      % go through each att state and model
      for pidx=1:paircount,
         for nlidx=1:length(nlnames),
            % get kernel
            tstrf=vstrf(nlidx,pidx,1,bootidx);
            tH=tstrf.h;
            tmS=tstrf.mS';
            tnltype=tstrf.nltype;
            tnlparms=tstrf.nlparms;
            attname=tstrf.parms.attname;
            
            % do the linear prediction
            estim=bstim(a0idx,:)-repmat(tmS,length(a0idx),1);
            tpred=estim * tH;
            
            if ~isempty(tnltype) & ~strcmp(tnltype,'none'),
               valattdc(a0idx,pidx,nlidx)=tnlparms(1);
               tnlparms(1)=0;
               valattpreds(a0idx,pidx,nlidx)=feval(tnltype,tnlparms,tpred);
            else
               valattpreds(a0idx,pidx,nlidx)=tpred;
            end
         end
      end
      
   end % for bootidx
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. TEST PREDICTIONS FOR EACH MODEL + RANDOMIZED ATT PREDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nStarting predictions...\n\n');

ESTIMATIONPHASE=0;
VALIDATIONPHASE=1;

% set up variables for saving xcs
nlcount=size(vstrf,1);   % ==2, mean and diff
attcount=size(vstrf,2);  % ==paircount

xc=zeros(params.noisecount+1,nlcount);
xccross=zeros(params.noisecount+1,nlcount,attcount);
xcboot=zeros(params.noisecount+1,nlcount,params.bootcount);
xcpair=zeros(params.noisecount+1,nlcount,paircount);

attshuff=zeros(attcount-1,params.bootcount,params.noisecount+1);
rprec=ones(blen,nlcount).*nan;
rangle=ones(blen,nlcount).*nan;

fprintf('pair:         ');
fprintf('%12s',pairnames{:});
fprintf('\n');
for noiseidx=1:params.noisecount+1,
   fprintf('noiseidx=%3d: ',noiseidx);
   
   % assume there's only ONE response dimension!
   rbinidx=1;
   
   for pidx=1:paircount,
      
       % find responses matching either state in this pair
      aokidx=find(~isnan(bresp(:,1,pairidx(pidx,1)+1)) | ...
                  ~isnan(bresp(:,1,pairidx(pidx,2)+1)));
      
      rpred=valattpreds(aokidx,pidx,1)+valattdc(aokidx,pidx,1);
      ract=bresp(aokidx,rbinidx,1);
      
      if noiseidx==1 & pidx==1,
         % get expected xc for randomly shuffled pred... to
         % determine whether this is a reasonable pred or not
         [xcpair(noiseidx,1,pidx),randxc,tt,p]=...
             randxcov(rpred,ract,0,200);
      elseif length(ract)>0 & std(ract)>0 & std(rpred)>0,
         xcpair(noiseidx,1,pidx)=xcov(rpred,ract,0,'coeff');
      end
      
      flips=~isnan(bresp(aokidx,1,pairidx(pidx,1)+1)) - ...
            ~isnan(bresp(aokidx,1,pairidx(pidx,2)+1));
      
      if noiseidx==1,
         rflips=flips;
      else
         % shuffle flips on randomized attention
         rflips=shuffle(flips);
      end
      
      rpred=rpred+valattdc(aokidx,pidx,2).*flips + ...
            valattpreds(aokidx,pidx,2).*rflips;
      xcpair(noiseidx,2,pidx)=xcov(rpred,ract,0,'coeff');
      
      fprintf('%6.3f',xcpair(noiseidx,:,pidx));
   end
   
   xc(noiseidx,:)=mean(xcpair(noiseidx,:,:),3);
   
   fprintf('\n');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. TEST EACH MODEL FOR SIGNIFICANTLY IMPROVED PREDICTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% place to record p values
pxc=zeros(size(xc,2),size(xc,3));
pxcpair=zeros(size(xcpair,2),size(xcpair,3));

% figure out p value for each model
for xx=1:size(xc,2),
   % add random amount to avoid p=1.0 for all non-modulated cells
   % this is a kludge, but it's better than heavy clumps at p<1.0
   if xx>1,
      rr=(rand(size(xc,1),1)-0.5).*(1000*eps);
   else
      rr=zeros(size(xc,1),1);
   end
   nn=length(find(xc(2:end,xx)+rr(2:end)>=xc(1,xx)+rr(1)));
   pxc(xx)=(nn+1)./(params.noisecount+1);
end

for xx=1:size(xcpair,2).*size(xcpair,3),
   nn=length(find(xcpair(2:end,xx)>=xcpair(1,xx)));
   pxcpair(xx)=(nn+1)./(params.noisecount+1);
end

fprintf('       p<    ');
fprintf('%6.3f',pxcpair(:));
fprintf('\n');

% assemble matrices for saving to sResults
predxc=cat(1,xc(1,:,:),mean(xc(2:end,:,:),1),std(xc(2:end,:,:),0,1));
predxcpair=cat(1,xcpair(1,:,:),mean(xcpair(2:end,:,:),1),...
                std(xcpair(2:end,:,:),0,1));


disp('using RJP''s single trial noise analysis:');
winsize=params.respfilterparms{2}-params.respfilterparms{1};
resphz=round(bresp(:,1).*winsize);
[mu,alpha,beta]=reversepoisson(resphz);
rmax=singletrialceiling(resphz,alpha,beta);
predinf=sqrt(predxc(1:2,:).^2./rmax.^2);
predpairinf=sqrt(predxcpair(1:2,:,:).^2./rmax.^2);

disp('saving the right things?');

% clean up big matrices
mSA2=mean(sSA2,4);
clear sSA2 timage tstrf tstrf1 tsSA2 tfitres tr tp tokidx tmS ...
   threshparm sn


