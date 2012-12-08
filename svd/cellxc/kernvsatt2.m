% kernvsatt2.m: remake of kerncompN.m to make sure it works and to
% make the code/analysis a bit cleaner.
%
% created SVD 8/29/03 (hacked coarsely from kerncomp5.m)
% inspired by comments from lab meeting 8/27/03
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
% PENDING:
% load information about targets/patches?  (a la kerncompX)
%
% POSSIBLE ISSUES OR IMPROVEMENTS:
% selectivity change?
% relation to targets?
% contrast response thingy--better way think about it?
%
disp('kernvsatt2.m : re-revised attention tester');

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
   %if strcmp(cellid,'m0128'),
   %   disp('cellid=m0128: cropping in extra few pix around stim');
   %   params.stimloadparms{1}=14;
   %end
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
      tresp=zeros(size(oresp,1),size(oresp,2),length(btarglist)).*nan;
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

if length(acountgood)<length(accounts),
   bresp=bresp(:,:,acountgood);
   btarglist=btarglist(acountgood);
end
mpatches=[mean(fpatches,2) fpatches(:,btarglist(2:end))];
targlist=btarglist(2:end);

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

if 0 & params.batch==104,
   anyokidx=find(sum(~isnan(bresp(:,end,1:end)),3));
else
   anyokidx=find(sum(~isnan(bresp(:,end,2:end)),3));
end

% choose a random set of fixations for each bootstrap. this should
% avoid bias in the non-randomized noiseidx=1
vcount=length(anyokidx);

ss={};
for attidx=2:attcount,
   aokidx=find(~isnan(bresp(:,end,attidx)));
   [vv,ss{attidx}]=sort(rand(size(aokidx)));
   ss{attidx}=aokidx(ss{attidx});
end
ss{1}=shuffle(setdiff(anyokidx,cat(1,ss{:})));

vidx=[];
for bootidx=1:params.bootcount,
   for attidx=1:attcount,
      acount=length(ss{attidx});
      ppidx=(round((bootidx-1)/params.bootcount*acount)+1):...
            round(bootidx/params.bootcount*acount);
      
      vidx=[vidx; ss{attidx}(ppidx)];
   end
end

% generate random attention state sets
sn=zeros(length(anyokidx),params.noisecount+1);
for noiseidx=1:params.noisecount+1;
   [tt,sn(:,noiseidx)]=sort(rand(size(anyokidx)));
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

OUTNLMODE=4;
nlnames={'glob none','glob dc','glob g','glob dcg','all att'};
nlcount=length(nlnames);

% use for non-cheating local att preds
globalgain=zeros(attcount-1,params.bootcount);
globaldc=zeros(attcount-1,params.bootcount);
resplen=size(bresp,1);
valattpreds=zeros(resplen,attcount-1,nlcount);
brespres=ones(size(bresp)).*nan;

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
      
      params.minsfsidx=4;
      
      resp0=resp;
      attcount0=attcount;
      
      resp=resp(:,1);
      
      % remove mean (and variance?) before fitting baseline
      % kernel. this should help, right?
      fprintf('att norm: ');
      %rs=nanstd(resp0(:,1));
      for attidx=2:attcount0,
         aa=find(~isnan(resp0(:,attidx)));
         resp(aa)=resp(aa)-mean(resp0(aa,attidx));
         %if std(resp(aa))./rs>0,
         %   resp(aa)=resp(aa)./std(resp(aa)).*rs;
         %end
         %fprintf('%d: (resp-%.3f)/%.2f ',attidx-1,...
         %        mean(resp0(aa,attidx)),std(resp0(aa,attidx))./rs);
         fprintf('(%d: %.3f) ',attidx-1,mean(resp0(aa,attidx)));
      end
      fprintf('\n');
      
      attcount=1;
      firstseg=1;
      rsize=size(resp);
      respcount=rsize(2); % KLOODGE!  HOW TO DEAL?
      params.iconside=size(stim,2);
      
      if params.repexclude,
         xcitercore;
      else
         xccore;
      end
      
      tbincount=size(mH,2);
      
      % clear temp variables so that xcfit doesn't get confused
      % from last bootidx's temps
      clear fdata strf tstrf
      
      % by default, altfit='xcfit';
      eval(params.altfit);
      
      attcount=size(resp0,2);
      for nlidx=1:4,
         for attidx=1:attcount-1,
            savestrf(nlidx,attidx)=strf;
         end
      end
      
      % generate a no attention prediction
      tstim=stim'-repmat(strf.mS,[1 size(stim,1)]);
      noattpred=kernpredict(strf.h,tstim,1,0);
      
      % select time bins with valid data
      fitresp=resp0(:,1);
      tokidx=find(~isnan(noattpred) & ~isnan(fitresp));
      epred=noattpred(tokidx);
      fitresp=fitresp(tokidx);
      
      % matrix to save predictions in each att condition
      nlinpred=zeros(size(noattpred,1),5);
      
      % index matrix to label each time bin with appropriate attidx
      attcode=zeros(size(tokidx));
      filler=ones(size(tokidx));
      for fitidx=1:attcount-1,
         attcode(find(~isnan(resp0(tokidx,fitidx+1))))=fitidx;
      end
      attcode(attcode==0)=attcount;
      
      % nlidx=1: fit dc and gain to all att states
      [dcgparms,beta0]=fitdcgain(epred,fitresp);
      nlinpred(tokidx,1)=dcgain(dcgparms,epred);
      for fitidx=1:attcount-1,
         savestrf(1,fitidx).nltype='dcgain';
         savestrf(1,fitidx).nlparms=dcgparms;
      end
      
      % nlidx=2:  global H, dc att
      [dcgparms,beta0]=fitdcgain(epred,fitresp,[attcode filler]);
      %[dcgparms,beta0]=fitdcgaindumb(epred,fitresp,[attcode filler]);
      %keyboard
      nlinpred(tokidx,2)=dcgain(dcgparms,epred,[attcode filler+max(attcode)]);
      for fitidx=1:attcount-1,
         savestrf(2,fitidx).nltype='dcgain';
         savestrf(2,fitidx).nlparms=[dcgparms(fitidx);dcgparms(end)];
      end
      
      % nlidx=3:  global H, g att
      [dcgparms,beta0]=fitdcgain(epred,fitresp,[filler attcode]);
      %[dcgparms,beta0]=fitdcgaindumb(epred,fitresp,[filler attcode]);
      nlinpred(tokidx,3)=dcgain(dcgparms,epred,[filler attcode+1]);
      for fitidx=1:attcount-1,
         savestrf(3,fitidx).nltype='dcgain';
         savestrf(3,fitidx).nlparms=[dcgparms(1); dcgparms(fitidx+1)];
      end
      
      % nlidx=4:  global H, dcg att
      for fitidx=1:attcount-1,
         aokidx=find(~isnan(resp0(tokidx,fitidx+1)));
         
         [dcgparms,beta0]=fitdcgain(epred(aokidx),fitresp(aokidx));
         %[dcgparms,beta0]=fitdcgaindumb(epred(aokidx),fitresp(aokidx));
         
         nlinpred(tokidx(aokidx),4)=dcgain(dcgparms,epred(aokidx));
         savestrf(4,fitidx).nltype='dcgain';
         savestrf(4,fitidx).nlparms=dcgparms;
         
         % save for display
         globaldc(fitidx,bootidx)=dcgparms(1);
         globalgain(fitidx,bootidx)=dcgparms(2);
         
         fprintf('boot=%d fit=%d dc=%.3f gn=%.3f\n',...
                 bootidx,fitidx,dcgparms);
      end
      
      if params.oddout,
         
         % measure local attention effect for one state at a time
         % versus all others averaged
         resp1=resp0(:,2:end);
         for ridx=1:size(resp1,2),
            gidx=find(~isnan(resp1(:,ridx)));
            resp1(gidx,ridx)=resp1(gidx,ridx)-nlinpred(gidx,4);
         end
         
         for battidx=1:size(bresp,3)-1,
            fprintf('doing odd-out local att test oddidx=%d:\n',battidx);
            % don't include baseline (all att) resp
            resp=cat(2,resp1(:,battidx),...
                     nansum(resp1(:,[1:battidx-1 battidx+1:end])')');
            rsize=size(resp);
            respcount=rsize(2);
            attcount=rsize(2);
            
            params.sfsstepbak=params.sfsstep;
            params.sfsstep=params.sfsstep-1;
            firstseg=1;
            
            xccore;
            
            % moosh so that attention dimension is the right one to
            % match xcfit's expectations
            mH=reshape(mH,spacecount,tbincount,params.sfscount,...
                    1,attcount);
            eH=reshape(eH,spacecount,tbincount,params.sfscount,...
                       1,attcount);
            mSall=reshape(mSall,spacecount,tbincount,attcount);
            clear fdata strf tstrf
            
            % by default, altfit='xcfit';
            eval(params.altfit);

            params.sfsstep=params.sfsstepbak;
            
            disp('finished fitting residuals. now fit attentional gains and stuff');
            nlidx=4+battidx;
            
            % output nl for separate att kernels. float dc/g fit
            % for each attidx. this is the differential prediction
            % on top of the global dcg prediction
            linpred=zeros(size(resp,1),attcount-1);
            for attidx=1:attcount,
               
               % predict with local att kernel
               tstim=stim'-repmat(strf(attidx).mS,[1 size(stim,1)]);
               linpred(:,attidx)=kernpredict(strf(attidx).h,tstim,1,0);
               
               % find only response bins that are in appropriate
               % att cond 
               atokidx=find(~isnan(resp(:,attidx)));
               afitresp=resp(atokidx,attidx);
               epred=linpred(atokidx,attidx);
               dcgparms=fitdcgain(epred,afitresp);
               
               nlinpred(atokidx,nlidx)=dcgain(dcgparms,epred);
               
               strf(attidx).nltype='dcgain';
               strf(attidx).nlparms=dcgparms;
            end
            
            savestrf(nlidx,battidx)=strf(1);
            for attidx=[1:battidx-1 battidx+1:size(bresp,3)-1],
               savestrf(nlidx,attidx)=strf(2);
            end
            
            nlinpred(:,nlidx)=nlinpred(:,nlidx)+nlinpred(:,4);
            nlnames{nlidx}=sprintf('allodd%d',battidx);
         end
         
         attcount=size(bresp,3); % reset to deal with kludgy skip above
         nlcount=size(savestrf,1);
         
      else
         
         % if not params.oddout -- basically the standard now for
         % fitting local srfs
         
         % don't include baseline (all att) resp
         resp=resp0(:,2:end);
         rsize=size(resp);
         respcount=rsize(2);
         attcount=size(resp,2);
         for attidx=1:attcount,
            gidx=find(~isnan(resp(:,attidx)));
            
            if max(gidx)>size(nlinpred,1),
               keyboard
            end
            
            resp(gidx,attidx)=resp(gidx,attidx)-nlinpred(gidx,4);
         end
         
         params.sfsstepbak=params.sfsstep;
         if ~(isfield(params,'matchattsfs') & params.matchattsfs),
            params.sfsstep=params.sfsstep-1;
         end
         
         firstseg=1;
         if params.repexclude,
            xcitercore;
         else
            xccore;
         end
         
         % moosh so that attention dimension is the right one to match
         % xcfit's expectations
         mH=reshape(mH,spacecount,tbincount,params.sfscount,...
                    1,attcount);
         eH=reshape(eH,spacecount,tbincount,params.sfscount,...
                    1,attcount);
         mSall=reshape(mSall,spacecount,tbincount,attcount);
         clear fdata strf tstrf
         
         if isfield(params,'matchattsfs') & params.matchattsfs, 
            disp('matchattsfs for local to baseline strfs');
            sfsidx=savestrf(1).parms.sfsfit;
            sigidx=savestrf(1).parms.sigfit;
            
            for attidx=1:attcount,
               strf(1,attidx).h=shrinkage(mH(:,:,sfsidx,:,attidx),...
                                          eH(:,:,sfsidx,:,attidx),...
                                          sigrange(sigidx));
               strf(1,attidx).parms=savestrf(1).parms;
               strf(1,attidx).mS=savestrf(1).mS;
               strf(1,attidx).powunbiased=powunbiased(:,:,sfsidx);
               strf(1,attidx).nlparms=[];
               strf(1,attidx).nltype='none';
               strf(1,attidx).name=sprintf('%s matched attidx=%d',...
                                           cellid,attidx);
           end
         else
            % by default, altfit='xcfit';
            eval(params.altfit);
         end
         
         attcount=size(bresp,3); % reset to deal with kludgy skip above
         nlcount=length(nlnames);
         
         savestrf(5,:)=strf;
         params.sfsstep=params.sfsstepbak;
         
         disp('done fitting residuals. now fit local dc gain : (0,1)?');
         
         % nlidx=5: output nl for separate att kernels. float dc/g
         % fit for each attidx
         % this is actually just the differential prediction on top of
         % the global dcg prediction
         
         linpred=zeros(size(resp,1),attcount-1);
         for attidx=1:attcount-1,
            
            % predict with local att kernel
            tstim=stim'-repmat(strf(attidx).mS,[1 size(stim,1)]);
            linpred(:,attidx)=kernpredict(strf(attidx).h,tstim,1,0);
            
            % find only response bins that are in appropriate att cond
            atokidx=find(~isnan(resp(:,attidx)));
            
            %[cxy,exy,tt,p]=randxcov(linpred(atokidx,attidx),...
            %                        resp(atokidx,attidx),0,500);
            %[cxy exy p]
            
            afitresp=resp(atokidx,attidx);
            epred=linpred(atokidx,attidx);
            dcgparms=fitdcgain(epred,afitresp);
            
            nlinpred(atokidx,5)=dcgain(dcgparms,epred);
            
            savestrf(5,attidx).nltype='dcgain';
            savestrf(5,attidx).nlparms=dcgparms;
         end
         nlinpred(:,5)=nlinpred(:,5)+nlinpred(:,4);
      end
      
      % subtract off global pred and fit each attention kernel on
      % its own
      expxc=zeros(nlcount,1);
      for ii=1:nlcount,
         expxc(ii)=xcov(nlinpred(tokidx,ii),resp0(tokidx,1),0,'coeff');
         
         for attidx=1:attcount-1,
            savestrf(ii,attidx).parms.attname=nlnames{ii};
         end
      end
      
      % save important stuff from this bootidx
      vexpxc=cat(4,vexpxc,expxc);
      vstrf(:,:,:,bootidx)=savestrf;
      
      %
      % predict validatation data with each kernel
      %
      % save in:
      %valattpreds=zeros(resplen,attcount,nlusecount);
      
      disp('saving val preds for each att state');
      a0idx=vidx((round((bootidx-1)/params.bootcount*vcount)+1):...
                 round(bootidx/params.bootcount*vcount));
      
      % go through each att state and model
      for nlidx=1:nlcount,
         for attidx=1:attcount-1,
            
            % get kernel
            tstrf=vstrf(nlidx,attidx,1,bootidx);
            tH=tstrf.h;
            tmS=tstrf.mS';
            tnltype=tstrf.nltype;
            tnlparms=tstrf.nlparms;
            attname=tstrf.parms.attname;
            
            % do the linear prediction
            estim=bstim(a0idx,:)-repmat(tmS,length(a0idx),1);
            tpred=estim * tH;
            
            if ~isempty(tnltype) & ~strcmp(tnltype,'none'),
               valattpreds(a0idx,attidx,nlidx)=feval(tnltype,tnlparms,tpred);
            else
               valattpreds(a0idx,attidx,nlidx)=tpred;
            end
         end
      end
      
      ADD_DCG=1;
      if ADD_DCG,
         % add dcg pred onto local residual pred
         for nlidx=5:nlcount,
            for att2=1:attcount-1,
               aidx=a0idx(find(~isnan(bresp(a0idx,1,att2+1))));
               valattpreds(aidx,:,nlidx)=valattpreds(aidx,:,nlidx)+...
                   repmat(valattpreds(aidx,att2,4),[1 attcount-1]);
            end
         end
      else
         disp('not adding dcg pred onto local pred');
         % subtract dcg pred from observed resp
         nlidx=4;
         for att2=1:attcount-1,
            aidx=a0idx(find(~isnan(bresp(a0idx,1,att2+1))));
            brespres(aidx,1,att2+1)=...
                bresp(aidx,1,att2+1)-valattpreds(aidx,att2,4);
            brespres(aidx,1,1)=brespres(aidx,1,att2+1);
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
nlcount=size(vstrf,1);

% force anyokidx to not include invalid att states, even though
% they might have been used for computing the baseline STRFs
anyokidx=find(sum(~isnan(bresp(:,end,2:end)),3));

% re-generate random attention state sets
sn=zeros(length(anyokidx),params.noisecount+1);
for noiseidx=1:params.noisecount+1;
   [tt,sn(:,noiseidx)]=sort(rand(size(anyokidx)));
end

% pred xc for pairs of att states
paircount=(attcount-1)*(attcount-2)/2;
pairidx=zeros(paircount,2);
pidx=0;
for a1idx=1:attcount-1,
   for a2idx=a1idx+1:attcount-1,
      pidx=pidx+1;
      pairidx(pidx,1)=a1idx;
      pairidx(pidx,2)=a2idx;
   end
end

xc=zeros(params.noisecount+1,nlcount);
xccross=zeros(params.noisecount+1,nlcount,attcount);
xcboot=zeros(params.noisecount+1,nlcount,params.bootcount);
xcpair=zeros(params.noisecount+1,nlcount,paircount);

attshuff=zeros(attcount-1,params.bootcount,params.noisecount+1);
rprec=ones(blen,nlcount).*nan;
rangle=ones(blen,nlcount).*nan;

for noiseidx=1:params.noisecount+1,
   fprintf('noiseidx=%3d: ',noiseidx);
   
   resp=ones(size(bresp))*nan;
   resp(anyokidx,:,1)=bresp(anyokidx,:,1);
   if noiseidx==1 | ~params.randcnf,
      resp(:,:,2:end)=bresp(:,:,2:end);
   else
      % pick random atts for each fixation. effectively the same as 0?
      for attidx=2:attcount,
         nn=anyokidx(sn(cumncount(attidx-1)+1:cumncount(attidx),noiseidx));
         resp(nn,:,attidx)=bresp(nn,:,1);
      end
   end
   
   % assume there's only ONE response dimension!
   rbinidx=1;
   
   for nlidx=1:nlcount,
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % JACKKNIFE PREDICTIONS
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % set up vectors to hold preds from different models
      rpred=zeros(blen,1);
      rpred0=zeros(blen,1);
      
      for attidx=2:attcount
         
         % figure out range to predict
         aidx=find(~isnan(resp(:,rbinidx,attidx)));
         
         rpred(aidx)=valattpreds(aidx,attidx-1,nlidx);
         %rpred(find(rpred<0))=0;
      end
      
      for attidx=2:attcount,
         % test single attcond preds only in correct att cond
         aidx=find(~isnan(bresp(:,rbinidx,attidx)));
         
         if ADD_DCG | nlidx<=4,
            ract=bresp(aidx,rbinidx,1);
         else
            ract=brespres(aidx,rbinidx,1);
         end
         
         % within each attention state, cc will be biased toward
         % higher values because dc/gain changes won't matter
         if length(ract)>0 & std(ract)>0,
            xccross(noiseidx,nlidx,attidx)=xcov(rpred(aidx),ract,0,'coeff');
         else
            xccross(noiseidx,nlidx,attidx)=0;
         end
      end
      
      for bootidx=1:params.bootcount
         a0idx=vidx((round((bootidx-1)/params.bootcount*vcount)+1):...
                    round(bootidx/params.bootcount*vcount));
         if ADD_DCG | nlidx<=4,
            ract=bresp(a0idx,rbinidx,1);
         else
            ract=brespres(a0idx,rbinidx,1);
         end
         
         if std(rpred(a0idx))>0 & std(ract)>0,
            xcboot(noiseidx,nlidx,bootidx)=xcov(rpred(a0idx),ract,0,'coeff');
         end
      end
      
      for pidx=1:paircount,
         tokidx=find(~isnan(bresp(:,1,pairidx(pidx,1)+1)) | ...
                     ~isnan(bresp(:,1,pairidx(pidx,2)+1)));
         if ADD_DCG | nlidx<=4,
            ract=bresp(tokidx,rbinidx,1);
         else
            ract=brespres(tokidx,rbinidx,1);
         end
         if length(ract)>0 & std(ract)>0,
            xcpair(noiseidx,nlidx,pidx)=xcov(rpred(tokidx),ract,0,'coeff');
         end
      end
      
      if ADD_DCG | nlidx<=4,
         ract=bresp(anyokidx,rbinidx,1);
      else
         ract=brespres(anyokidx,rbinidx,1);
      end
      
      if noiseidx==1 & nlidx==1,
         % get expected xc for randomly shuffled pred... to
         % determine whether this is a reasonable pred or not
         [xc(noiseidx,nlidx),randxc,tt,p]=...
             randxcov(rpred(anyokidx),ract,0,100);
      elseif length(ract)>0 & std(ract)>0 & std(rpred(anyokidx))>0,
         xc(noiseidx,nlidx)=xcov(rpred(anyokidx),ract,0,'coeff');
      end
      
      % not doing full xc.  average over booted xc
      %xc(noiseidx,attidx)=mean(xcboot(noiseidx,attidx,:));
      
      fprintf(' %.3f',xc(noiseidx,nlidx));
      
      if noiseidx==1,
         rprec(:,nlidx)=rpred;
      end
   end
   
   fprintf('\n');
end



if 0

realattidx=zeros(size(bresp,1),1);
for attidx=2:attcount,
   tokidx=find(~isnan(bresp(:,1,attidx)));
   realattidx(tokidx)=attidx-1;
end

for noiseidx=1:params.noisecount+1,
   fprintf('noiseidx=%3d: ',noiseidx);
   
   resp=ones(size(bresp))*nan;
   resp(anyokidx,:,1)=bresp(anyokidx,:,1);
   if noiseidx==1 | ~params.randcnf,
      resp(:,:,2:end)=bresp(:,:,2:end);
   else
      % pick random atts for each fixation. effectively the same as 0?
      for attidx=2:attcount,
         nn=anyokidx(sn(cumncount(attidx-1)+1:cumncount(attidx),noiseidx));
         resp(nn,:,attidx)=bresp(nn,:,1);
      end
   end
   
   % assume there's only ONE response dimension!
   rbinidx=1;
   
   counter=zeros(blen,nlcount);
   for nlidx=1:nlcount,
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % JACKKNIFE PREDICTIONS
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      % set up vectors to hold preds from different models
      rpred=zeros(blen,1);
      rpred0=zeros(blen,1);
      
      % generate pred for each bootstrapped segment
      for bootidx=1:params.bootcount;
         a0idx=vidx((round((bootidx-1)/params.bootcount*vcount)+1):...
                    round(bootidx/params.bootcount*vcount));
         
         % pull out stimulus
         for attidx=2:attcount
            
            % figure out range to predict
            aidx=a0idx(find(~isnan(resp(a0idx,rbinidx,attidx))));
            counter(aidx,nlidx)=counter(aidx,nlidx)+1;
            
            % get kernel
            tstrf=vstrf(nlidx,attidx-1,rbinidx,bootidx);
            tH=tstrf.h;
            tmS=tstrf.mS';
            tnltype=tstrf.nltype;
            tnlparms=tstrf.nlparms;
            attname=tstrf.parms.attname;
            
            % do the linear prediction
            estim=bstim(aidx,:)-repmat(tmS,length(aidx),1);
            tpred=estim * tH;
            
            if OUTNLMODE==4,
               % fix non-local preds to be the same
               if strcmp(attname,'all att'),
                  if ~isempty(tnltype) & ~strcmp(tnltype,'none'),
                     tpred=feval(tnltype,tnlparms,tpred);
                  end
                  tpred=tpred+rprec(aidx,nlidx-1);
               elseif strcmp(attname,'glob dcg'),
                  for attidx00=1:attcount-1,
                     ttt=find(realattidx(aidx)==attidx00);
                     nlparms00=[vstrf(nlidx,attidx00,rbinidx,...
                                   bootidx).nlparms(1) tnlparms(2)];
                     tpred(ttt)=feval(tnltype,nlparms00,tpred(ttt));
                  end
               else
                  if ~isempty(tnltype) & ~strcmp(tnltype,'none'),
                     tpred=feval(tnltype,tnlparms,tpred);
                  end
               end
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
            
            if 0,
               figure(3);
               plot(rpred(anyokidx),'g')
               figure(4);
               plot(tstrf.h);
               [bootidx attidx tstrf.nlparms']
               pause
            end
            
            rpred(aidx)=tpred;
         end
         
         % rectify... this seems like an ok thing to do. since
         % actual responses are always >=0 and this is unbiased
         %rpred(find(rpred<0))=0;
         
         
         ract=bresp(a0idx,rbinidx,1);
         if std(rpred(a0idx))>0 & std(ract)>0,
            xcboot(noiseidx,nlidx,bootidx)=xcov(rpred(a0idx),ract,0,'coeff');
         end
      end
      
      % evaluate against actual resp
      for attidx=1:attcount,
         tokidx=find(~isnan(bresp(:,1,attidx)));
         ract=bresp(tokidx,rbinidx,1);
         
         % within each attention state, cc will be biased toward
         % higher values because dc/gain changes won't matter
         if length(ract)>0 & std(ract)>0,
            xccross(noiseidx,nlidx,attidx)=xcov(rpred(tokidx),ract,0,'coeff');
         else
            xccross(noiseidx,nlidx,attidx)=0;
         end
      end
      
      for pidx=1:paircount,
         tokidx=find(~isnan(bresp(:,1,pairidx(pidx,1)+1)) | ...
                     ~isnan(bresp(:,1,pairidx(pidx,2)+1)));
         ract=bresp(tokidx,rbinidx,1);
         if length(ract)>0 & std(ract)>0,
            xcpair(noiseidx,nlidx,pidx)=xcov(rpred(tokidx),ract,0,'coeff');
         end
      end
      
      ract=bresp(anyokidx,rbinidx,1);
      if noiseidx==1 & nlidx==1,
         % get expected xc for randomly shuffled pred... to
         % determine whether this is a reasonable pred or not
         [xc(noiseidx,nlidx),randxc,tt,p]=...
             randxcov(rpred(anyokidx),ract,0,100);
      elseif length(ract)>0 & std(ract)>0 & std(rpred(anyokidx))>0,
         xc(noiseidx,nlidx)=xcov(rpred(anyokidx),ract,0,'coeff');
      end
      
      % not doing full xc.  average over booted xc
      xc(noiseidx,attidx)=mean(xcboot(noiseidx,attidx,:));
      
      fprintf(' %.3f',xc(noiseidx,nlidx));
      
      if noiseidx==1,
         rprec(:,nlidx)=rpred;
         rangle(:,nlidx)=rpred0;
      end
   end
   
   fprintf('\n');
   
   % update queue if active
   dbsetqueue;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. TEST EACH MODEL FOR SIGNIFICANTLY IMPROVED PREDICTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if OUTNLMODE==3,
   % make real pred the avg of control preds for next att cond up
   xc(1,4)=mean(xc(2:end,5));
   xc(1,5)=mean(xc(2:end,6));
end

% place to record p values
pxc=zeros(size(xc,2),size(xc,3));
pxccross=zeros(size(xccross,2),size(xccross,3));
pxcboot=zeros(size(xcboot,2),size(xcboot,3));
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

for xx=1:size(xccross,2).*size(xccross,3),
   % flipped sign of inequality because this is actually MSE -- NOT xc
   nn=length(find(xccross(2:end,xx)<=xccross(1,xx)));
   pxccross(xx)=(nn+1)./(params.noisecount+1);
end
for xx=1:size(xcboot,2).*size(xcboot,3),
   nn=length(find(xcboot(2:end,xx)>=xcboot(1,xx)));
   pxcboot(xx)=(nn+1)./(params.noisecount+1);
end
for xx=1:size(xcpair,2).*size(xcpair,3),
   nn=length(find(xcpair(2:end,xx)>=xcpair(1,xx)));
   pxcpair(xx)=(nn+1)./(params.noisecount+1);
end

fprintf('       p<     ');
fprintf(' %.3f',pxc);
fprintf('\n');

% assemble matrices for saving to sResults
predxc=cat(1,xc(1,:,:),mean(xc(2:end,:,:),1),std(xc(2:end,:,:),0,1));
predxccross=cat(1,xccross(1,:,:),mean(xccross(2:end,:,:),1),...
                std(xccross(2:end,:,:),0,1));
predxcpair=cat(1,xcpair(1,:,:),mean(xcpair(2:end,:,:),1),...
                std(xcpair(2:end,:,:),0,1));


fprintf('fitting sigmoids:\nattidx=');
njack=params.bootcount;
asigparms=zeros(4,njack,attcount);
for attidx=1:attcount,
   fprintf('%d ',attidx);
   
   if attidx==1,
      tokidx=find(sum(~isnan(bresp(:,1,2:end)),3)>0);
   else
      tokidx=find(~isnan(bresp(:,1,attidx)));
   end
   
   jackstep=length(tokidx)/njack;
   
   for jj=1:njack,
      % use same jackknife sets as for fitting 
      useidx=vidx([1:round((jj-1)/params.bootcount*vcount) ...
                 round(jj/params.bootcount*vcount+1):end]);
      useidx=useidx(find(~isnan(bresp(useidx,1,attidx))));
      ract=bresp(useidx,1,1);
      rpred=rprec(useidx,1);
      
      asigparms(:,jj,attidx)=fitsigmoid(rpred,ract,0);
   end
end
fprintf('\n');
sigparms=squeeze(mean(asigparms,2));
esigparms=squeeze(std(asigparms,1,2)).*sqrt(njack-1);
nsigparms=sigparms./repmat(sigparms(:,1),[1 attcount]);
nesigparms=esigparms./repmat(sigparms(:,1),[1 attcount]);

disp('using RJP''s single trial noise analysis:');
winsize=params.respfilterparms{2}-params.respfilterparms{1};
resphz=round(bresp(:,1).*winsize);
[mu,alpha,beta]=reversepoisson(resphz);
rmax=singletrialceiling(resphz,alpha,beta);
predinf=sqrt(predxc(1:2,:).^2./rmax.^2);

disp('saving the right things?');

% clean up big matrices
mSA2=mean(sSA2,4);
clear sSA2 timage tstrf tstrf1 tsSA2 tfitres tr tp tokidx tmS threshparm


