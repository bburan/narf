% kerncomp4.m : compare kernels from different attentional states to
% see how similar they are

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

disp('kerncomp4.m : back to a more traditional strf');

% P-VALUES FOR SIGNIFICANCE IN VARIOUS TESTS
PFILT=0.010;  % p value (any attention) to be considered stimulus-modulated
PATT=0.008;   % p<0.008 any pair has p<ATTA attention modulation
PATTA=0.05;   % p for attention modulation

% load targets/patches and convert to same format as used for
% analysis

if params.runclassid==10, % ie, dms
   z=load(params.respfiles{1});
   bigpatches=z.target;
   [fmask,crop]=movfmask(size(z.target,1),0,size(z.target,1));
   targpatches=movresize(floor(z.target),params.stimloadparms{1},fmask,crop);
   ftargs=feval(params.stimfiltercmd,targpatches,...
                params.stimfilterparms{:});
   
   mpatches=[mean(ftargs,2),ftargs];
   fpatches=ftargs;
   
else   %runclass==4 -- fvvs
   loadfreetargs;
end

%
% LOAD AND CONVERT ALL STIM AND RESPONSE FILES TO APPROPRIATE FORMAT
%

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
      validtargs=unique(z.targets(find(z.phase>-1 & ~z.result)));
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

clear z

% optionally only do analysis on att states with biggest mean
% difference in response
if 0,
   respidx=min([size(bresp,2) 7]);
   nresp=squeeze(nansum(~isnan(bresp(:,respidx,:))));
   mresp=squeeze(nanmean(bresp(:,respidx,:)));
   attmin=find(mresp(2:end)==min(mresp(2:end)))+1;
   attmax=find(mresp(2:end)==max(mresp(2:end)))+1;
   bresp=bresp(:,:,[1 attmin attmax]);
   %attlong=find(nresp>100);
   %bresp=bresp(:,:,attlong);
end


% figure out size and mean of response
blen=size(bresp,1);
respcount=size(bresp,2);
attcount=size(bresp,3);

spacecount=size(bstim,2);
bsmean=mean(bstim,1);
bstim=bstim-repmat(bsmean,[blen,1]);


%
% BEGIN STRF CALCULATION
%
% 0. compute STRF in the PC domain
disp('0. Testing for linear modulation...');

DODAMPING=0;

DAMPFACTOR=0.8;
bootcount=20;
spacelims=128;
%spacelims=35;
xccount=1;
xcidx=1;
spaceuses=1:spacelims;
predlims=spacelims;

% if active cutoff selection is on, this
% will be replaced after the first iteration
predspace=spaceuses;

if params.randexp,
   kerncount=params.noisecount+1;
else
   kerncount=1;
end

% initialize kernels and other parameters for the various
% attentional models plus randomization over the hierarchically
% introduced parameter

% dc offset for each attentional state
c1=zeros(respcount,attcount,kerncount,bootcount);

% global gain for each attentional kernel
a1=zeros(respcount,attcount,kerncount,bootcount);

% baseline
bH=zeros(spacelims,respcount,kerncount,bootcount);

% DC: dc shift
dH=zeros(spacelims,respcount,kerncount,bootcount);

% GL: global gain
gH=zeros(spacelims,respcount,kerncount,bootcount);

% LO: local gain
lH=zeros(spacelims,respcount,attcount,kerncount,bootcount);

% prediction scores for real and randomized attention
% conditions, exploratory and confirmatory response data
paircount=(attcount-1)*(attcount-2)/2;
xc=zeros(params.noisecount+1,4,xccount);
xcp=zeros(params.noisecount+1,4,paircount,xccount);
xccross=zeros(params.noisecount+1,4,attcount,xccount);

% generate PC domain basis set for projecting all Hs
aokidx=find(sum(~isnan(bresp(:,end,2:end)),3));
stim=bstim(aokidx,:);
sm=mean(stim,1);
stim=stim-repmat(sm,[size(stim,1),1]);
sSAfull=stim'*stim ./ length(aokidx);
[ua,sa,va]=svd(sSAfull);
eigstim=bstim*ua;

% copy resp and subtract mean from each attentional state in
% order to support randomization over gain/local effects only
% (ie, randomized kernels have appropriate dc shifts)
brespnodc=bresp;
for attidx=2:attcount,
   for respidx=1:respcount,
      aokidx=find(~isnan(brespnodc(:,respidx,attidx)));
      brespnodc(aokidx,respidx,attidx)=...
          brespnodc(aokidx,respidx,attidx)-...
          mean(brespnodc(aokidx,respidx,attidx));
      brespnodc(aokidx,respidx,1)=brespnodc(aokidx,respidx,attidx);
   end
end
meanadjust=bresp(:,:,1)-brespnodc(:,:,1);

% figure out baseline/dc/global/local kernels to produce
% residuals for fitting in later steps
% global decorr kernel
respidx=1;
sst=[1./sa(find(sa>0)); zeros(length(find(diag(sa)==0)),1)];
B=diag(sst(spaceuses)) * ua(:,spaceuses)';

% figure out optimal global gain for att all
% to pred current attidx responses
okidx=find(sum(~isnan(bresp(:,end,2:end)),3));
slen=length(okidx);
r=bresp(okidx,respidx,1)-mean(bresp(okidx,respidx,1));
rnodc=brespnodc(okidx,respidx,1)-mean(brespnodc(okidx,respidx,1));
stim=bstim(okidx,:);
sm=mean(stim,1);
stim=stim-repmat(sm,[slen,1]);
estim=eigstim(okidx,spaceuses)-repmat(sm * ua(:,spaceuses),[slen 1]);

bSR=stim'*r./slen;
bHpc=B*bSR;

dSR=stim'*rnodc./slen;
dHpc=B*dSR;

g1=zeros(respcount,attcount);
gainadjust=zeros(size(bresp,1),1);
brespnoglob=brespnodc;
for attidx=1:attcount,
   if attidx==1,
      lokidx=okidx;
   else
      lokidx=find(sum(~isnan(bresp(:,end,attidx)),3));
   end
   slen=length(lokidx);
   stim=bstim(lokidx,:);
   sm=mean(stim,1);
   estim=eigstim(lokidx,spaceuses)-repmat(sm * ua(:,spaceuses),[slen 1]);
   
   r0=brespnodc(lokidx,respidx,1)-mean(brespnodc(lokidx,respidx,1));
   r1=estim * dHpc;
   d1=sum(r1.^2);
   if d1>0,
      g1(respidx,attidx)=sum(r0.*r1)./d1;
   end
   if ~g1(respidx,attidx);
      g1(respidx,attidx)=1;
   end
   if attidx>1,
      brespnoglob(lokidx,respidx,1)=...
          brespnodc(lokidx,respidx,1)./g1(respidx,attidx);
      brespnoglob(lokidx,respidx,attidx)=...
          brespnodc(lokidx,respidx,attidx)./g1(respidx,attidx);
   end
   if attidx>1,
      gainadjust(lokidx)=g1(attidx);
   end
end
%keyboard
   
% figure out number of valid fixations in each attentional condition
ncount=squeeze(sum(~isnan(bresp(:,end,:)),1));
cumncount=[0;cumsum(ncount(2:end))];
nrange=find(sum(~isnan(bresp(:,end,2:end)),3));

% choose a random set of fixations for each bootstrap. this should
% avoid bias in the non-randomized noiseidx=1
vidx=nrange;
vcount=length(vidx);
[vv,ss]=sort(rand(size(vidx)));
vidx=vidx(ss);

% generate random attention state sets
sn=zeros(length(nrange),params.noisecount+1);
for noiseidx=2:params.noisecount+1;
   [tt,sn(:,noiseidx)]=sort(rand(size(nrange)));
end

fprintf('DODAMPING=%d  randexp=%d  randcnf=%d length(spaceuses)=%d\n',...
        DODAMPING,params.randexp,params.randcnf,length(spaceuses));

for noiseidx=1:kerncount,
   fprintf('\nnoiseidx=%d: ',noiseidx);
   
   resp=ones(size(bresp))*nan;
   resp(nrange,:,1)=bresp(nrange,:,1);
   if noiseidx==1,
      resp(:,:,2:end)=bresp(:,:,2:end);
   else
      % pick random atts
      for attidx=2:attcount,
         nn=nrange(sn(cumncount(attidx-1)+1:cumncount(attidx),noiseidx));
         resp(nn,:,attidx)=bresp(nn,:,1);
      end
   end
   
   %keyboard
   smean=zeros(spacecount,respcount,attcount,bootcount);
   
   for bootidx=1:bootcount;
      fprintf('.');
      
      % generate exploratory response matrix that is resp
      % with cnf data removed
      eresp=resp;
      ppidx=vidx(round((bootidx-1)/bootcount*vcount+1):...
                 round(bootidx/bootcount*vcount));
      eresp(ppidx,:,:)=nan;
      
      % find DC offset for this bootstrap
      ermean=reshape(nanmean(eresp),respcount,attcount);
      ermean(find(isnan(ermean)))=0;
      c1(:,:,noiseidx,bootidx,xcidx)=ermean;
      
      for attidx=1:attcount,
         
         % find spatial correlation for stim in attidx
         
         % non-zero frames for this att/resample
         aokidx=find(~isnan(eresp(:,end,attidx)));
         
         % set up B to transform strf from attention state idx
         % into PC domain of all att, decorrelating as needed.
         if length(aokidx)>0,
            
            % ssa=u*s*v'
            % ssa_pc=ua'* (u * s * v') * ua
            % ssa^-1=v*s^-1*u'
            % ssa_pc^-1 = ua'* v* s^-1 * u'* ua
            if 0 & attidx==1,
               % basically just have global correlations.
               sst=[1./sa(find(sa>0)); ...
                    zeros(length(find(diag(sa)==0)),1)];
               B=diag(sst(spaceuses)) * ua(:,spaceuses)';
            else
               stim=bstim(aokidx,:);
               sm=mean(stim,1);
               stim=stim-repmat(sm,[size(stim,1),1]);
               tSA=stim'*stim ./ length(aokidx);
               
               if 0,
                  [ut,st,vt]=svd(tSA);
                  sst=[1./st(find(st>0)); ...
                       zeros(length(find(diag(st)==0)),1)];
                  B=ua(:,spaceuses)' * vt * diag(sst) * ut';
               else
                  % alt, assume that all kernels have same
                  % decorr matrix
                  st=diag(ua' * tSA *ua);
                  sst=[1./st(find(st>0)); ...
                       zeros(length(find(st<=0)),1)];
                  B=diag(sst(spaceuses)) * ua(:,spaceuses)';
               end
            end
            
            for respidx=1:respcount,
               okidx=find(~isnan(eresp(:,respidx,attidx)));
               slen=length(okidx);
               if slen>0,
                  
                  stim=bstim(okidx,:);
                  smean(:,respidx,attidx,bootidx)=mean(stim,1)';
                  stim=stim-repmat(smean(:,respidx,attidx,bootidx)',...
                                   [slen,1]);
                  
                  if attidx==1,
                     
                     % fit baseline kernel for this att
                     % state... should be the same for every one
                     r=bresp(okidx,respidx,1);
                     r=r-mean(r);
                     
                     SR=stim'*r./slen;
                     bH(:,respidx,noiseidx,bootidx)=B*SR;
                     
                     % fit dc mod kernel for this att state
                     r=bresp(okidx,respidx,1);
                     
                     % subtract dc of current channels
                     for aa=2:attcount,
                        ii=find(~isnan(eresp(okidx,respidx,aa)));
                        r(ii)=r(ii)-mean(r(ii));
                     end
                     
                     % do the RC
                     r=r-mean(r);
                     
                     SR=stim'*r./slen;
                     dH(:,respidx,noiseidx,bootidx)=B*SR;
                     
                     % find real dc-subtracted kernel for
                     % fitting global gain
                     r=brespnodc(okidx,respidx,1);
                     r=r-mean(r);
                     SR=stim'*r./slen;
                     gH(:,respidx,noiseidx,bootidx)=B*SR;
                  end
                  
                  % figure out optimal global gain
                  % to pred current attidx responses
                  
                  estim=eigstim(okidx,spaceuses)-...
                        repmat(smean(:,respidx,attidx,bootidx)' ...
                               * ua(:,spaceuses), [slen 1]);
                  r=brespnodc(okidx,respidx,1);
                  r=r-mean(r);
                  r1=estim * gH(spaceuses,respidx,noiseidx,bootidx);
                  d1=sum(r1.^2);
                  if d1>0,
                     a1(respidx,attidx,noiseidx,bootidx,xcidx)=...
                         sum(r.*r1)./d1;
                  else
                     fprintf('z');
                     a1(respidx,attidx,noiseidx,bootidx,xcidx)=1;
                  end
                  
                  % fit local kernel
                  r=brespnoglob(okidx,respidx,1);
                  r=r-mean(r);
                  SR=stim'*r./slen;
                  lH(:,respidx,attidx,noiseidx,bootidx)=B*SR;
                  if sum(isnan(lH(:,respidx,attidx,noiseidx,bootidx))),
                     disp('nan lH');
                     keyboard
                  end
               end
            end
         end
         if exist('BATQUEUEID','var') & BATQUEUEID>0,
            dbsetqueue(BATQUEUEID,1);
         end
      end
   end
end

noiseidx=1;

% do prediction test noatt/att comparison
respidx=1;
ract=bresp(:,respidx,1);
ractnodc=brespnodc(:,respidx,1);
ractnoglob=brespnoglob(:,respidx,1);

% figure out optimal damping... currently not used for anything
   r0=zeros(blen,1);
   
   Hm=mean(bH(:,respidx,1,:),4);
   Hs=std(bH(:,respidx,1,:),1,4) .* sqrt(bootcount);
   Hs(find(Hs==0))=1;
   
   ddr=0.3:0.1:1.5;
   xct=zeros(length(ddr),1);
   agoodidx=find(sum(~isnan(bresp(:,end,2:attcount)),3));
   
   for dd=1:length(ddr),
      Hdamp=(abs(Hm)./(Hs.*ddr(dd)));
      Hdamp=(1-Hdamp.^(-2));
      Hdamp=Hdamp.*(Hdamp>0);
      Hdamp(find(isnan(Hdamp)))=0;
      
      for bootidx=1:bootcount,
         a0idx=vidx(round((bootidx-1)/bootcount*vcount+1):...
                    round(bootidx/bootcount*vcount));
         
         % pull out stimulus
         estim=eigstim(a0idx,spaceuses)- ...
               repmat(smean(:,respidx,1,bootidx)' * ...
                      ua(:,spaceuses), [length(a0idx) 1]);
         
         % prepare kernels for prediction with damping
         H0=bH(predspace,respidx,1,bootidx).*Hdamp(predspace,1);
         
         % no attention prediction
         r0(a0idx)=estim * H0;
         r0(a0idx)=r0(a0idx) + c1(respidx,1,noiseidx,bootidx,xcidx);
      end
      xct(dd)=xcov(r0(agoodidx),ract(agoodidx),0,'coeff');
   end
   
   % generate optimal Hdamp for the remainder of predictions
   ddmax=min(find(xct==max(xct)));
   Hdamp=(abs(Hm)./(Hs.*ddr(ddmax)));
   Hdamp=(1-Hdamp.^(-2));
   Hdamp=Hdamp.*(Hdamp>0);
   Hdamp(find(isnan(Hdamp)))=0;
   if sum(Hdamp)==0,
      disp('lame situation: Hdamp = all zeros.');
   end
   Hdamp=repmat(Hdamp,[1 attcount]);
   
   fprintf('\nmaximum DAMPFACTOR=%.1f\n',ddr(ddmax));

% determine optimal number of eigenvectors to include in
% kernel using all-att condition
%keyboard
if ~DODAMPING,
   nidx=1;  % use first kernels, ie the real ones
   r0=zeros(blen,spacelims);
   rnoatt=zeros(blen,spacelims);
   rdcatt=zeros(blen,spacelims);
   rglobalatt=zeros(blen,spacelims);
   rfullatt=zeros(blen,spacelims);
   xctb=zeros(spacelims,bootcount,4);
   for bootidx=1:bootcount,
      a0idx=vidx(round((bootidx-1)/bootcount*vcount+1):...
                 round(bootidx/bootcount*vcount));
      
      % pull out stimulus
      estim=eigstim(a0idx,spaceuses)- ...
            repmat(smean(:,respidx,1,bootidx)' * ...
                   ua(:,spaceuses), [length(a0idx) 1]);
      
      % prepare kernels for prediction no damping
      H0=bH(spaceuses,respidx,nidx,bootidx);
      H1=dH(spaceuses,respidx,nidx,bootidx);
      H2=gH(spaceuses,respidx,nidx,bootidx);
      H=squeeze(lH(spaceuses,respidx,2:end,nidx,bootidx));
      
      % no attention prediction
      r0(a0idx,:)=estim .* repmat(H0',[length(a0idx),1]);
      %r0(a0idx)=estim * H0;
      rnoatt(a0idx,:)=r0(a0idx,:) + c1(respidx,1,nidx,bootidx,xcidx);
      
      rdcatt(a0idx,:)=estim .* repmat(H1',[length(a0idx),1]);
      rglobalatt(a0idx,:)=estim .* repmat(H1',[length(a0idx),1]);
      
      for attidx=2:attcount,
         % find prediction fixations for this bootstrap and attidx
         aidx=a0idx(find(~isnan(resp(a0idx,respidx,attidx))));
         
         % dc offset prediction
         rdcatt(aidx)=rdcatt(aidx) + ...
             c1(respidx,attidx,nidx,bootidx,xcidx);
         
         % global gain prediction
         rglobalatt(aidx)=rglobalatt(aidx) .* ...
             a1(respidx,attidx,nidx,bootidx,xcidx);
         
         % local gain prediction
         estim=eigstim(aidx,predspace)-...
               repmat(smean(:,respidx,attidx,bootidx)' ...
                      * ua(:,predspace), [length(aidx) 1]);
         rfullatt(aidx,:)=estim .* repmat(H(:,attidx-1)',[length(aidx),1]);;
      end
      
      for ss=1:spacelims,
         if std(ract(a0idx))>0 & std(abs(sum(r0(a0idx,1:ss))))>0,
            xctb(ss,bootidx,1)=xcov(sum(rnoatt(a0idx,1:ss),2),ract(a0idx),...
                                  0,'coeff');
            xctb(ss,bootidx,2)=xcov(sum(rdcatt(a0idx,1:ss),2),ract(a0idx),...
                                    0,'coeff');
            xctb(ss,bootidx,3)=xcov(sum(rglobalatt(a0idx,1:ss),2),...
                                    ractnodc(a0idx),0,'coeff');
            xctb(ss,bootidx,4)=xcov(sum(rfullatt(a0idx,1:ss),2),...
                                    ractnoglob(a0idx),0,'coeff');
         end
      end
   end
      
   agoodidx=find(sum(~isnan(bresp(:,end,2:attcount)),3));
   xct=zeros(spacelims,4);
   for ss=1:spacelims,
      xct(ss,1)=xcov(sum(rnoatt(agoodidx,1:ss),2),ract(agoodidx),...
                   0,'coeff');
      xct(ss,2)=xcov(sum(rdcatt(agoodidx,1:ss),2),ract(agoodidx),...
                   0,'coeff');
      xct(ss,3)=xcov(sum(rglobalatt(agoodidx,1:ss),2),ractnodc(agoodidx),...
                   0,'coeff');
      xct(ss,4)=xcov(sum(rfullatt(agoodidx,1:ss),2),ractnoglob(agoodidx),...
                   0,'coeff');
   end
   
   % smooth with boxcar to reduce local maxima
   xct=rconv2(xct,[1/3 1/3 1/3]');
   xctb=rconv2(xctb,[1/3 1/3 1/3]');
   
   %emax=min(xct==max(xct));
   emax=min(find((xct(:,1)-min(xct(:,1)))>=...
                 (max(xct(:,1))-min(xct(:,1))).*0.95));
   bemax=zeros(1,bootcount);
   for bootidx=1:bootcount,
      if sum(~isnan(xctb(:,bootidx)))>0,
         bemax(bootidx)=min(find((xctb(:,bootidx,1)-min(xctb(:,bootidx,1)))>= ...
                         (max(xctb(:,bootidx,1))-min(xctb(:,bootidx,1))).*0.95));
      end
   end
   
   predlims=emax;
   %predlims=round(mean(bemax));
   %predlims=ceil(0.6 .* emax);
   if predlims==1,
      predlims=2;
   end
   
   predspace=1:predlims;
   
   % figure out scaling factor
   ra=sum(r0(agoodidx,1:predlims),2);
   ra=ra-mean(ra);
   rb=ract(agoodidx);
   rb=rb-mean(rb);
   
   d1=sum(rb.^2);
   if d1>0,
      scf=sum(ra.*rb)./d1;
   else
      scf=1;
   end
   
   % apply scaling factor to all kernels so that they'll be about
   % the right size for computing mse. this gain adjustment has
   % nothing to do with attention, and thus won't mess up xc
   % attention tests.
   %bH=bH.*scf;
   %dH=dH.*scf;
   %gH=gH.*scf;
   %lH=lH.*scf;
   
   fprintf('90%% maximum pred at predlims=%d eigs (emax=%d bemax=%d)\n',...
           predlims,emax,round(mean(bemax)));
end

for noiseidx=1:params.noisecount+1,
   fprintf('noiseidx=%d: ',noiseidx);
   
   resp=ones(size(bresp))*nan;
   resp(nrange,:,1)=bresp(nrange,:,1);
   if noiseidx==1 | ~params.randcnf,
      resp(:,:,2:end)=bresp(:,:,2:end);
   else
      % pick random atts
      for attidx=2:attcount,
         nn=nrange(sn(cumncount(attidx-1)+1:cumncount(attidx),noiseidx));
         resp(nn,:,attidx)=bresp(nn,:,1);
      end
   end
   if params.randexp,
      nidx=noiseidx;
   else
      nidx=1;  % always same exp data
   end
   
   % randomize cnf data with each noise kernel
   % no, always pred the same data
   %resp=ones(size(bresp))*nan;
   %resp(nrange,:,:)=bresp(nrange,:,:);
   % nidx=noiseidx;
   
   % set up vectors to hold preds from different models
   r0=zeros(blen,1);
   rnoatt=zeros(blen,1);
   rdcatt=zeros(blen,1);
   rglobalatt=zeros(blen,1);
   rfullatt=zeros(blen,1);
   
   % compute contribution of each bootstrap component to the preds.
   for bootidx=1:bootcount,
      if DODAMPING==2,
         Hm=mean(lH(:,respidx,:,1,...
                    [1:bootidx-1 bootidx+1:bootcount]),5);
         Hs=std(lH(:,respidx,:,1,...
                   [1:bootidx-1 bootidx+1:bootcount]),1,5)...
            .*sqrt(bootcount-1);
         
         % slightly cheating method for generating Hdamp.
         %Hm=mean(lHpc(:,:,:,1,:),5);
         %Hs=std(lHpc(:,:,:,1,:),1,5).*sqrt(bootcount);
         %Hs=repmat(std(lHpc(:,:,1,1,:),1,5).*sqrt(bootcount),
         %     [1 1 attcount]);
         
         %
         % shrinkage filter. is this better? seems to help a
         % lot! but it's kind of cheating.
         % 
         Hs(find(Hs==0))=1;
         Hdamp=(abs(Hm)./(Hs.*DAMPFACTOR));
         Hdamp=(1-Hdamp.^(-2));
         Hdamp=Hdamp.*(Hdamp>0);
         Hdamp(find(isnan(Hdamp)))=0;
         
         % problem--damping reduces power of kernel, abs
         % variance of prediction falls wrt actual response.
         %
         % BAD SOLUTION:
         % rescale to keep power of pred response about right.
         %for attidx=1:attcount,
         %   if sum(Hdamp(:,attidx))>0,
         %      Hdamp(:,attidx)=Hdamp(:,attidx)./ ...
         %          (sum(Hm(:,attidx).^2.*Hdamp(:,attidx).*diag(sa)) ...
         %           ./ sum(Hm(:,attidx).^2.*diag(sa)));
         %   end
         %end
      end
      
      % moved bootidx loop earlier so that Hdamp is not
      % biased by confirmatory set.  is this all i need to do
      % to make it legitimately unbiased?
      
      % find cnf fixations for this bootstrap
      a0idx=vidx(round((bootidx-1)/bootcount*vcount+1):...
                 round(bootidx/bootcount*vcount));
      
      % pull out stimulus
      %stim=bstim(a0idx,:);
      %stim=stim-repmat(smean(:,respidx,1,bootidx)',...
      %                 [size(stim,1),1]);
      %estim=stim*ua(:,predspace);
      estim=eigstim(a0idx,predspace)- ...
            repmat(smean(:,respidx,1,bootidx)' * ...
                   ua(:,predspace), [length(a0idx) 1]);
      
      % prepare kernels for prediction
      if ~DODAMPING,
         % no damping
         H0=bH(predspace,respidx,nidx,bootidx);
         H1=dH(predspace,respidx,nidx,bootidx);
         H2=gH(predspace,respidx,nidx,bootidx);
         H=squeeze(lH(predspace,respidx,2:end,nidx,bootidx));
      else
         % damping on
         H0=bH(predspace,respidx,nidx,bootidx).*Hdamp(predspace,1);
         H1=dH(predspace,respidx,nidx,bootidx).*Hdamp(predspace,1);
         H2=gH(predspace,respidx,nidx,bootidx).*Hdamp(predspace,1);
         H=squeeze(lH(predspace,respidx,2:end,nidx,bootidx)) .* ...
           Hdamp(predspace,2:end);
      end
      
      % no attention prediction
      r0(a0idx)=estim * H0;
      rnoatt(a0idx)=r0(a0idx) + c1(respidx,1,nidx,bootidx,xcidx);
      
      rdcatt(a0idx)=estim * H1;
      rglobalatt(a0idx)=estim * H2;
      
      for attidx=2:attcount,
         % find prediction fixations for this bootstrap and attidx
         aidx=a0idx(find(~isnan(resp(a0idx,respidx,attidx))));
         
         % dc offset prediction
         rdcatt(aidx)=rdcatt(aidx) + ...
             c1(respidx,attidx,nidx,bootidx,xcidx);
         
         % global gain prediction
         rglobalatt(aidx)=rglobalatt(aidx) .* ...
             a1(respidx,attidx,nidx,bootidx,xcidx);
         
         % local gain prediction
         estim=eigstim(aidx,predspace)-...
               repmat(smean(:,respidx,attidx,bootidx)' ...
                      * ua(:,predspace), [length(aidx) 1]);
         rfullatt(aidx)=estim * H(:,attidx-1);
      end
   end
   
   % compensate for dc effect (built in)
   rglobalatt=rglobalatt+meanadjust;
   rfullatt=rfullatt.*gainadjust+meanadjust;
   
   % any attention state
   agoodidx=find(sum(~isnan(bresp(:,end,2:attcount)),3));
   
   % count how many PCs are significantly non-zero for each
   % attentional state (for actual local attention strfs)
   if noiseidx==1 & xcidx==1,
      nondampedcount=squeeze(sum(Hdamp(:,respidx,:)>0,1));
   end
   
   % for local attention only choose local atttention states
   % that had any significant power in their STRF
   if DODAMPING==2,
      goodatt=find(nondampedcount(2:end)>0)'+1;
   else
      goodatt=2:attcount; % old way, no excluded attidx's
   end
   alocalidx=find(sum(~isnan(resp(:,end,goodatt)),3));
   
   if 0 & xcidx==2,
      % add back dc shifts to make xc scores comparable
      rnoatt=rnoatt+meanadjust;
      rdcatt=rdcatt+meanadjust;
      rglobalatt=rglobalatt+meanadjust;
      rfullatt=rfullatt+meanadjust;
   end
   
   % rectify?  doesn't seem to help worth dog poo poo
   if 0,
      rnoatt(find(rnoatt<0))=0;
      rdcatt(find(rdcatt<0))=0;
      rglobalatt(find(rglobalatt<0))=0;
      rfullatt(find(rfullatt<0))=0;
   end
   
   if length(agoodidx)>0,
      xc(noiseidx,1,xcidx)=xcov(rnoatt(agoodidx),ract(agoodidx),...
                                0,'coeff');
      xc(noiseidx,2,xcidx)=xcov(rdcatt(agoodidx),ract(agoodidx),...
                                0,'coeff');
      xc(noiseidx,3,xcidx)=xcov(rglobalatt(agoodidx),...
                                ract(agoodidx),0,'coeff');
      if length(alocalidx)>0,
         xc(noiseidx,4,xcidx)=xcov(rfullatt(alocalidx),...
                                   ract(alocalidx),0,'coeff');
      end
      fprintf('xc%d: %.3f %.3f %.3f %.3f\n',xcidx,xc(noiseidx,:,xcidx));
   else
      disp('zero len, xc=0');
   end
   
   
   if noiseidx==1,
      fprintf('ncount: %d %d %d %d\n',ncount(2:end));
      fprintf('(valid att conds = %d/%d [%d %d %d %d %d',...
              length(goodatt),attcount-1,nondampedcount);
      fprintf('/%d])\n',spacelims);
   end
   % pairwise prediction test... can you find an attention
   % pair that shows more/less modulation?
   pidx=0;
   pairlist=zeros(paircount,2);
   for p1=2:attcount-1,
      for p2=(p1+1):attcount,
         pidx=pidx+1;
         pairlist(pidx,:)=[p1 p2];
         aattidx=find(sum(~isnan(bresp(:,end,[p1 p2])),3));
         xcp(noiseidx,1,pidx,xcidx)=xcov(rnoatt(aattidx),...
                                         ract(aattidx),0,'coeff');
         xcp(noiseidx,2,pidx,xcidx)=xcov(rdcatt(aattidx),...
                                         ract(aattidx),0,'coeff');
         xcp(noiseidx,3,pidx,xcidx)=xcov(rglobalatt(aattidx),...
                                         ract(aattidx),0,'coeff');
         xcp(noiseidx,4,pidx,xcidx)=xcov(rfullatt(aattidx),...
                                         ract(aattidx),0,'coeff');
      end
   end
   
   for p1=1:attcount,
      if p1==1,
         aattidx=agoodidx;
      else
         aattidx=find(~isnan(bresp(:,end,p1)));
      end

      %xccross(noiseidx,1,p1,xcidx)=xcov(rnoatt(aattidx),...
      %                                  ract(aattidx),0,'coeff');
      %xccross(noiseidx,2,p1,xcidx)=xcov(rdcatt(aattidx),...
      %                                  ract(aattidx),0,'coeff');
      %xccross(noiseidx,3,p1,xcidx)=xcov(rglobalatt(aattidx),...
      %                                  ract(aattidx),0,'coeff');
      %xccross(noiseidx,4,p1,xcidx)=xcov(rfullatt(aattidx),...
      %                                  ract(aattidx),0,'coeff');
      xccross(noiseidx,1,p1,xcidx)=...
          sqrt(mean((rnoatt(aattidx)-ract(aattidx)).^2))./std(ract(aattidx));
      xccross(noiseidx,2,p1,xcidx)=...
          sqrt(mean((rdcatt(aattidx)-ract(aattidx)).^2))./std(ract(aattidx));
      xccross(noiseidx,3,p1,xcidx)=...
          sqrt(mean((rglobalatt(aattidx)-ract(aattidx)).^2))./std(ract(aattidx));
      xccross(noiseidx,4,p1,xcidx)=...
          sqrt(mean((rfullatt(aattidx)-ract(aattidx)).^2))./std(ract(aattidx));
   end
   if sum(isnan(xc(noiseidx,:,xcidx))),
      keyboard
   end
   if mod(noiseidx,50)==1,
      dbsetqueue;
   end
end

predxc=cat(1,xc(1,:,:),mean(xc(2:end,:,:),1),std(xc(2:end,:,:),0,1));
predxccross=cat(1,xccross(1,:,:),mean(xccross(2:end,:,:),1),...
                std(xccross(2:end,:,:),0,1));

%keyboard

%
% ALL ATT SIGNIFICANCE TESTS
%

% test for pred differences
pxc=ones(1,size(xc,2),2);
pxcp=zeros(1,size(xcp,2),size(xcp,3),2);
pxccross=zeros(1,size(xccross,2),size(xccross,3),2);
if attcount>2,   % only if there's more than one att state
   for xx=1:size(xc,2),
      nn=length(find(xc(2:end,xx)>=xc(1,xx)));
      pxc(1,xx,1)=(nn+1)./(params.noisecount+1);
   end
   
   % test for DC tuning differences
   dd=reshape(mean(c1(:,2:end,:,:),4),attcount-1,kerncount);
   ds=reshape(std(c1(:,2:end,:,:),0,4),attcount-1,kerncount);
   ds(find(ds==0))=1;
   dd=dd-repmat(mean(dd,1),[attcount-1 1]);
   dd=dd./ds;
   dd=sqrt(mean(dd.^2,1));
   nn=length(find(dd>=dd(1)));
   pxc(1,2,2)=nn./(params.noisecount+1);
   
   % test for GLOBAL gain tuning differences
   dd=reshape(mean(a1(:,2:end,:,:),4),attcount-1,kerncount);
   ds=reshape(std(a1(:,2:end,:,:),0,4),attcount-1,kerncount);
   ds(find(ds==0))=1;
   dd=dd-repmat(mean(dd,1),[attcount-1 1]);
   dd=dd./ds;
   dd=sqrt(mean(dd.^2,1));
   nn=length(find(dd>=dd(1)));
   pxc(1,3,2)=nn./(params.noisecount+1);
   
   % test for LOCAL gain tuning differences
   dd=reshape(mean(lH(predspace,:,2:end,:,:),5),...
              predlims,attcount-1,kerncount);
   ds=reshape(std(lH(predspace,:,2:end,:,:),0,5),...
              predlims,attcount-1,kerncount);
   ds(find(ds==0))=1;
   dd=dd-repmat(mean(dd,2),[1 attcount-1 1]);
   dd=dd./ds;
   dd=sqrt(mean(mean(dd.^2,1),2));
   nn=length(find(dd>=dd(1)));
   pxc(1,4,2)=nn./(params.noisecount+1);
   
   %
   % PAIRED ATT SIGNIFICANCE TESTS
   %
   % test for pred diffs
   if paircount>0,
      for xx=1:size(xcp,2)*size(xcp,3),
         nn=length(find(xcp(2:end,xx)>=xcp(1,xx)));
         pxcp(xx)=(nn+1)./(params.noisecount+1);
      end
   end
   
   % "XC" is actually MSE in this case, just for normalization purposes.
   for xx=1:size(xccross,2)*size(xccross,3),
      nn=length(find(xccross(2:end,xx)<=xccross(1,xx)));
      pxccross(xx)=(nn+1)./(params.noisecount+1);
   end
end

fprintf('XC:    %.3f %.3f %.3f %.3f\n',predxc(1,:));
fprintf('Noise: %.3f %.3f %.3f %.3f\n',predxc(2,:));
fprintf('PXC:   %.3f %.3f %.3f %.3f\n',pxc);
disp('');
fprintf('CRXC:  %.3f %.3f %.3f %.3f\n',predxccross(1,:,:));
fprintf('Noise: %.3f %.3f %.3f %.3f\n',predxccross(2,:,:));

%targsim=zeros(attdocount,spacelim,2);
%attcc=zeros(spacelim,attuse,3);
%p=0;
pstrf=ones(spacelims,respcount,attcount);


clear x1 x2 r2 r2std resp stim eresp presp
clear data rbinset y1 y2 ystd *idx tSA timage sSAfull
clear z sn timage tresp tstim tt stim estim eigstim
clear fmask alphamask rfullatt rglobalatt
