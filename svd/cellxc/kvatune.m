% function res=kvatune(cellid,batch,showextra)
%
% load results from resfile in sRunData and display tuning
% info ... resfile generated from kernvsatt or xcdms
%
% showextra - default = 'none'
%
% created SVD 10/31/04 - hacked from dmsres.m
%
function res=kvatune(cellid,batch,showextra)

dbopen;

if exist('batch','var') & batch>0,
   fprintf('%s: (cell %s bat %d)\n',mfilename,cellid,batch);
   goodbatch=zeros(1,length(batch));
   batchcount=0;
   for ii=1:length(batch),
      sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
           ' AND batch=',num2str(batch(ii))];
      trd=mysql(sql);
      if ~isempty(trd),
         goodbatch(ii)=1;
         batchcount=batchcount+1;
         rundata(batchcount)=trd;
      end
   end
   batch=batch(find(goodbatch));
   resfile=[rundata(1).respath,rundata(1).kernfile,'.gz'];
else
   % no batch specified, assume it's a file.
   resfile=cellid;
   if ~exist(resfile,'file'),
      batchcount=0;
   else
      batchcount=1;
   end
end

if batchcount==0,
   disp('no entries found in db!');
   if nargout>0,
      r=0;
   end
   return
end

if ~exist('showextra','var'),
   showextra='none';
end
showextra=strsep(showextra,':');
if length(showextra)>1,
   prefframe=showextra{2};
   showextra=showextra{1};
end

global GCOLORMAP
GCOLORMAP=redblue;

rcsetstrings;

if ~exist(resfile,'file'),
   USERESFILE=1;
   resfile=[rundata(1).respath,rundata(1).resfile,'.gz'];
else
   USERESFILE=0;
end

if ~exist(resfile,'file'),
   fprintf('%s not found!\n',resfile);
   res.cellid=cellid;
   res.cc=nan;
   res.bad=1;
   return
end

fprintf('Loading: %s\n',resfile);
zload(resfile);

cellid=params.cellid;
batch=params.batch;

if batch==77,
   USERESFILE=1;
end

if ~isfield(params,'bootcount'),
   if isfield(params,'fitboot') & params.fitboot>0,
      params.bootcount=params.fitboot;
   else
      params.bootcount=1;
   end
elseif params.bootcount==0,
   params.bootcount=1;
end

% raw kernels
if USERESFILE,
   
   % for V1 data (cellxcmaster output)
   disp('processing V1 cell');
   H=zeros([size(strf(params.nlidxsave,1).h,1),params.bootcount]);
   mS=zeros([size(strf(params.nlidxsave,1).mS,1),params.bootcount]);
   nlparms=zeros(1,params.bootcount,attcount);
   for bootidx=1:params.bootcount,
      H(:,bootidx)=sum(strf(params.nlidxsave,bootidx).h,2) ./ ...
          length(strf(params.nlidxsave,bootidx).tempresp);
      mS(:,bootidx)=strf(params.nlidxsave,bootidx).mS;
      nlparms(bootidx)=-strf(params.nlidxsave,bootidx).nlparms;
   end
   mH=mean(H,2);
   for bootidx=1:params.bootcount,
      if H(:,bootidx)'*mH < 0,
         H(:,bootidx)=-H(:,bootidx);
      end
   end
   
   blen=0;
   randxc=predres(end).prederr(params.nlidxsave);
   %meanactresp=nanmean(predres(end).act_resp{1});
   
   meanactresp=nanmean(resp(:)).*72;
   choserange=params.nlidxsave:strfcount:size(predres(end).mod_psth{1},3);
   meanpredresp=nanmean(predres(end).mod_psth{1}(:,choserange));
   scale2Hz=(meanactresp./meanpredresp);
else
   
   % for V4 data (dmsatt output)
   disp('processing V4 cell');
   H=zeros([size(vstrf(1).h,1),params.bootcount,attcount]);
   mS=zeros([size(vstrf(1).mS,1),params.bootcount,attcount]);
   nlparms=zeros(2,params.bootcount,attcount);
   for attidx=1:size(vstrf,1),
      H(:,:,attidx)=cat(2,vstrf(attidx,1,1,:).h);
      mS(:,:,attidx)=cat(2,vstrf(attidx,1,1,:).mS);
      nlparms(:,:,attidx)=cat(2,vstrf(attidx,1,1,:).nlparms);
   end
   
   meanactresp=mean(nlparms(1,:,1)).*60;
   scale2Hz=repmat(60,1,params.bootcount);
end

% if extra dc dimension was specified in pfft, set those coeffs to zero
kernfmt=params.kernfmt;
if strcmp(kernfmt,'fft') | strcmp(kernfmt,'lin2'),
   phasecount=4;
elseif length(kernfmt)>4 & ...
      (kernfmt(end-1)=='+' | kernfmt(end-1)=='-'),
   phasecount=str2num(kernfmt(end));
else
   phasecount=1;
end

chancount=size(H,1)./phasecount;

showSFOR=1;
if showSFOR & ~strcmp(params.kernfmt,'space'),
   Xmax=sqrt(chancount.*2);
   kernfmt='pfftgr';
   [cfilt,cfiltconj]=gencfilt(Xmax,Xmax);
else
   Xmax=params.stimloadparms{3};
   kernfmt=params.kernfmt;
end

SHOWALLKERNS=0;
if SHOWALLKERNS,
   % show jackknifed kernels
   figure(1);
   showkern(H,kernfmt);
   
   kshowcount=min([params.bootcount 12]);
   for ii=1:kshowcount,
      for attidx=1:attcount,
         subplot(attcount,kshowcount,(attidx-1).*kshowcount+ii),
         
         if ii==1,
            title(sprintf('%s %d/%d',cellid,...
                          vstrf(attidx,1,1,ii).parms.sfsfit,...
                          vstrf(attidx,1,1,ii).parms.sigfit));
         else
            title(sprintf('%d/%d',vstrf(attidx,1,1,ii).parms.sfsfit,...
                          vstrf(attidx,1,1,ii).parms.sigfit));
         end
      end
   end
   drawnow
end

%
% FIND BEST NCART MATCH TO KERNEL
%
% load ncart stimuli
ncartfile=sprintf('/auto/k5/david/tmp/kvaparms/ncart.batch%d.mat',batch);
ncartfile2=sprintf('/auto/k5/david/tmp/kvaparms/ncart.big.batch%d.mat',batch);

if exist(ncartfile,'file'),
   fprintf('Loading ncarts: %s\n',ncartfile);
   load(ncartfile);
else
   % 'batch' has already been set to appropriate value
   genncarts;
end

BIGNCARTS=0;
if BIGNCARTS,
   
   fprintf('Using big ncart test set: %s\n',ncartfile);
   load(ncartfile2);
   
   wcount=3;
   cart=zeros(Xmax,Xmax,dcount,dcount2,dcount3,wcount,wcount);
   pol=zeros(Xmax,Xmax,dcount,dcount2,dcount3,wcount,wcount);
   hyp=zeros(Xmax,Xmax,dcount,dcount2,dcount3,wcount,wcount);
   
   mm=movfmask(size(bcart,1),0.51,40,1);
   bcart=bcart.*repmat(mm,[1 1 dcount,dcount2,dcount3]);
   bpol=bpol.*repmat(mm,[1 1 dcount,dcount2,dcount3]);
   bhyp=bhyp.*repmat(mm,[1 1 dcount,dcount2,dcount3]);
   for wx=1:wcount,
      for wy=1:wcount,
         xrange=(1:Xmax)+floor((size(bcart,1)-Xmax)./(wcount-1)*(wx-1));
         yrange=(1:Xmax)+floor((size(bcart,1)-Xmax)./(wcount-1)*(wy-1));
         cart(:,:,:,:,:,wx,wy)=bcart(xrange,yrange,:,:,:);
         pol(:,:,:,:,:,wx,wy)=bpol(xrange,yrange,:,:,:);
         hyp(:,:,:,:,:,wx,wy)=bhyp(xrange,yrange,:,:,:);
      end
   end
   
   dcount3=dcount3*wcount*wcount;
   cart=reshape(cart,Xmax*Xmax,dcount,dcount2,dcount3);
   pol=reshape(pol,Xmax*Xmax,dcount,dcount2,dcount3);
   hyp=reshape(hyp,Xmax*Xmax,dcount,dcount2,dcount3);
   
   cartprange=repmat(cartprange,[1 wcount*wcount]);
   polprange=repmat(polprange,[1 wcount*wcount]);
   hypprange=repmat(hypprange,[1 wcount*wcount]);
   
   clear bcart bpol bhyp   
else
   disp('Skipping position invariance test');
end

if ~exist('nrmov','var'),
   
   % load natrev file if it wasn't included in the cart set file
   natrevfile=['/auto/data/archive/stimarchive/REESE17/imsm/',...
               'natrev-apr-96.index60.1.pixel.imsm'];
   
   % match filter to approximately what was used for analysis
   nstimloadparms=params.stimloadparms;
   %if nstimloadparms{1}==0,
      nstimloadparms{1}=30;
   %end
   
   nrmov=loadimfile(natrevfile,1,18000,nstimloadparms{:});
end
nrmovlen=size(nrmov,3);

if ~exist('nok','var'),
   % remove frames with below threshold contrast
   
   % this operation now performed in genncarts to make everything
   % run faster (ie, ncart file is smaller, and only need to run it
   % once). of course, it means you have to re-run genncarts if you
   % want to change this threshold now.
   
   nmm=mean(reshape(nrmov,Xmax*Xmax,nrmovlen));
   nss=std(reshape(nrmov,Xmax*Xmax,nrmovlen));
   nss(find(nss<1e-13))=1;
   nok=find(nss>15);
   
   % only save relatively high contrast images
   nrmov=nrmov(:,:,nok);
   nmmrange=nmm(nok);
   nssrange=nss(nok);
   
   fprintf('Reduced nat test set from %d',nrmovlen);
   nrmovlen=size(nrmov,3);
   fprintf(' to %d frames\n',nrmovlen);
end

% normalize all the test stimuli to match mean and variance of
% stimuli actually used to measure the STRF
NORMEACHFRAME=1;
DOSPECTNORM=0;
DOEDSTIM=1;

if DOEDSTIM,
   % don't crop out ed stim
   nstimloadparms{1}=nstimloadparms{3};
   edmov=loadimfile('/auto/k5/david/sample/stim/edstim.imsm',...
                    0,0,nstimloadparms{:});
   edmovlen=size(edmov,3);
   %disp('only using contour stim in edmov');
   %edmov=edmov(:,:,1:edmovlen/3);
   %edmovlen=size(edmov,3);
end

% figure out mean and variance in space domain for fit
% data. normalize grating and natural image sets to match.
mm=0;
ss=0;
flentot=0;
pteststim=zeros(size(nrmov(:,:,1)));

nstimloadparms=params.stimloadparms;
nstimfilterparms=params.stimfilterparms;

mstim=zeros(Xmax);
ttot=0;
for fidx=1:filecount,
   
   % deal with online stim aperture sizing
   if (strcmp(params.stimloadcmd,'loadimfile') |...
       strcmp(params.stimloadcmd,'loadimscaled')) & ...
         length(nstimloadparms)>0 & nstimloadparms{1}==0,
      nstimloadparms{1}=round(nstimloadparms{3} .* ...
                              params.stimcrfs(fidx)./params.stimwindowcrf);
   end
   
   flen=imfileinfo(params.stimfiles{fidx},1);
   if flen>(4000./filecount) & ismember(params.batch,[77,78,89]),
      flen=round(4000./filecount);
   elseif flen>20000./filecount,
      flen=round(20000./filecount);
   end
   
   teststim=feval(params.stimloadcmd,params.stimfiles{fidx},...
                  1,flen,nstimloadparms{:});
   mstim=mstim+sum(teststim,3);
   ttot=ttot+sum(teststim(Xmax/2,Xmax/2,:));
   if DOSPECTNORM,
      ptemp=fft(fft(teststim,[],1),[],2);
      ptemp=sqrt(mean(abs(ptemp).^2,3));
      %ptemp=(mean(abs(ptemp),3));
      pteststim=pteststim+ptemp.*params.resplen(fidx);
   end
   
   teststim=reshape(teststim,Xmax*Xmax,size(teststim,3));
   mm=mm+mean(teststim(:)).*params.resplen(fidx);
   mpf=mean(teststim,1);
   sstd=std(teststim);
   sspf=(mean((teststim-repmat(mpf,Xmax*Xmax,1)).^2));
   spf=sqrt(mean((teststim-repmat(mpf,Xmax*Xmax,1)).^2))./(mpf+(mpf==0));
   ss=ss+mean(sstd) .* params.resplen(fidx);
end
mm=mm./sum(params.resplen);
ss=(ss./sum(params.resplen));
pteststim=pteststim./sum(params.resplen);

fprintf('fit set stats: mm=%.1f ss=%.1f\n',mm,ss);

mstim=mstim./ttot;


% currently testing several different means of normalizing natural
% image and grating sets
fprintf('normalizing test stimuli (normeachframe=%d)\n',NORMEACHFRAME);

if DOSPECTNORM,
   fprintf('DOSPECTNORM=%d\n',DOSPECTNORM);

   mtemp=fft(fft(nrmov,[],1),[],2);
   mnrmovtest=sqrt(mean(abs(mtemp).^2,3));
   %mnrmovtest=(mean(abs(mtemp),3));
   spectadj=pteststim./mnrmovtest;
   mtemp=mtemp.*repmat(spectadj,[1 1 nrmovlen]);
   mtemp=real(ifft(ifft(mtemp,[],1),[],2));
   nrmovspace=mtemp;
   
   ts=size(cart);
   tso=[Xmax Xmax prod(ts(2:end))];
   mtemp=fft(fft(reshape(cart,tso),[],1),[],2);
   mcartpow=sqrt(mean(abs(mtemp).^2,3));
   spectadj=pteststim./mcartpow;
   mtemp=mtemp.*repmat(spectadj,[1 1 tso(3)]);
   mtemp=real(ifft(ifft(mtemp,[],1),[],2));
   cart=reshape(mtemp,ts);
   
   ts=size(pol);
   tso=[Xmax Xmax prod(ts(2:end))];
   mtemp=fft(fft(reshape(pol,tso),[],1),[],2);
   mpolpow=sqrt(mean(abs(mtemp).^2,3));
   spectadj=pteststim./mpolpow;
   mtemp=mtemp.*repmat(spectadj,[1 1 tso(3)]);
   mtemp=real(ifft(ifft(mtemp,[],1),[],2));
   pol=reshape(mtemp,ts);
   
   ts=size(hyp);
   tso=[Xmax Xmax prod(ts(2:end))];
   mtemp=fft(fft(reshape(hyp,tso),[],1),[],2);
   mhyppow=sqrt(mean(abs(mtemp).^2,3));
   spectadj=pteststim./mhyppow;
   mtemp=mtemp.*repmat(spectadj,[1 1 tso(3)]);
   mtemp=real(ifft(ifft(mtemp,[],1),[],2));
   hyp=reshape(mtemp,ts);
   
   clear mtemp
else
   nrmovspace=nrmov;
end


% NORMEACHFRAME==2: don't do any normalization?
if NORMEACHFRAME==1,
   
   nrn=(0:Xmax/2+1)+round(Xmax/4);
   nrnl=length(nrn);
   nmm=repmat(mean(reshape(nrmovspace(nrn,nrn,:),nrnl^2,1,nrmovlen)),Xmax,Xmax);
   nss=repmat(std(reshape(nrmovspace(nrn,nrn,:),nrnl^2,1,nrmovlen)),Xmax,Xmax);
   nss(nss==0)=1;
   
   % normalize mean & variance of each frame individually
   % extra round of normalization to deal with hitting grayscale rails
   %nrmovspace=((nrmovspace-nmm)./nss).*ss+mm;
   %nrmovspace(nrmovspace<0)=0;nrmovspace(nrmovspace>255)=255;
   %nmm=repmat(mean(reshape(nrmovspace,Xmax*Xmax,1,nrmovlen)),Xmax,Xmax);
   %nss=repmat(std(reshape(nrmovspace,Xmax*Xmax,1,nrmovlen)),Xmax,Xmax);
   
   cmm=mean(cart(:));
   pmm=mean(pol(:));
   hmm=mean(hyp(:));
   css=mean(std(cart(:,:)));
   pss=mean(std(pol(:,:)));
   hss=mean(std(hyp(:,:)));
   %cmm=repmat(mean([cart(:,:)]),Xmax*Xmax,1);
   %pmm=repmat(mean([pol(:,:)]),Xmax*Xmax,1);
   %hmm=repmat(mean([hyp(:,:)]),Xmax*Xmax,1);
   %css=repmat(std([cart(:,:)]),Xmax*Xmax,1);
   %pss=repmat(std([pol(:,:)]),Xmax*Xmax,1);
   %hss=repmat(std([hyp(:,:)]),Xmax*Xmax,1);
   
   % do the normalization for each frame
   nrmovspace=((nrmovspace-nmm)./nss).*ss+mm;
   cart(:,:)=((cart(:,:)-cmm)./css).*ss+mm;
   pol(:,:)=((pol(:,:)-pmm)./pss).*ss+mm;
   hyp(:,:)=((hyp(:,:)-hmm)./hss).*ss+mm;
   
   if DOEDSTIM,
      emm=mean(edmov(:));
      ess=mean(std(edmov(:,:)));
      edmov(:,:)=((edmov(:,:)-emm)./ess).*ss+mm;
   end
end


if 1,
   % do this no matter what since NORMEACH GETS screwed up by
   % splitting groups
   
   nmm=mean(reshape(nrmovspace(nrn,nrn,:),nrnl^2*nrmovlen,1));
   nss=mean(std(reshape(nrmovspace(nrn,nrn,:),nrnl^2,nrmovlen)));
   
   cmm=mean(cart(:));
   pmm=mean(pol(:));
   hmm=mean(hyp(:));
   css=mean(std(cart(:,:)));
   pss=mean(std(pol(:,:)));
   hss=mean(std(hyp(:,:)));
   if DOEDSTIM,
      emm=mean(edmov(:));
      ess=mean(std(edmov(:,:)));
   end
end

if 1,
   % adjust mm and ss of carts to make them vary as much as natural images
   disp('splitting into high and low contrast groups');

   meanadj=[1 1   2 2];
   stdadj= [1 1.4 1.0 1.4]; 
   %meanadj=[1 2];
   %stdadj= [1 1];
   %meanadj=[1 1 1 1];
   %stdadj= [1 0.75 1.5 2.25]; 
   %meanadj=[1 1.1 1.2 1.3];
   %stdadj= [1 1 1 1];
   sc=length(meanadj);
   
   % double the length of the test natural image set up to 20000 frames
   nlenmax=5000;
   if nrmovlen>nlenmax,
      newidx=round(linspace(1,min([nrmovlen 12000]),nlenmax));
      nrmovspace=nrmovspace(:,:,newidx);
      nrmovlen=nlenmax;
   end
   nok=repmat(nok(newidx),[1 sc]);
   nrmovspace=repmat(nrmovspace,[1 1 1 sc]);
   nrmovlen=nrmovlen*sc;
   if DOEDSTIM,
      edmov=repmat(edmov,[1 1 1 sc]);
   end
   
   for ii=2:length(meanadj),
      nrmovspace(:,:,:,ii)=...
          (nrmovspace(:,:,:,ii)-nmm).*stdadj(ii)+nmm*meanadj(ii);
      if DOEDSTIM,
         edmov(:,:,:,ii)=(edmov(:,:,:,ii)-emm).*stdadj(ii)+emm*meanadj(ii);
      end
   end
   nrmovspace=reshape(nrmovspace,Xmax,Xmax,nrmovlen);
   if DOEDSTIM,
      edmovlen=edmovlen*sc;
      edmov=reshape(edmov,Xmax,Xmax,edmovlen);
      emm=mean(edmov(:));
      ess=mean(std(edmov(:,:)));
   end
   
   % create four reps along the "phase" dimension so that it
   % multiplexes phase,luminance and contrast
   dcount3=dcount3*sc;
   
   for ii=2:length(meanadj),
      cart(:,:,:,:,ii)=(cart(:,:,:,:,1)-cmm)*stdadj(ii)+cmm*meanadj(ii);
      pol(:,:,:,:,ii)=(pol(:,:,:,:,1)-pmm)*stdadj(ii)+pmm*meanadj(ii);
      hyp(:,:,:,:,ii)=(hyp(:,:,:,:,1)-hmm)*stdadj(ii)+hmm*meanadj(ii);
   end
   cart=reshape(cart,Xmax*Xmax,dcount,dcount2,dcount3);
   cartprange=repmat(cartprange,[1 4]);
   
   pol=reshape(pol,Xmax*Xmax,dcount,dcount2,dcount3);
   polprange=repmat(polprange,[1 4]);
   
   hyp=reshape(hyp,Xmax*Xmax,dcount,dcount2,dcount3);
   hypprange=repmat(hypprange,[1 4]);
   
   % recalc mean and std across expanded stim sets
   nmm=mean(reshape(nrmovspace(nrn,nrn,:),nrnl^2*nrmovlen,1));
   nss=mean(std(reshape(nrmovspace(nrn,nrn,:),nrnl^2,nrmovlen)));
   
   cmm=mean(cart(:));
   pmm=mean(pol(:));
   hmm=mean(hyp(:));
   css=mean(std(cart(:,:)));
   pss=mean(std(pol(:,:)));
   hss=mean(std(hyp(:,:)));
end

if NORMEACHFRAME<2,
   % normalize test stimuli
   nrmovspace=((nrmovspace-nmm)./nss).*ss+mm;
   cart(:,:)=((cart(:,:)-cmm)./css).*ss+mm;
   pol(:,:)=((pol(:,:)-pmm)./pss).*ss+mm;
   hyp(:,:)=((hyp(:,:)-hmm)./hss).*ss+mm;
   if DOEDSTIM,
      edmov(:,:)=((edmov(:,:)-emm)./ess).*ss+mm;
   end
end

% threshold to match grayscale limits
nrmovspace(nrmovspace<0)=0;nrmovspace(nrmovspace>255)=255;
cart(cart<0)=0; cart(cart>255)=255;
pol(pol<0)=0; pol(pol>255)=255;
hyp(hyp<0)=0; hyp(hyp>255)=255;

nmm=mean(reshape(nrmovspace(nrn,nrn,:),nrnl^2*nrmovlen,1));
nss=mean(std(reshape(nrmovspace(nrn,nrn,:),nrnl^2,nrmovlen)));
fprintf('final  nat: nmm=%.1f nss=%.1f\n',nmm,nss);
fprintf('final cart: cmm=%.1f css=%.1f\n',...
        mean(cart(:)),mean(std(cart(:,:))));
fprintf('final  pol: pmm=%.1f pss=%.1f\n',...
        mean(pol(:)),mean(std(pol(:,:))));
fprintf('final  hyp: hmm=%.1f hss=%.1f\n',...
        mean(hyp(:)),mean(std(hyp(:,:))));
if DOEDSTIM,
   edmov(edmov<0)=0; edmov(edmov>255)=255;
   fprintf('final  edm: emm=%.1f ess=%.1f\n',...
           mean(edmov(:)),mean(std(edmov(:,:))));
end

if 0,
   hist(mean(reshape(nrmovspace,Xmax*Xmax,nrmovlen)))
   hist(std(reshape(nrmovspace,Xmax*Xmax,nrmovlen)))
   
   hist(mean(cart(:,:))) 
   hist(mean(pol(:,:))) 
   hist(mean(hyp(:,:))) 
   hist(std(cart(:,:))) 
   hist(std(pol(:,:))) 
   hist(std(hyp(:,:))) 
end

% convert normalized test stimuli into linearized domain for doing
% predictions with the SRF

% batch 82,99 - bg=46
% batch 77 - bg=20

if ~isempty(params.stimfiltercmd),
   if batch==77,
      mstim=ones(size(mstim));
      mz=0;
   else
      % kludgy step to match window to the window on the fit stim
      mstim=mean(cat(3,mstim,flipud(mstim),fliplr(mstim),mstim'),3);
      mz=mstim(1,1);
      mstim=(mstim-mz)./max(mstim(:)-mz);
      mz=teststim(1);
   end
   
   nrmov=feval(params.stimfiltercmd,...
               (nrmovspace-mz).*repmat(mstim,[1 1 size(nrmovspace,3)])+mz,...
               nstimfilterparms{:});
   
   ts=[size(cartpfft,1) dcount dcount2 dcount3];
   tso=[Xmax Xmax prod(ts(2:end))];
   cartpfft=feval(params.stimfiltercmd,...
                  reshape(cart-mz,tso).*repmat(mstim,[1 1 tso(3)])+mz,...
                  nstimfilterparms{:});
   cartpfft=reshape(cartpfft,ts);
   polpfft=feval(params.stimfiltercmd,...
                 reshape(pol-mz,tso).*repmat(mstim,[1 1 tso(3)])+mz,...
                 nstimfilterparms{:});
   polpfft=reshape(polpfft,ts);
   hyppfft=feval(params.stimfiltercmd,...
                 reshape(hyp-mz,tso).*repmat(mstim,[1 1 tso(3)])+mz,...
                 nstimfilterparms{:});
   hyppfft=reshape(hyppfft,ts);
   
   if DOEDSTIM,
      edpfft=feval(params.stimfiltercmd,edmov,nstimfilterparms{:});
   else
      edpfft=[];
      edmovlen=0;
   end
end

if 0,
   figure(3);
   clf
   plot(mean(nrmov,2));
   hold on
   plot(mS(:,1),'k--');
   plot(mean(polpfft(:,:),2),'g');
   hold off
end

% subtract mean from pfft stim, since that's always necessary for
% any pred/corr thing:
nrmov=nrmov-repmat(mS(:,1),[1 nrmovlen]);
ts=[size(cartpfft,1) dcount dcount2 dcount3];
cartpfft=cartpfft-repmat(mS(:,1),[1 ts(2:end)]);
polpfft=polpfft-repmat(mS(:,1),[1 ts(2:end)]);
hyppfft=hyppfft-repmat(mS(:,1),[1 ts(2:end)]);
edpfft=edpfft-repmat(mS(:,1),[1 edmovlen]);

% predict response to each grating stimulus
cartmatch=zeros(dcount,dcount2,dcount3,params.bootcount);
polmatch=zeros(dcount,dcount2,dcount3,params.bootcount);
hypmatch=zeros(dcount,dcount2,dcount3,params.bootcount);

cct=size(H,1);
%Hnorm=mean(H(:,:,1),2);
Hnorm=H(:,:,1);
Hnorm=(Hnorm-repmat(mean(Hnorm),[cct 1]));

HNcount=size(Hnorm,2);

cartmatchcc=zeros(dcount*dcount2*dcount3,HNcount);
polmatchcc=zeros(dcount*dcount2*dcount3,HNcount);
hypmatchcc=zeros(dcount*dcount2*dcount3,HNcount);
natmatchcc=zeros(nrmovlen,HNcount);
edmatchcc=zeros(edmovlen,HNcount);

fprintf('computing ccs');
for bootidx=1:HNcount,
   fprintf('.');
   if sum(abs(Hnorm(:,bootidx)))>0,
      Hnorm(:,bootidx)=Hnorm(:,bootidx)./sqrt(Hnorm(:,bootidx)'*Hnorm(:,bootidx));
   
      cartmatchcc(:,bootidx)=...
          ((cartpfft(:,:)-repmat(mean(cartpfft(:,:)),[cct 1]))'*Hnorm(:,bootidx))./...
          sqrt(sum((cartpfft(:,:)-repmat(mean(cartpfft(:,:)),[cct 1])).^2)');
      polmatchcc(:,bootidx)=...
          ((polpfft(:,:)-repmat(mean(polpfft(:,:)),[cct 1]))'*Hnorm(:,bootidx))./...
          sqrt(sum((polpfft(:,:)-repmat(mean(polpfft(:,:)),[cct 1])).^2)');
      hypmatchcc(:,bootidx)=...
          ((hyppfft(:,:)-repmat(mean(hyppfft(:,:)),[cct 1]))'*Hnorm(:,bootidx))./...
          sqrt(sum((hyppfft(:,:)-repmat(mean(hyppfft(:,:)),[cct 1])).^2)');
      natmatchcc(:,bootidx)=...
          (nrmov-repmat(mean(nrmov),[cct 1]))'*Hnorm(:,bootidx) ./ ...
          sqrt(sum((nrmov-repmat(mean(nrmov),[cct 1])).^2)');
      edmatchcc(:,bootidx)=...
          (edpfft-repmat(mean(edpfft),[cct 1]))'*Hnorm(:,bootidx) ./ ...
          sqrt(sum((edpfft-repmat(mean(edpfft),[cct 1])).^2)');
   end
end
fprintf('\n');
cartmatchcc=reshape(cartmatchcc,dcount,dcount2,dcount3,HNcount);
polmatchcc=reshape(polmatchcc,dcount,dcount2,dcount3,HNcount);
hypmatchcc=reshape(hypmatchcc,dcount,dcount2,dcount3,HNcount);

fprintf('testing grat resp...\n');
for d1=1:dcount,
   for d2=1:dcount2,
      for d3=1:dcount3,
         ctest=cartpfft(:,d1,d2,d3);
         ptest=polpfft(:,d1,d2,d3);
         htest=hyppfft(:,d1,d2,d3);
         
         for bootidx=1:params.bootcount,
            cartmatch(d1,d2,d3,bootidx)=(ctest)'*H(:,bootidx,1);
            polmatch(d1,d2,d3,bootidx)=(ptest)'*H(:,bootidx,1);
            hypmatch(d1,d2,d3,bootidx)=(htest)'*H(:,bootidx,1);
         end
      end
   end
end
cartmatch=cartmatch+repmat(reshape(nlparms(1,:,1),[1 1 1 params.bootcount]),...
                           [dcount,dcount2,dcount3 1]);
polmatch=polmatch+repmat(reshape(nlparms(1,:,1),[1 1 1 params.bootcount]),...
                         [dcount,dcount2,dcount3 1]);
hypmatch=hypmatch+repmat(reshape(nlparms(1,:,1),[1 1 1 params.bootcount]),...
                         [dcount,dcount2,dcount3 1]);

% predict response to each natural image
natmatch=zeros(nrmovlen,params.bootcount);
natmatchpos=zeros(nrmovlen,params.bootcount);
natmatchneg=zeros(nrmovlen,params.bootcount);
edmatch=zeros(edmovlen,params.bootcount);

fprintf('testing nat resp...\n');
for bootidx=1:params.bootcount,
   natmatch(:,bootidx)=nrmov'*H(:,bootidx,1);
   
   Hp=H(:,bootidx,1).*(H(:,bootidx,1)>0);
   Hn=H(:,bootidx,1).*(H(:,bootidx,1)<0);
   natmatchpos(:,bootidx)=(nrmov)'*Hp;
   natmatchneg(:,bootidx)=(nrmov)'*Hn;

   edmatch(:,bootidx)=edpfft'*H(:,bootidx,1);
end
natmatch=natmatch+repmat(nlparms(1,:,1),[nrmovlen 1]);
edmatch=edmatch+repmat(nlparms(1,:,1),[edmovlen 1]);

if (std(natmatchpos)+std(natmatchneg))>0,
   pnrat=(std(natmatchpos)-std(natmatchneg))./...
         (std(natmatchpos)+std(natmatchneg));
else
   pnrat=zeros(size(std(natmatchpos)+std(natmatchneg)));
end

% this is an old, unsuccessful attempt at more realstic mean firing
% rate measurements.
if USERESFILE,
   tr=[reshape(cartmatch,dcount*dcount2*dcount3,params.bootcount);
       reshape(polmatch,dcount*dcount2*dcount3,params.bootcount);
       reshape(hypmatch,dcount*dcount2*dcount3,params.bootcount);
       natmatch];
   scale2Hz=repmat(72,[1 size(tr,2)]); % meanactresp./mean(tr);
   clear tr
end

% adjust spike rate to something resembling actual rates. hopefully
% doing so in a way that doesn't create big variance between
% jackknifed strfs
for bootidx=1:params.bootcount,
   cartmatch(:,:,:,bootidx)=cartmatch(:,:,:,bootidx).*scale2Hz(bootidx);
   polmatch(:,:,:,bootidx)=polmatch(:,:,:,bootidx).*scale2Hz(bootidx);
   hypmatch(:,:,:,bootidx)=hypmatch(:,:,:,bootidx).*scale2Hz(bootidx);
   natmatch(:,bootidx)=natmatch(:,bootidx).*scale2Hz(bootidx);
   edmatch(:,bootidx)=edmatch(:,bootidx).*scale2Hz(bootidx);
end

% figure out mean and stderr on each of the predicted responses
if 1
   cartmatch0=cartmatch;
   polmatch0=polmatch;
   hypmatch0=hypmatch;
   natmatch0=natmatch;
end

cartmatchm=mean(mean(mean(cartmatch,1),2),3);
polmatchm=mean(mean(mean(polmatch,1),2),3);
hypmatchm=mean(mean(mean(hypmatch,1),2),3);
natmatchm=mean(natmatch,1);
edmatchm=mean(edmatch,1);
cartmatcherr=std(cartmatch-repmat(cartmatchm,[dcount,dcount2,dcount3]),...
                 1,4) .*sqrt(params.bootcount-1);
polmatcherr=std(polmatch-repmat(polmatchm,[dcount,dcount2,dcount3]),...
                1,4).*sqrt(params.bootcount-1);
hypmatcherr=std(hypmatch-repmat(polmatchm,[dcount,dcount2,dcount3]),...
                1,4).*sqrt(params.bootcount-1);
natmatcherr=std(natmatch-repmat(natmatchm,[size(natmatch,1) 1]),...
                1,2).*sqrt(params.bootcount-1);
edmatcherr=std(edmatch-repmat(edmatchm,[size(edmatch,1) 1]),...
                1,2).*sqrt(params.bootcount-1);
cartmatch=mean(cartmatch,4);
polmatch=mean(polmatch,4);
hypmatch=mean(hypmatch,4);
natmatch=mean(natmatch,2);
edmatch=mean(edmatch,2);

cartmatchccm=mean(mean(mean(cartmatchcc,1),2),3);
polmatchccm=mean(mean(mean(polmatchcc,1),2),3);
hypmatchccm=mean(mean(mean(hypmatchcc,1),2),3);
natmatchccm=mean(natmatchcc,1);
edmatchccm=mean(edmatchcc,1);
cartmatchccerr=std(cartmatchcc-repmat(cartmatchccm,[dcount,dcount2,dcount3]),...
                 1,4) .*sqrt(params.bootcount-1);
polmatchccerr=std(polmatchcc-repmat(polmatchccm,[dcount,dcount2,dcount3]),...
                 1,4).*sqrt(params.bootcount-1);
hypmatchccerr=std(hypmatchcc-repmat(polmatchccm,[dcount,dcount2,dcount3]),...
                1,4).*sqrt(params.bootcount-1);
natmatchccerr=std(natmatchcc-repmat(natmatchccm,[size(natmatchcc,1) 1]),...
                1,2).*sqrt(params.bootcount-1);
edmatchccerr=std(edmatchcc-repmat(edmatchccm,[size(edmatchcc,1) 1]),...
                1,2).*sqrt(params.bootcount-1);
if 1
   cartmatchcc0=cartmatchcc;
   polmatchcc0=polmatchcc;
   hypmatchcc0=hypmatchcc;
   natmatchcc0=natmatchcc;
end
cartmatchcc=median(cartmatchcc,4);
polmatchcc=median(polmatchcc,4);
hypmatchcc=median(hypmatchcc,4);
natmatchcc=median(natmatchcc,2);
edmatchcc=median(edmatchcc,2);

% compute entropy stuff
mmmin=min([cartmatch(:); polmatch(:); hypmatch(:); natmatch(:); edmatch(:)]);
mmmax=max([cartmatch(:); polmatch(:); hypmatch(:); natmatch(:); edmatch(:)]);

bins=floor(mmmin):ceil(mmmax);
%bins=linspace(min(natmatch(:)),max(natmatch(:)),30+1);
nn=hist(cartmatch(:),bins);
nn=nn./sum(nn);
hcart=-sum(nn.*log2(nn+(nn==0)));
nn=hist([polmatch(:);hypmatch(:)],bins);
nn=nn./sum(nn);
hpol=-sum(nn.*log2(nn+(nn==0)));
nn=hist([cartmatch(:);polmatch(:);hypmatch(:)],bins);
nn=nn./sum(nn);
hhyp=-sum(nn.*log2(nn+(nn==0)));
nn=hist(natmatch(:),bins);
nn=nn./sum(nn);
hnat=-sum(nn.*log2(nn+(nn==0)));

%
% figure out peak responses to each stim class

bootidx=1;

cm=cartmatch(:,:,:,bootidx);
cmin=min(find(cm==min(cm(:))));
cmax=min(find(cm==max(cm(:))));
[cm1,cm2,cm3]=ind2sub([dcount dcount2 dcount3],cmax);

pm=polmatch(:,:,:,bootidx);
pmin=min(find(pm==min(pm(:))));
pmax=min(find(pm==max(pm(:))));
[pm1,pm2,pm3]=ind2sub([dcount dcount2 dcount3],pmax);

hm=hypmatch(:,:,:,bootidx);
hmin=min(find(hm==min(hm(:))));
hmax=min(find(hm==max(hm(:))));
[hm1,hm2,hm3]=ind2sub([dcount dcount2 dcount3],hmax);

nmin=min(find(natmatch(:,bootidx)==min(natmatch(:,bootidx))));
nmax=min(find(natmatch(:,bootidx)==max(natmatch(:,bootidx))));

[dummy,nmset]=sort(-natmatch(:,bootidx));
noptsetlen=30;
nmset=nmset(1:noptsetlen);
if exist('prefframe','var'),
   fprintf('prefframe=%d\n',prefframe);
else
   prefframe=1;
end

nm1=nmset(prefframe);

[edmax,edset]=sort([-edmatch(:,bootidx) edmatcherr(:)]);

edmaxerr=edmax(1,2);
edmax=-edmax(1,1);

edsetlen=30;
edset=edset(1:edsetlen);
ed1=edset(1);

edocount=8;
edsepcount=7;
edpolcount=3;
[em1,em2,em3]=ind2sub([edocount edsepcount edpolcount*sc],ed1);


[edccmax,edcc1]=sortrows([-edmatchcc(:) edmatchccerr(:)]);
edccmaxerr=edccmax(1,2);
edccmin=-edccmax(end,1);
edccmax=-edccmax(1,1);

[ecc1,ecc2,ecc3]=ind2sub([edocount edsepcount edpolcount*4],edcc1(1));

edmatch=reshape(edmatch,[edocount edsepcount edpolcount*4]);

[cartmax ii]=sortrows([cartmatch(:) cartmatcherr(:)]);
%cartmaxerr=median(cartmax(end-3:end,1));
%cartmax=median(cartmax(end-3:end,1));
cartmaxerr=cartmax(end,2);
cartmax=cartmax(end,1);
cmfull=squeeze(cartmatch0(cm1,cm2,cm3,:));

[ncartmax ii]=sortrows([polmatch(:) polmatcherr(:); ...
                    hypmatch(:) hypmatcherr(:)]);
%ncartmaxerr=median(ncartmax(end-7:end,2));
%ncartmax=median(ncartmax(end-7:end,1));
ncartmaxerr=median(ncartmax(end-1:end),2);
ncartmax=median(ncartmax(end-1:end,1));
if ii(end)>length(polmatch),
   nmfull=squeeze(hypmatch0(hm1,hm2,hm3,:));
else
   nmfull=squeeze(polmatch0(pm1,pm2,pm3,:));
end

[cartccmax ii]=sortrows([cartmatchcc(:) cartmatchccerr(:)]);
%cartccmaxerr=median(cartccmax(end-3:end,2));
%cartccmax=median(cartccmax(end-3:end,1));
cartccmaxerr=cartccmax(end,2);
cartccmin=cartccmax(1,1);
cartccmax=cartccmax(end,1);
cmccfull=squeeze(cartmatchcc0(cm1,cm2,cm3,:));

ncartccmax=[polmatchcc(:) polmatchccerr(:); hypmatchcc(:) hypmatchccerr(:)];
[ncartccmax ii]=sort(ncartccmax);
%ncartccmax=median(ncartccmax(end-7:end));
%ncartccmaxerr=median(ncartccmaxerr(ncrank(end-7:end)));
ncartccmaxerr=median(ncartccmax(end-1:end,2));
ncartccmin=median(ncartccmax(1:2,1));
ncartccmax=median(ncartccmax(end-1:end,1));
if ii(end)>length(polmatchcc),
   nmccfull=squeeze(hypmatchcc0(hm1,hm2,hm3,:));
else
   nmccfull=squeeze(polmatchcc0(pm1,pm2,pm3,:));
end


% 1/cartcount picks inner edge of outer bin for natural stim responses
cartcount=length(cartmatch(:));
gratcount=length([cartmatch(:); polmatch(:); hypmatch(:)]);

%cartpct=1./cartcount *4;  % since the four phases are kind of fake
cartpct=1./cartcount;
gratpct=1./gratcount;

fprintf('cartcount=%d, cartpct=%.3f\n',cartcount,cartpct);

[natccmax ii]=sortrows([natmatchcc(:) natmatchccerr(:)]);
natccmaxerr=median(natccmax(round(length(natmatchcc).*(1-cartpct)):end,2));
natccmin=median(natccmax(1:round(length(natmatchcc).*cartpct),1));
natccmax=median(natccmax(round(length(natmatchcc).*(1-cartpct)):end,1));
natmccfull=natmatchcc0(ii(end),:)';

synthmax=[cartmatch(:) cartmatcherr(:); polmatch(:) polmatcherr(:); ...
          hypmatch(:) hypmatcherr(:)];
[synthmax,ncrank]=sortrows(synthmax);
synthmaxerr=synthmax(end,2);
synthmax=synthmax(end,1);

natfracbetter=mean(natmatch>synthmax);
natovercart=mean(natmatch>cartmatch(cmax));

[nsort ii]=sortrows([natmatch natmatcherr]);
natmin1=median(nsort(1:round(length(nsort).*cartpct),1));
natmax1=median(nsort(round(length(nsort).*(1-cartpct)):end,1));
natmax1err=median(nsort(round(length(nsort).*(1-cartpct)):end,2));
natmin2=median(nsort(1:round(length(nsort).*gratpct),1));
natmax2=median(nsort(round(length(nsort).*(1-gratpct)):end,1));
natmax2err=median(nsort(round(length(nsort).*(1-gratpct)):end,2));
natmaxerr=nsort(end,2);
natmfull=natmatch0(ii(end),:)';

testfull=[cmfull nmfull natmfull];
mmm=mean(testfull,2);
testfull=testfull./repmat(mmm,[1 size(testfull,2)]);
testfullerr=std(testfull,1,1).*sqrt(params.bootcount-1).*mean(mmm)
cartmaxerr=testfullerr(1);
ncartmaxerr=testfullerr(2);
natmax1err=testfullerr(3);

testfull=[cmccfull nmccfull natmccfull];
mmm=nanmean(testfull')';
testfull=testfull./repmat((mmm+(mmm==0)),[1 size(testfull,2)]);
testfullerr=std(testfull,1,1).*sqrt(params.bootcount-1).*mean(mmm)
cartccmaxerr=testfullerr(1);
ncartccmaxerr=testfullerr(2);
natccmaxerr=testfullerr(3);

fprintf(' cart max: %5.1f +/- %.1f  cc: %.2f +/- %.2f\n',...
        cartmax,cartmaxerr,cartccmax,cartccmaxerr);
fprintf('ncart max: %5.1f +/- %.1f  cc: %.2f +/- %.2f\n',...
        ncartmax,ncartmaxerr,ncartccmax,ncartccmaxerr);
if DOEDSTIM,
   fprintf('edmov max: %5.1f +/- %.1f  cc: %.2f +/- %.2f\n',...
           edmax,edmaxerr,edccmax,edccmaxerr);
end
fprintf('  nat max: %5.1f +/- %.1f  cc: %.2f +/- %.2f\n',...
        natmax1,natmax1err,natccmax,natccmaxerr);


titles{1}=sprintf('%s: m %.1f Hz ncc=%.0f (%.0f) crt=%.0f nc=%.0f',...
                  cellid,meanactresp,...
                  acos(natccmax).*180/pi,acos(max(natmatchcc)).*180/pi,...
                  acos(cartccmax).*180/pi,acos(ncartccmax).*180/pi);
titles{2}=sprintf('cart: %.1f +%.1f (%.1f,%.1f,%.0f)',...
                  cartmax,cartmaxerr,...
                  cartfrange(cm1),cartorange(cm2),cartprange(cm3));
titles{3}=sprintf('polar: %.1f +%.1f (%.0f,%.2f,%.0f)',...
                  polmatch(pm1,pm2,pm3,bootidx),...
                  polmatcherr(pm1,pm2,pm3,bootidx),...
                  polprange(pm3),polcrange(pm1),polrrange(pm2));
titles{4}=sprintf('hyp: %.1f +%.1f (%.0f,%.0f,%.2f)',...
                  hypmatch(hm1,hm2,hm3,bootidx),...
                  hypmatcherr(hm1,hm2,hm3,bootidx),...
                  hyporange(hm1),hypprange(hm3),hypfrange(hm2));
titles{5}=sprintf('nat: %.1f +%.1f (%d,%.3f)',...
                  natmax1,natmax1err,nm1,natfracbetter);
titles{6}=sprintf('ed: %.1f +%.1f (%d,%d,%d)',edmax,edmaxerr,em1,em2,em3);

% fourier domain optimal gratings
fcopt=cartpfft(:,cm1,cm2,cm3);
fpopt=polpfft(:,pm1,pm2,pm3);
fhopt=hyppfft(:,hm1,hm2,hm3);
fnopt=nrmov(:,nm1);
fnset=nrmov(:,nmset);
if DOEDSTIM,
   feopt=edpfft(:,edset(1));
   feset=edpfft(:,edset);
else
   feopt=[];
end

%mean(H(:,:,1),2)
tk=cat(3,mean(H(:,:,1),2),fcopt,fpopt,fhopt,fnopt,feopt);

% space domain optimal gratings at high resolution
if BIGNCARTS,
   copt=reshape(cart(:,cm1,cm2,cm3),Xmax,Xmax)./255;
   popt=reshape(pol(:,cm1,cm2,cm3),Xmax,Xmax)./255;
   hopt=reshape(hyp(:,cm1,cm2,cm3),Xmax,Xmax)./255;
else
   
   [t,copt]=ncart_pfft([cartfrange(cm1),cartorange(cm2),cartprange(cm3)],...
      Xmax*4,'cart',params.stimfiltercmd,params.stimfilterparms,[.5 .3]);
   [t,popt]=ncart_pfft([polprange(pm3),polcrange(pm1)./4,polrrange(pm2)],...
      Xmax.*4,'pol',params.stimfiltercmd,params.stimfilterparms,[.5 .3]);
   [t,hopt]=ncart_pfft([hyporange(hm1),hypprange(hm3),hypfrange(hm2)./4],...
      Xmax.*4,'hyp',params.stimfiltercmd,params.stimfilterparms,[.5 .3]);
   
   copt(find(copt>1))=1; copt(find(copt<0))=0;
   popt(find(popt>1))=1; popt(find(popt<0))=0;
   hopt(find(hopt>1))=1; hopt(find(hopt<0))=0;
end

noptset=loadimframes(natrevfile,nok(nmset),96,0,64);

if NORMEACHFRAME,
   nrnl2=size(noptset,1);
   nmm=repmat(mean(reshape(noptset,nrnl2*nrnl2,1,noptsetlen)),nrnl2,nrnl2);
   nss=repmat(std(reshape(noptset,nrnl2*nrnl2,1,noptsetlen)),nrnl2,nrnl2);
   nss(find(nss==0))=1;
else
   % keep nmm,nss same as above
end

noptset=(noptset-nmm)./nss;
noptset=(noptset.*0.2)+0.5;
noptset(noptset<0)=0;
noptset(noptset>1)=1;
nopt=noptset(:,:,prefframe);
if DOEDSTIM,
   eopt=edmov(:,:,edset(1))./255;
end

% display optimal grating from each class and optimal natural image
figure(1);
clf
showkern(repmat(tk,[1,2]),kernfmt,[Xmax Xmax],titles,1);

subplot(6,2,1);

if 1,
   bbdata=[cartmax ncartmax natmax1];
   bbdataerr=[cartmaxerr ncartmaxerr natmax1err];
   errorbar(bbdata,bbdataerr,'k+');
   hold on
   h=bar(bbdata);
   hold off
   set(h,'FaceColor',[0 0 0]);
   
   axis square
   
elseif 1
   [x,y]=arc(0,0,0,90,1);
   plot(x,y,'k-','LineWidth',1);
   hold on
   plot([0 0],[0 1],'k-','LineWidth',1);
   plot([0 sqrt(1-natccmax.^2)],[0 natccmax],'k-','LineWidth',2);
   plot([0 sqrt(1-max(natmatchcc).^2)],[0 max(natmatchcc)],'k--','LineWidth',1);
   plot([0 sqrt(1-cartccmax.^2)],[0 cartccmax],'k-','LineWidth',2);
   plot([0 sqrt(1-ncartccmax.^2)],[0 ncartccmax],'k-','LineWidth',2);
   hold off
   axis square
   axis off
else
   tsfe=pfft2sf(H(:,1),params.kernfmt);
   tsfo=tsfe;
   tsfo(cfilt)=i.*tsfo(cfilt);
   tsfo(cfiltconj)=-i.*tsfo(cfiltconj);
   
   fhann=fftshift(abs(fft2(hanning2(Xmax))).^2);
   tsfo=cconv2(fhann,tsfo,1);
   tsfe=cconv2(fhann,tsfe,1);
   
   ospace=fftshift(real(ifft2(fftshift(tsfo))));
   espace=fftshift(real(ifft2(fftshift(tsfe))));
   hspace=[ospace espace];
   hspace=hspace-min(hspace(:));
   hspace=hspace./max(hspace(:));
   
   %imagesc(repmat(hspace,[1 1 3])); axis image; axis off;
end   
   
title(titles{1});
   
subplot(6,2,3);
imagesc(repmat(copt,[1 1 3])); axis image; axis off;
title(titles{2});
subplot(6,2,5);
imagesc(repmat(popt,[1 1 3])); axis image; axis off;
title(titles{3});
subplot(6,2,7);
imagesc(repmat(hopt,[1 1 3])); axis image; axis off;
title(titles{4});
subplot(6,2,9);
imagesc(repmat(nopt,[1 1 3])); axis image; axis off;
title(titles{5});
subplot(6,2,11);
imagesc(repmat(eopt,[1 1 3])); axis image; axis off;
title(titles{6});

fullpage portrait

drawnow

% 
% COUNT OPTIMAL NUMBER OF GABORS IN KERNEL
% 

% remove flipped dc component
if length(params.stimfilterparms)>=4 & params.stimfilterparms{4}>0,
   H(find(cfilt==(Xmax-1)*Xmax+1),:)=0;
end

if ~strcmp(params.kernfmt,'space'),
   tsf=pfft2sf(H,params.kernfmt);
else
   tsf=reshape(H(:,:,1),sqrt(size(H,1)),sqrt(size(H,1)),size(H,2));
end
etsf=std(tsf,0,3).*sqrt(params.bootcount);
tsf=mean(tsf,3);
%tsf=shrinkage(tsf,etsf,0.5);
%tsf=tsf(:,:,1:10);
kcount=size(tsf,3);

Nmax=6;
fitmtx=zeros(spacecount,3,Nmax);
gabsim=zeros(blen,Nmax).*nan;
allbeta={};

if strcmp(params.kernfmt,'space'),
   Nmax=6;
   fitmtx=zeros(spacecount,5,Nmax);
   ccraw=ones(kcount,Nmax+2);
   gabsim=zeros(blen,Nmax).*nan;
   cc=repmat(mean(ccraw,1),[2 1]);
   if ~USERESFILE,
      cc(1)=predxc(1);
   else
      cc(1)=predxc(end,params.nlidxsave);
   end
elseif ~strcmp(params.kernfmt,'space'),
   Nmax=6;
   fitmtx=zeros(spacecount,5,Nmax);
   ccraw=ones(kcount,Nmax+2);
   gabsim=zeros(blen,Nmax).*nan;
   
   fprintf('bootidx=');
   for bootidx=1:kcount,
      fprintf(' %d',bootidx);
      
      if ~USERESFILE,
         strfsim=(H(:,bootidx)'*(bstim'-repmat(mS(:,bootidx),1,blen)))' + ...
                 repmat(nlparms(1,bootidx,1),[blen 1]);
      else
         strfsim=[];
      end
      
      for N=1:Nmax,
         fprintf('.');
         if N<3,
            % pos only
            ttsf=tsf(:,:,bootidx).*(tsf(:,:,bootidx)>0);
            [beta,beta0]=fitgaussfp(ttsf,mod(N-1,2)+1);
            allbeta{bootidx,N}=beta;
         elseif N<5
            % neg only
            ttsf=tsf(:,:,bootidx).*(tsf(:,:,bootidx)<0);
            [beta,beta0]=fitgaussfp(ttsf,mod(N-1,2)+1);
            allbeta{bootidx,N}=beta;
         else
            ttsf=tsf(:,:,bootidx);
            beta=[allbeta{bootidx,N-4}; allbeta{bootidx,N-2}];
            allbeta{bootidx,N}=beta;
         end
         
         tfit=reshape(gaussfpN(beta,Xmax),Xmax,Xmax);
         tH=tfit(cfilt);
         if strcmp(params.kernfmt,'pfft+4'),
            tH=repmat(tH,[4 1]);
         end
         
         if ~USERESFILE,
            gabsim(:,N)=(tH'*(bstim(:,:)'-repmat(mS(:,1),1,blen)))' + ...
                repmat(nlparms(1,bootidx,1),[blen 1]);
            ccraw(bootidx,N+2)=xcov(strfsim,gabsim(:,N),0,'coeff');
         end
         
         if bootidx==1,
            fitmtx(:,1,N)=repmat(ttsf(cfilt),[phasecount,1]);
            fitmtx(:,2,N)=repmat(tfit(cfilt),[phasecount,1]);
            fitmtx(:,3,N)=fitmtx(:,1,N)-fitmtx(:,2,N);
         end
      end
   end
   
   cc=repmat(mean(ccraw,1),[2 1]);
   if ~USERESFILE,
      cc(1)=predxc(1);
   else
      cc(1)=predxc(end,params.nlidxsave);
   end
   
elseif 0,
   ccraw=ones(kcount,Nmax+2);
   fprintf('bootidx=');
   for bootidx=1:kcount,
      fprintf(' %d',bootidx);
      
      strfsim=(H(:,bootidx)'*(bstim'-repmat(mS(:,bootidx),1,blen)))' + ...
              repmat(nlparms(1,bootidx,1),[blen 1]);
      
      for N=1:Nmax,
         fprintf('.');
         [beta,beta0]=fitgaussfp(tsf(:,:,bootidx,1),N);
         allbeta{bootidx,N}=beta;
         
         tfit=reshape(gaussfpN(beta),Xmax,Xmax);
         
         tH=tfit(cfilt);
         gabsim(:,N)=(tH'*(bstim(:,:)'-repmat(mS(:,1),1,blen)))' + ...
             repmat(nlparms(1,bootidx,1),[blen 1]);
         
         if bootidx==1,
            ttsf=tsf(:,:,bootidx);
            fitmtx(:,1,N)=ttsf(cfilt);
            fitmtx(:,2,N)=tfit(cfilt);
            fitmtx(:,3,N)=fitmtx(:,1,N)-fitmtx(:,2,N);
         end
         
         ccraw(bootidx,N+2)=xcov(strfsim,gabsim(:,N),0,'coeff');
      end
   end
   cc=repmat(mean(ccraw,1),[2 1]);
   cc(1)=predxc(1);
   
elseif 0,
   strfsim=ones(size(bresp,1),1).*nan;
   fprintf('bootidx=');
   for bootidx=1:kcount,
      fprintf(' %d',bootidx);
      
      if kcount==1,
         ppidx=vidx;
      else
         ppidx=vidx((round((bootidx-1)/kcount*vcount)+1):...
                    round(bootidx/kcount*vcount));
      end
      
      strfsim(ppidx)=(H(:,bootidx)'* ...
          (bstim(ppidx,:)'-repmat(mS(:,bootidx),1,length(ppidx))))' + ...
          repmat(nlparms(1,bootidx,1),[length(ppidx) 1]);
      
      for N=1:Nmax,
         fprintf('.');
         [beta,beta0]=fitgaussfp(tsf(:,:,bootidx,1),N);
         tfit=reshape(gaussfpN(beta),Xmax,Xmax);
         
         tH=tfit(cfilt);
         gabsim(ppidx,N)=(tH'*(bstim(ppidx,:)'-...
                               repmat(mS(:,1),1,length(ppidx))))' + ...
             repmat(nlparms(1,bootidx,1),[length(ppidx) 1]);
         
         if bootidx==1,
            ttsf=tsf(:,:,bootidx);
            fitmtx(:,1,N)=ttsf(cfilt);
            fitmtx(:,2,N)=tfit(cfilt);
            fitmtx(:,3,N)=fitmtx(:,1,N)-fitmtx(:,2,N);
         end
      end
   end
   
   ggidx=find(~isnan(bresp(:,1)+strfsim));
   cc=[xcov(bresp(ggidx),rprec(ggidx,1),0,'coeff') ...
       xcov(bresp(ggidx),strfsim(ggidx),0,'coeff') ...
       xcov(bresp(ggidx),gabsim(ggidx,1),0,'coeff') ...
       xcov(bresp(ggidx),gabsim(ggidx,2),0,'coeff') ...
       xcov(bresp(ggidx),gabsim(ggidx,3),0,'coeff') ...
       xcov(bresp(ggidx),gabsim(ggidx,4),0,'coeff') ...
       xcov(bresp(ggidx),gabsim(ggidx,5),0,'coeff') ...
       xcov(bresp(ggidx),gabsim(ggidx,6),0,'coeff') ...
      ];
   cc=[cc; cc./cc(2)]
   ccraw=cc;
end
fprintf('\n');


if ~strcmp(showextra,'none');
   figure(2);
   clf
end

if strcmp(showextra,'gabor');
   titles={};
   for N=1:Nmax,
      titles{N}=sprintf('%s - %d gabors (ccfrac=%.3f)',cellid,N,cc(2,N+2));
   end
   showkern(fitmtx,kernfmt,[Xmax Xmax],titles);
   
   beta=[allbeta{1,5};allbeta{1,6}];
   rowcount=size(fitmtx,3);
   colcount=size(fitmtx,2);
   for N=1:Nmax,
      [fe,fo]=gaborcurve(beta((N-1)*5+(1:5)));
      fe=fe-fo(1);
      if max(abs(fe(:)))>0,
         fe=fe./max(abs(fe(:)))./2 + 0.5;
      else
         fe=fe+0.5;
      end
      fo=fo-fo(1);
      if max(abs(fo(:)))>0,
         fo=fo./max(abs(fo(:)))./2 + 0.5;
      else
         fo=fo+0.5;
      end
      
      subplot(rowcount,colcount,(N-1)*colcount+4);
      imagesc(repmat(fe,[1 1 3]));
      axis image; axis off;
      subplot(rowcount,colcount,(N-1)*colcount+5);
      imagesc(repmat(fo,[1 1 3]));
      axis image; axis off;
   end
   fullpage portrait
   
elseif strcmp(showextra,'curves'),
   % show marginal tuning curves
   disp('computing tuning curves');
   obincount=16;
   sfbincount=8;
   ortuning=zeros(obincount,params.bootcount);
   sftuning=zeros(sfbincount,params.bootcount);
   orslice=zeros(obincount,params.bootcount);
   sfslice=zeros(sfbincount,params.bootcount);
   seprat=zeros(params.bootcount,1);
   cm=zeros(params.bootcount,1);
   cmstd=zeros(params.bootcount,1);
   
   for bootidx=1:params.bootcount
      tsfIR=pfft2sf(H(:,bootidx,1),params.kernfmt);
      [tsfgr,obins,sfbins]=sf2gr(tsfIR,obincount,sfbincount);
      psfgr=tsfgr.*(tsfgr>0);
      
      if bootidx==1,
         tsf0=pfft2sf(mean(H(:,:,1),2));
         tsf0=sf2gr(tsf0,obincount,sfbincount);
         mor=mean(tsf0,2);
      end
      
      [u,s,v]=svd(tsfgr);
      if sum(diag(s))>0,
         seprat(bootidx)=s(1)./sum(diag(s));
      end
      s(1)=1;
      if 1,
         if u(:,1)'*mor >0,
            por=u(:,1).*s(1);
            psf=v(:,1).*s(1);
         else
            por=-u(:,1).*s(1);
            psf=-v(:,1).*s(1);
         end
         mor=mor+por;
      else
         por=mean(tsfgr,2);
         psf=mean(tsfgr',2);
      end
      
      % adjust circstats to deal with extra pi half of circle
      [cm(bootidx),cstd(bootidx)]=circstats(por(:,1));
      cstd(bootidx)=cstd(bootidx)./2 .*180/pi; 
      cm(bootidx)=cm(bootidx); % ./2 .*180/pi;
      
      ortuning(:,bootidx)=por(:,1);
      sftuning(:,bootidx)=psf(:,1);
      
      peaksfidx=min(find(psf(:,1)==max(psf(:,1))));
      orslice(:,bootidx)=tsfgr(:,peaksfidx);
      peakoridx=min(find(por(:,1)==max(por(:,1))));
      sfslice(:,bootidx)=tsfgr(peakoridx,:)';
   end
   [mcm,ecm]=circstats(ones(size(cm)),cm);
   mcm=mod(mcm.*180/pi/2,180);
   ecm=ecm.*180/pi/2.*sqrt(params.bootcount-1);
   mcstd=mean(cstd);
   ecstd=std(cstd,1).*sqrt(params.bootcount-1);
   mseprat=mean(seprat);
   eseprat=std(seprat,1).*sqrt(params.bootcount-1);
   
   mor=mean(ortuning,2);
   
   if 1,
      tor=ortuning-repmat(mean(ortuning,1),[size(ortuning,1) 1]);
      dor=diff([tor; tor(1,:)]);
      dor2=diff([tor(end,:); tor]);
      eor=(std(dor,1,2)+std(dor2,1,2))./2.*sqrt(params.bootcount-1);
   else
      eor=std(ortuning-repmat(mean(ortuning,1),[size(ortuning,1) 1]),1,2).*sqrt(params.bootcount-1);
      %eor=std(ortuning,1,2).*sqrt(params.bootcount-1);
   end
   
   msf=mean(sftuning,2);
   esf=std(sftuning-repmat(mean(sftuning,1),[size(sftuning,1) 1]),1,2).*sqrt(params.bootcount-1);
   
   peaksfidx=find(msf==max(msf));
   morslice=tsfgr(:,peaksfidx);
   
   tsf=pfft2sf(mean(H(:,:,1),2),kernfmt);
   tsf=flipud(sf2gr(tsf,obincount,sfbincount)');
   tsfpred=flipud((mor*msf')');
   ttsf=cat(3,tsf(:),tsf(:),tsfpred(:),tsfpred(:));
   
   showkern(ttsf,'space',[length(sfbins) length(obins)],{},1);
   colormap(redblue);
   
   subplot(4,1,1);
   title(sprintf('%s predcc=%.2f',cellid,predxc(1)));
   
   subplot(4,1,4);
   title(sprintf('sepidx=%.2f',mseprat));
   
   subplot(4,1,2);
   errorbar(sfbins,msf,esf,'k-','LineWidth',2);
   %plot(sfbins,msf,'k-','LineWidth',2);
   %plot(sfbins,mean(sfslice,2),'k-','LineWidth',2);
   hold on
   %plot([sfbins(1) sfbins(1)],max(msf)+mean(esf).*[1 -1],'k-','LineWidth',2);
   sfedge=(sfbins(2)-sfbins(1))./2;
   plot([sfbins(1)-sfedge sfbins(end)+sfedge],[0 0],'k--');
   hold off
   set(gca,'YTickLabel',[]);
   axis([sfbins(1)-sfedge sfbins(end)+sfedge ...
         min([msf-mean(esf); -mean(esf)])*1.1 ...
         max([msf+mean(esf); mean(esf)]).*1.1]);
   xlabel('spatial freq (cyc/RF)');
   ylabel('rel response');
   
   subplot(4,1,3);
   errorbar(obins,mor,eor,'k-','LineWidth',2);
   %plot(obins,mor,'k-','LineWidth',2);
   %plot(obins,mean(orslice,2),'k-','LineWidth',2);
   hold on
   %plot([obins(1) obins(1)],max(mor)+mean(eor).*[1 -1],'k-','LineWidth',2);
   oedge=(obins(2)-obins(1))./2;
   plot([obins(1)-oedge obins(end)+oedge],[0 0],'k--','LineWidth',1);
   hold off
   set(gca,'YTickLabel',[]);
   axis([obins(1)-oedge obins(end)+oedge ...
         min([mor-mean(eor); -mean(eor)])*1.1 ...
         max([mor+mean(eor); mean(eor)]).*1.1]);
   title(sprintf('cm=%.1f sem=%.1f cstd=%.1f',mcm,ecm,mcstd));
   xlabel('orientation (deg)');
   ylabel('rel response');
   
   set(gcf,'PaperOrientation','portrait','PaperPosition',[3 2.5 2.5 6]);
  
   %keyboard
   
   % display alternative "best" natural stimuli
   figure(3);
   ttitles={sprintf('%s best pred nat stim',cellid)};
   gam=1.0;
   nshow=reshape(noptset,64*64,30).^(1./gam);
   offset=0.55;
   %offset=repmat(mean(nshow,1),64*64,1);
   showstim(nshow-offset,...
            'space',[64 64],3,10,ttitles);
   %showkern(reshape(noptset,64*64,10,3)-127,'space',[64 64],ttitles);
   fullpage portrait
   
   figure(4);
   ttitles={sprintf('%s best pred nat stim pfftgr',cellid)};
   showkern(reshape(fnset,size(fnset,1),10,3),kernfmt,[64 64],ttitles);
   fullpage portrait
   
elseif strcmp(showextra,'grat'),
   % sample cart and ncart gratings, contrast scaled to match
   % strength of response. thus the higher contrast stimuli are
   % "better"
   
   cartbest=cartmatch(:,:,cm3);
   [xx,cartbest]=sort(-cartbest(:));
   cartbest=cartbest([1:5 end-4:end]);
   cpstim=cartpfft(:,:,:,cm3);
   cpstim=cpstim(:,cartbest(1:10));
   ncbest=cat(3,polmatch(:,:,:),hypmatch(:,:,:));
   [xx,ncbest]=sort(-ncbest(:));
   ncbest=ncbest([1:3:13 end-12:3:end]);
   ncstim=cat(2,polpfft(:,:),hyppfft(:,:));
   ncstim=ncstim(:,ncbest(1:10));
   polbest=polmatch(:,:,:);
   [xx,polbest]=sort(-polbest(:));
   ppstim=polpfft(:,polbest(1:10));
   hypbest=hypmatch(:,:,:);
   [xx,hypbest]=sort(-hypbest(:));
   hpstim=hyppfft(:,hypbest(1:10));
   [xx,nmset]=sort(-natmatch(1:length(natmatch)/sc));
   nmset=nmset([1:2:30 end-29:2:end]);
   npstim=reshape(nrmov(:,nmset),spacecount,10,3);
   
   spix=Xmax*3;
   sampstim=zeros(spix,spix,10,8);
   for ii=1:10,
      [x1,x2]=ind2sub([dcount dcount2],cartbest(ii));
      [t,copt]=ncart_pfft([cartfrange(x1),cartorange(x2),cartprange(cm3)],...
           spix,'cart',params.stimfiltercmd,params.stimfilterparms,[mm ss]);
      [x1,x2,x3]=ind2sub([dcount dcount2 dcount3*2],ncbest(ii));
      if x3<=dcount3,
         [t,ncopt]=ncart_pfft([polprange(x3),polcrange(x1)./3,polrrange(x2)],...
           spix,'pol',params.stimfiltercmd,params.stimfilterparms,[mm ss]);
      else
         [t,ncopt]=ncart_pfft([hyporange(x1),hypprange(x3-dcount3),hypfrange(x2)./3],...
              spix,'hyp',params.stimfiltercmd,params.stimfilterparms,[mm ss]);
      end
      
      [x1,x2,x3]=ind2sub([dcount dcount2 dcount3],polbest(ii));
      [t,popt]=ncart_pfft([polprange(x3),polcrange(x1)./3,polrrange(x2)],...
           spix,'pol',params.stimfiltercmd,params.stimfilterparms,[mm ss]);
      %[x1,x2,x3]=ind2sub([dcount dcount2 dcount3],hypbest(ii));
      %[t,hopt]=ncart_pfft([hyporange(x1),hypprange(x3),hypfrange(x2)./3],...
      %     spix,'hyp',params.stimfiltercmd,params.stimfilterparms,[mm ss]);
      
      sampstim(:,:,ii,3)=copt;
      sampstim(:,:,ii,4)=ncopt;
      sampstim(:,:,ii,5)=popt;
   end
   
   nshow=loadimframes(natrevfile,nok(nmset),96,0,spix);

   nshow=nshow./255;
   gamma=3;
   nshow=nshow.^(1/gamma);
   nshow=nshow.*255;
   
   nmm=repmat(mean(reshape(nshow,spix*spix,1,30)),spix,spix);
   nss=repmat(std(reshape(nshow,spix*spix,1,30)),spix,spix);
   nss(find(nss==0))=1;
   nshow=((nshow-nmm)./nss .* ss+mm);
   sampstim(:,:,:,6:8)=reshape(nshow./1.4,spix,spix,10,3);
   
   sampstim(find(sampstim>mm*2.5))=mm*2.2; sampstim(find(sampstim<0))=0;
   %showkern(reshape(sampstim,spix*spix,10,size(sampstim,4))-mm.*1.1,'space');
   showstim(reshape(sampstim,spix*spix,10*size(sampstim,4))-mm.*1.1);
   
   subplot(4,1,1);
   axis off
   
   subplot(4,2,1);
   
   cartm=flipud(sort(cartmatch(:)));
   polm=flipud(sort(polmatch(:)));
   hypm=flipud(sort(hypmatch(:)));
   ncm=flipud(sort([hypmatch(:);polmatch(:)]));
   ncm=ncm(round(linspace(1,length(ncm),length(polmatch(:)))));
   natm=flipud(sort(natmatch(:)));
   natm=natm(round(linspace(10,length(natm)-10,length(polmatch(:)))));
   showrange1=1:length(cartm);
   showrange2=length(cartm)-50:length(cartm);
   
   rankpct=showrange1./length(cartmatch(:)).*100;
   plot(rankpct,cartm(showrange1),'g','LineWidth',1);
   hold on
   %plot(polm(showrange1),'r','LineWidth',1);
   %plot(hypm(showrange1),'b','LineWidth',1);
   plot(rankpct,ncm(showrange1),'r','LineWidth',1);
   plot(rankpct,natm(showrange1),'b','LineWidth',1);
   %plot(55+(1:length(showrange2)),cartm(showrange2),'g--');
   %plot(55+(1:length(showrange2)),polm(showrange2),'r');
   %plot(55+(1:length(showrange2)),hypm(showrange2));
   %plot(55+(1:length(showrange2)),natm(showrange2),'g');
   hold off
   %legend('cart','pol','hyp','nat',-1);
   legend('cart','ncart','nat',-1);
   a=axis;
   axis([0 rankpct(end) a(3:4)])
   %axis([0 rankpct(end) a(3)-(a(4)-a(3)).*0.2 a(4)+(a(4)-a(3)).*0.2])
   
   set(gcf,'PaperOrientation','portrait','PaperPosition',[0.5 2 7.5 6]);
   
   %keyboard
   
   figure(3);
   showkern(cat(3,cpstim.*nan,cpstim.*nan,cpstim,ncstim,ppstim,npstim),'pfftgr')
   set(gcf,'PaperOrientation','portrait','PaperPosition',[0.5 2 7.5 6]);
   
elseif strcmp(showextra,'grat2'),
   % sample cart and ncart gratings, contrast scaled to match
   % strength of response. thus the higher contrast stimuli are
   % "better"
   tnp1=permute(squeeze(cartmatch(:,:,cm3)),[4 1 3 2]);
   tnp2=permute(squeeze(polmatch(:,:,pm3)),[4 1 3 2]);
   tnp3=permute(squeeze(hypmatch(:,:,hm3)),[4 1 3 2]);
   mm=min([tnp1(:);tnp2(:);tnp3(:)]);
   
   tnc1=permute(squeeze(cart(:,:,:,cm3)),[1 2 4 3]);
   tnc1=tnc1.*repmat(tnp1-mm,[size(tnc1,1) 1 1 1]);
   
   tnc2=permute(squeeze(pol(:,:,:,pm3)),[1 2 4 3]);
   tnc2=tnc2.*repmat(tnp2-mm,[size(tnc2,1) 1 1 1]);
   
   tnc3=permute(squeeze(hyp(:,:,:,hm3)),[1 2 4 3]);
   tnc3=tnc3.*repmat(tnp3-mm,[size(tnc3,1) 1 1 1]);
   
   tnc=cat(2,tnc1,tnc2,tnc3);
   showstim(tnc(:,:),'space',[1 1].*sqrt(size(tnc,1)),dcount2,dcount*3);
   
   set(gcf,'PaperOrientation','portrait','PaperPosition',[0.75 2 7 3.5]);
   
elseif strcmp(showextra,'ed');
   
   % display best ed stim
   estim=edmov(:,:,edset);
   estim=reshape(estim,Xmax*Xmax,10,length(edset)/10);
   
   % [edocount edsepcount edpolcount*sc]
   
   estim=edmov(:,:,em1+(0:(edsepcount-1))*edocount+edocount*edsepcount*(em3-1));
   em=edmatch(em1+(0:(edsepcount-1))*edocount+edocount*edsepcount*(em3-1));
   em=em-min(em(:));
   %em=edmatch(1:edocount*edsepcount)-min(edmatch(1:edocount*edsepcount));
   estim=estim.*repmat(permute(em,[2 3 1]),[Xmax Xmax]);
   %estim=reshape(estim,Xmax*Xmax,edocount,edsepcount);
   %estim=permute(estim,[1 3 2]);
   estim=reshape(estim,Xmax*Xmax,size(estim,3));
   
   showkern(estim,'space')
   %showstim(estim(:,:),'space',[Xmax Xmax],edocount,edsepcount);
   keyboard
   
elseif strcmp(showextra,'cart');
   % sample carts
   tnp=permute(squeeze(cartmatch(:,:,cm3)),[4 1 3 2]);
   tnp=tnp-min(tnp(:));
   tnc=permute(squeeze(cart(:,:,:,cm3)),[1 2 4 3]);
   tnc=tnc.*repmat(tnp,[size(tnc,1) 1 1 1]);
   tnc=tnc(:,:);
   
   showstim(tnc,'space',[1 1].*sqrt(size(tnc,1)),dcount2,dcount);
   
elseif strcmp(showextra,'cartpfft');
   tnp=permute(squeeze(cartmatch(:,:,cm3)),[4 1 3 2]);
   tnp=tnp-min(tnp(:));
   tnc=permute(squeeze(cartpfft(:,:,:,cm3)),[1 2 4 3]);
   tnc=(tnc-repmat(mS(:,1),size(tnp))).*repmat(tnp,[size(tnc,1) 1 1 1]);
   tnc=tnc(:,:);
   
   showstim(tnc,params.kernfmt,[Xmax Xmax],dcount,dcount2);
elseif strcmp(showextra,'pol');
   % sample pols
   tnp=permute(squeeze(polmatch(:,:,pm3)),[4 1 3 2]);
   tnp=tnp-min(tnp(:));
   tnc=permute(squeeze(pol(:,:,:,pm3)),[1 2 4 3]);
   tnc=tnc.*repmat(tnp,[size(tnc,1) 1 1 1]);
   tnc=tnc(:,:);
   
   showstim(tnc,'space',[1 1].*sqrt(size(tnc,1)),dcount2,dcount);
   
elseif strcmp(showextra,'polpfft');
   tnp=permute(squeeze(polmatch(:,:,pm3)),[4 1 3 2]);
   tnp=tnp-min(tnp(:));
   tnc=permute(squeeze(polpfft(:,:,:,pm3)),[1 2 4 3]);
   tnc=(tnc-repmat(mS(:,1),size(tnp))).*repmat(tnp,[size(tnc,1) 1 1 1]);
   tnc=tnc(:,:);
   
   showstim(tnc,params.kernfmt,[Xmax Xmax],dcount,dcount2);
elseif strcmp(showextra,'hyp');
   % sample hyps
   tnp=permute(squeeze(hypmatch(:,:,hm3)),[4 1 3 2]);
   tnp=tnp-min(tnp(:));
   tnc=permute(squeeze(hyp(:,:,:,hm3)),[1 2 4 3]);
   tnc=tnc.*repmat(tnp,[size(tnc,1) 1 1 1]);
   tnc=tnc(:,:);
   
   showstim(tnc,'space',[1 1].*sqrt(size(tnc,1)),dcount,dcount2);
   
elseif strcmp(showextra,'hyppfft');
   tnp=permute(squeeze(hypmatch(:,:,hm3)),[4 1 3 2]);
   tnp=tnp-min(tnp(:));
   tnc=permute(squeeze(hyppfft(:,:,:,hm3)),[1 2 4 3]);
   tnc=(tnc-repmat(mS(:,1),size(tnp))).*repmat(tnp,[size(tnc,1) 1 1 1]);
   tnc=tnc(:,:);
   
   showstim(tnc,params.kernfmt,[Xmax Xmax],dcount,dcount2);
end

% determine how many pixels wide the RC window is:
if params.stimloadparms{1}==0,
   stimdiamrat=params.stimwindowcrf./params.stimcrfs(1);
else
   stimdiamrat=params.stimloadparms{3}./params.stimloadparms{1};
end

if 1,
   % kludge b/c v1 cells aren't in nsl database
   windowpix=stimdiamrat.*120;
   windowdeg=0.8788;
   fprintf('stim window used (deg) = %.2f\n',windowdeg);
   
elseif 1
   sql=['SELECT * FROM sCellFile',...
        ' WHERE cellid=''',cellid,'''',...
        ' AND concat(stimpath,stimfile)=''',params.stimfiles{1},''''];
   
   cellfiledata=mysql(sql);
   
   windowpix=stimdiamrat.*cellfiledata(1).stimwindowsize;
   windowdeg=rfdiam(cellid,windowpix);
   fprintf('stim window used (deg) = %.2f\n',windowdeg);
else
   windowpix=stimdiamrat.*128;
   windowdeg=9.0026;
   fprintf('GUESSING!!!!  stim window used (deg) = %.2f\n',windowdeg);
end


drawnow

% if nargout==0, don't return anything in res
clear res

if nargout>0,
   
   % if requested output useful results to res 
   res.cellid=cellid;
   res.bad=0;
   res.windowdeg=windowdeg;
   res.cc=cc;
   res.H=H;
   res.ccraw=ccraw;
   res.meanactresp=meanactresp;
   res.meanpredresp=mean(natmatch(:));
   res.cartmin=cartmatch(cmin);
   res.cartmax=[cartmax max(cartmatchcc(:))] ;
   res.cartparms=[cartfrange(cm1),cartorange(cm2),cartprange(cm3)];
   res.polmin=polmatch(pmin);
   res.polmax=[polmatch(pm1,pm2,pm3,1) max(polmatchcc(:))];
   res.polparms=[polprange(pm3),polcrange(pm1),polrrange(pm2)];
   res.hypmin=hypmatch(hmin);
   res.hypmax=[hypmatch(hm1,hm2,hm3,1) max(hypmatchcc(:))];
   res.hypparms=[hyporange(hm1),hypprange(hm3),hypfrange(hm2)];
   res.natparms=nmax;
   res.edparms=[em1 em2 em3];
   res.edccparms=[ecc1 ecc2 ecc3];
   res.edmatch=edmatch;
   res.pnrat=pnrat;
   res.natfracbetter=natfracbetter;
   res.natovercart=natovercart;
   res.natmin=[natmin1 natmin2 nsort(1)];
   res.natmax=[natmax1 natmax2 nsort(end)];
   res.natmaxerr=[natmax1err natmax2err natmaxerr];
   
   res.cartvs=[res.cartmax(1) ncartmax natmax1 edmax];
   res.cartvserr=[cartmaxerr ncartmaxerr natmax1err edmaxerr];
   res.cartvscc=[cartccmax ncartccmax natccmax edccmax];
   res.cartvsccmin=[cartccmin ncartccmin natccmin edccmin];
   res.cartvsccerr=[cartccmaxerr ncartccmaxerr natccmaxerr edccmaxerr];
   res.stimclassstd=[cartccmax ncartccmax natccmax edccmax];
   res.synthvs=[synthmax natmax2];
   res.synthvserr=[synthmaxerr natmax2err];
   
   res.sparseness=[sparseness(cartmatch(:)-mmmin) ...
                   sparseness(polmatch(:)-mmmin) ...
                   sparseness(hypmatch(:)-mmmin) ...
                   sparseness(natmatch(:)-mmmin)];
   res.sparseplus=[sparseness(cartmatch(:)-mmmin) ...
                   sparseness([polmatch(:);hypmatch(:)]-mmmin) ...
                   sparseness([cartmatch(:);polmatch(:);hypmatch(:)]-mmmin) ...
                   sparseness(natmatch(:)-mmmin)];
   
   res.entropy=[hcart hpol hhyp hnat];
   
   res.randxc=randxc;
   res.beta1=cat(2,allbeta{:,end-1});
   res.beta2=cat(2,allbeta{:,end});
end

%keyboard

return

