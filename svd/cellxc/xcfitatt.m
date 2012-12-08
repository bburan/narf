% xcfitatt.m
%
% we have two things to fit: 
% 1. the svd regularization parameter (sfsidx)
% 2. the shrinkage filter fudge factor (sfiltsigma)
%
% specs (so this is interchangable with other fit routines):
% inputs: mH, eH (mean and std error of kernels)
%         params (bunch of crap listed in cellxcnodb.m)
%         stim,resp (if fitfrac=0, exploratory fit data, used for
%                   fitting other parms as well. otherwise, data is
%                   loaded using params.times(2))
%
% strf structure:
% strf( ).h    :  linear kernel
% strf( ).mS   :  mean stimulus used in fit for subtracting before pred
% strf( ).name :  info about what this is
% strf( ).parms.sfsidx   : best fit idx
% strf( ).parms.sigidx   : best fit idx
% strf( ).parms.kernfmt  : kernfmt for showkern
% strf( ).parms.iconside : stim size info for showkern
% strf( ).nlparms : parms for output nl
% strf( ).nltype  : string desc of nl
% strf( ).xc      : results of fit xc test.
%
% strf array can be many dimensional (eg, nl X att)
%

disp('xcfitatt.m:');

params.matchattsfs=getparm(params,'matchattsfs',0);

fparams=params;
fparams.nrandxcov=0;

%
% figure out dimensions of relevant matrices
%

spacecount=size(mH,1);
if diff(fparams.maxlag)>0,
   tbincount=min([fparams.maxlag(2)+1 14]);
else
   tbincount=size(mH,2);
end
sfscount=size(mH,3);
respcount=size(mH,4);
attcount=size(mH,5);
lagrange=(1-fparams.maxlag(1)):(tbincount-fparams.maxlag(2)+10);

sfsidxtouse=(1:fparams.sfscount)';  % test all sfs cutoffs

% define a range of sigmas for combining resampled estimates (via
% shrinkage filter)
if length(fparams.sffiltsigma)>1,
   sigrange=fparams.sffiltsigma(:)';
   
elseif fparams.respfmtcode==0,
   
   %
   % USED FOR V1 ANALYSIS
   % (Adjusted down by a factor of sqrt(resampcount) after eH
   % computation was corrected to be actual standard error.
   %
   sigrange=exp(linspace(log(0.8),log(1.6),fparams.sffiltsigma));
else
   
   % USED FOR PFTH-DRIVEN ANALYSIS
   sigrange=exp(linspace(log(0.8),log(1.5),fparams.sffiltsigma));
end

if 0, % run fast  for debugging
   sigrange=exp([0 1]);
   sfsidxtouse=0:2:10;
   sfscount=length(sfsidxtouse);
end
sigcount=length(sigrange);

if ~exist('OUTNLMODE','var'),
   OUTNLMODE=0;
end
if OUTNLMODE==0,
   nlnames={'no att','dc att','gain att','dc/g att',...
            'hng att','all att'};
   localnls=6;
elseif OUTNLMODE==1,
   nlnames={'no att','x10','x90','D','A','all att'};
   localnls=6;
elseif OUTNLMODE==2,
   nlnames={'no att','dc att','gain att','thr att',...
            'hng att','all att'};
   localnls=6;
elseif OUTNLMODE==3,
   nlnames={'glob none','glob dc','glob dcg',...
            'dc att','gain att','all att'};
   localnls=[4 5 6];
end
fprintf('OUTNLMODE=%d\n',OUTNLMODE);
nlcount=length(nlnames);

% for PFTH preds, separate score for each latency bin
if fparams.respfmtcode==0,
   latcount=respcount;
else
   latcount=tbincount;
end
% based on
%predxc=ones(strfcount,attcount,latcount) * nan;
xc0=zeros(sfscount,sigcount,nlcount);
xc=zeros(sfscount,sigcount,nlcount);
xcatt=zeros(sfscount,sigcount,attcount-1,nlcount);
nmse=zeros(sfscount,sigcount,nlcount);
mutinfo=zeros(sfscount,sigcount,nlcount,attcount,latcount);
sfsfit=zeros(attcount-1,nlcount);
sigfit=zeros(attcount-1,nlcount);
expxc=zeros(attcount-1,nlcount);

threshparm=zeros(sfscount,sigcount,latcount,nlcount,attcount);
threshfit=zeros(latcount,nlcount,attcount);

clear tstrf

%samplecount=6; %number of evenly distributed sample kernels to save
%hfsample=zeros(spacecount,tbincount,samplecount,attcount);

attidx=1;
if exist('fdata','var'),
   
   % do nothing
   
elseif ~exist('VCELLXC','var') | VCELLXC==1 | fparams.fitfrac>0,
   % load stim and response segments to fit. reload the stimulus
   % and response for each attentional state
   
   [fdata.stim,fdata.resp]=...
       xcloadstimresp(fitfile(:,attidx),fitstartframe(:,attidx),...
                      fitstopframe(:,attidx),fparams);
elseif size(resp,1)>5000,
   rgoodidx=find(~isnan(resp(:,1)));
   rminidx=rgoodidx(1);
   if length(rgoodidx)>=5000,
      rmaxidx=rgoodidx(5000);
   else
      rmaxidx=rgoodidx(end);
   end
   fprintf('trimming exp data to frames %d-%d for fitting\n',...
              rminidx,rmaxidx);
   fdata.resp=resp(rminidx:rmaxidx,:);
   fdata.stim=stim(rminidx:rmaxidx,:);
else
   fdata.stim=stim;
   fdata.resp=resp;
end

fprintf('Testing fit (sfsidx/%d):',sfscount);
respsave=fdata.resp;
nlparmset={};
nlnameset={};

for sfsidx=1:sfscount,
   fprintf('(%d)',sfsidx);
   
   % does this do anything?
   for respidx=1:respcount,
      
      % apply variable shrinkage 
      for sigidx=1:sigcount,
         
         tstrf1=[];
         for attidx=1:attcount,
            % current kernel is shrinkage-filtered from current (sfs,r,attidx)
            % only use causal time lags for pred!
            thf=shrinkage(...
               mH(:,(1:tbincount)-fparams.maxlag(1),sfsidx,respidx,attidx),...
               eH(:,(1:tbincount)-fparams.maxlag(1),sfsidx,respidx,attidx),...
               sigrange(sigidx));
            
            tstrf1(1,attidx).h=thf;
            tstrf1(1,attidx).nlparms=0;
            tstrf1(1,attidx).nltype='none';
            tstrf1(1,attidx).mS=mSall(:,respidx,attidx);
            tstrf1(1,attidx).parms.kernfmt=fparams.kernfmt;
            tstrf1(1,attidx).parms.iconside=fparams.iconside;
            tstrf1(1,attidx).name=...
                sprintf('%s sfs=%d sig=%d attidx=%d',fparams.cellid,...
                        sfsidx,sigidx,attidx);
         end
         
         fdata.resp=respsave(:,1);
         
         % generate straight linear preds from each kernel. predict
         % entire estimation data set using shrinkaged
         % kernels. this method seems to work better than choosing
         % the cutoff with excluded parts from each jackknife.
         tfitres=xcval(tstrf1(:),fparams,fdata);
         linpred=reshape(tfitres.mod_psth{1},length(fdata.resp),attcount);
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % fit att-dependent output NLs
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         tstrf=cat(1,...
                   repmat(tstrf1(:,1),nlcount-length(localnls),attcount-1),...
                   repmat(tstrf1(1,2:attcount),length(localnls),1));
         
         % record names in appropriate kernels
         for nlidx=1:nlcount,
            for fitidx=1:attcount-1,
               tstrf(nlidx,fitidx).name=...
                   sprintf('%s sfs=%d sig=%d att=%d %s',...
                           cellid,sfsidx,sigidx,fitidx,nlnames{nlidx});
            end
         end
         
         respforfit=respsave;
         
         kvafitnl;
         for nlidx=1:nlcount,
            if std(linpred(tokidx,1))>0,
               xc0(sfsidx,sigidx,nlidx)=xcov(linpred(tokidx,1),fitresp,...
                                             0,'coeff');
            end
            
            if std(nlinpred(tokidx,nlidx))>0,
               xc(sfsidx,sigidx,nlidx)=xcov(nlinpred(tokidx,nlidx),...
                                            fitresp,0,'coeff');
            end
            
            for fitidx=1:attcount-1,
               attidx=fitidx+1;
               atokidx=find(~isnan(respsave(tokidx,attidx)));
               if std(nlinpred(tokidx(atokidx),nlidx))>0 & ...
                     std(fitresp(atokidx))>0,
                  xcatt(sfsidx,sigidx,fitidx,nlidx)=...
                      xcov(nlinpred(tokidx(atokidx),nlidx),...
                           fitresp(atokidx),0,'coeff');
               end
               
               nlparmset{sfsidx,sigidx,nlidx,fitidx}=tstrf(nlidx,fitidx).nlparms;
               nlnameset{sfsidx,sigidx,nlidx,fitidx}=tstrf(nlidx,fitidx).nltype;
            end
         end
      end
   end
   if not(isempty(BATQUEUEID)),
      % record latest (sfsidx,attidx) predicted
      dbsetqueue(BATQUEUEID,1);
   end
end

fprintf('\n');
fdata.resp=respsave;

clear strf
% generate optimal kernels and save in hf

% set up to smooth txc
if size(xc,1)>=8 & size(xc,2)>=3
   [XX,YY]=meshgrid(-2:2,-2:2);
   sx=0.5; sy=1.0;
   % fvvs
   % sx=1.5; sy=2.0;
   g=exp(-XX.^2./(2*sx.^2) - YY.^2./(2*sy.^2));
   g=g./sum(g(:));
else
   g=1;
end

for nlidx=1:nlcount,
   
   txc=xc(:,:,nlidx);
   txc=rconv2(txc,g);
   maxidx=min(find(txc==max(txc(:))));
   [amaxsfsidx,amaxsigidx]=ind2sub([sfscount,sigcount],maxidx);
   
   if ismember(nlidx,localnls) & amaxsfsidx>sfsfit(1,min(localnls)-1),
      amaxsfsidx=sfsfit(1,min(localnls)-1);
      amaxsigidx=sigfit(1,min(localnls)-1);
   end
   
   % go through each att state separately so that you can have
   % different cutoffs for the different local kernels. has the
   % potential of helping?
   for fitidx=1:attcount-1,
      txc=xcatt(:,:,fitidx,nlidx);
      txc=rconv2(txc,g);
      
      maxidx=min(find(txc==max(txc(:))));
      [maxsfsidx,maxsigidx]=ind2sub([sfscount,sigcount],maxidx);
      
      % this is over-ridden below to force cutoff to be the same
      %for each att state
      sfsfit(fitidx,nlidx)=maxsfsidx;
      sigfit(fitidx,nlidx)=maxsigidx;
      
      expxc(fitidx,nlidx)=xcatt(maxsfsidx,maxsigidx,fitidx,nlidx);
   end
   
   if 1,
      % choose sfsidx,sigidx based on xc for whole prediction
      % (including all att states) ie, not a separate cutoff for
      % each local att kernel
      sfsfit(:,nlidx)=amaxsfsidx;
      sigfit(:,nlidx)=amaxsigidx;
   else
      % choose min sfs for expression (conservative/smoothed)
      mm=min(find(sfsfit(:,nlidx)==min(sfsfit(:,nlidx))));
      sfsfit(:,nlidx)=sfsfit(mm,nlidx);
      sigfit(:,nlidx)=sigfit(mm,nlidx);
   end
   
   for fitidx=1:attcount-1,
      maxsfsidx=sfsfit(fitidx,nlidx);
      maxsigidx=sigfit(fitidx,nlidx);
      
      fprintf('nlidx=%d fitidx=%d (sfs,sig)=(%d,%d)\n',...
              nlidx,fitidx,maxsfsidx,maxsigidx);
      
      if ismember(nlidx,localnls),
         % local kernel for each att state
         attidx=fitidx+1;
      else
         % same kernel for all att states
         attidx=1;
      end
      
      % fill strf with best fit kernel according to sigfit and sfsfit
      th=shrinkage(mH(:,(1:tbincount)-fparams.maxlag(1),maxsfsidx,1,attidx),...
                   eH(:,(1:tbincount)-fparams.maxlag(1),maxsfsidx,1,attidx),...
                   sigrange(maxsigidx));
      
      strf(nlidx,fitidx).h=th;
      strf(nlidx,fitidx).parms.sfsfit=maxsfsidx;
      strf(nlidx,fitidx).parms.sigfit=maxsigidx;
      strf(nlidx,fitidx).parms.kernfmt=fparams.kernfmt;
      strf(nlidx,fitidx).parms.iconside=fparams.iconside;
      strf(nlidx,fitidx).parms.tbinms=fparams.tbinms;
      strf(nlidx,fitidx).parms.respfmtcode=fparams.respfmtcode;
      strf(nlidx,fitidx).mS=mSall(:,1,attidx);
      
      strf(nlidx,fitidx).nltype=nlnameset{maxsfsidx,maxsigidx,nlidx,fitidx};
      strf(nlidx,fitidx).nlparms=nlparmset{maxsfsidx,maxsigidx,nlidx,fitidx};
      strf(nlidx,fitidx).name=sprintf('%s sfs=%d sig=%d att=%d %s',...
                                      cellid,maxsfsidx,maxsigidx,...
                                      fitidx,nlnames{nlidx});
      strf(nlidx,fitidx).parms.attname=nlnames{nlidx};
   end
end

clear tresp tmod_psth tgoodidx tgoodlen tstimloadparms tt txc thf
clear r0 r1 smm smd sms scf sfsidx ridx nlidx maxidx fidx
clear d1 XX YY 

return


