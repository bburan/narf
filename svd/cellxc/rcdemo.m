% rcdemo.m
%
% sample reverse correlation routine using a model simple cell with
% dynamic natural image (DNAT) stimulus. for the model cell RC can
% be carried out in pixel domain (kernfmt='space') because the cell is
% simple.
%
% alternatively, demo can be run on actual DNAT data. this gives a
% worse result in the pixel domain than the fourier power domain
% (kernfmt='pfft').
%
% REQUIREMENTS:
% program should be run with matlab in the appropriate directory
% containing two subdirectories: data/ and code/ which contain
% appropriate support files.
%
% Basic outline of natural stim reverse correlation process:
% 1. load data. stored in one or more pairs of stimulus-response
%    files.
% 2. pre-process.  this includes applying an input nonlinearity (if
%    desired) and reshaping the data matrices to match the standard
%    format.
% 3. computer stimulus-response cross correlation and stimulus
%    autocorrelation.
% 4. remove stimulus bias for final strf
% 5. display results
%
% created SVD 8/1/02
% modified SVD 1/28/03 - added R221A (MODELDATA=0) and kernfmt options
%

%%
%% high-level parameters
%%

MODELDATA=1;      % if 1 use model cell. if 0, use real data
% toggle commenting to choose input nonlinearity
kernfmt='space';   % this kernel will be computed in space domain
%kernfmt='pfft';   % this kernel will be computed in fourier power domain

% add paths with necessary code
if exist('loadimfile') < 2,
   addpath code/mutils/
   addpath code/gen/
   addpath code/cellxc/
end

%%
%% 1. Load the stimulus and response data
%%

if MODELDATA,
   % load stimulus and response data for model simple cell shown a
   % series of dynamic natural image patches
   
   cellid='model';
   cellpath='data/';
   % load movie (and downsample it for computational tractibility)
   stim=loadimfile([cellpath,'test.natrev.andros.30_pix.imsm'],...
                  0,0,16,33,16,0,1);
   resp=respload([cellpath,'model.nr.simp.resp.mat']);
   
   % only need first column of response (averaged over trials)
   resp=resp(:,1);

else
   % load stimulus and response data for a real cell shown a
   % series of dynamic natural image patches
   % (stored in multiple files due to experimental procedures)

   cellid='R221A';
   cellpath='data/';
   stimfiles={'test.natrev.andros.32_pix.imsm',...
              'test.natrev.forestferns.32_pix.imsm'...
             };
   respfiles={'R221A.natrev.andros.32_pix.1.all.d-.psth.neg_flag',...
              'R221A.natrev.forestferns.32_pix.1.all.d-.psth.neg_flag'...
             };
   
   % load each file and append to the end of big stimulus and
   % response matrices
   stim=[];
   resp=[];
   for ii=1:length(stimfiles),
      % load, crop, and downsample the movie to 16 x 16 pixels
      % matched to 2 times the receptive field diameter
      if strcmp(kernfmt,'space'),
         tstim=loadimfile([cellpath,stimfiles{ii}],0,0,32,33,16,0,1);
      else
         % apply a 1 CRF gaussian window to each stim frame before downsample
         tstim=loadimfile([cellpath,stimfiles{ii}],0,0,32,8,16,0,1);
      end
      
      tresp=respload([cellpath,respfiles{ii}]);
      % only need first column of response (averaged over trials)
      tresp=tresp(:,1);
      
      % append to full stim and resp matrices
      stim=cat(3,stim,tstim(:,:,1:length(tresp)));
      resp=[resp;tresp];
   end
end

%%
%% 2. Pre-process the stimulus data
%%

% figure out dimensions of stimulus
movx=size(stim,1);   % size of downsampled movie in pixels
movlen=size(stim,3); % length of movie

% if we're using a non-linear input transform, apply it to the stimulus here:
if strcmp(kernfmt,'pfft'),
   % fourier power transformation, outpt is space X time
   stim=movpower(stim,0,0,0,1,0);
elseif strcmp(kernfmt,'space'),
   % space domain: reshape movie to 2 dimensions: time X space
   stim=reshape(stim,movx*movx,movlen);
else
   fprintf('invalid kernfmt=%d. aborting.\n',kernfmt);
   return
end

% transpose to get time X space
stim=stim';
spacecount=size(stim,2);

%%
%% 3. perform the stimulus/response cross correlation, also get the
%% stimulus autocorrelation
%%
maxlag=[0 10]; % range of time lags to compute STRF (14 ms bins)

% these are extracted and unwarped from movxc() and normalize0()

% invalid time bins are set to nan in resp. these need to be
% excluded from the cross correlation
rgoodidx=find(~isnan(resp));

% number of valid time bins
n=length(rgoodidx);

% calculate mean response and stimulus
mR=sum(resp(rgoodidx))' ./ n;
mS=sum(stim(rgoodidx,:))' ./ n;

% calculate stimulus-response cross correlation
fprintf('Computing stimulus-response cross correlation');
SR=zeros(spacecount,diff(maxlag)+1);
for tt=maxlag(1):maxlag(2),
   fprintf('.');
   trg=rgoodidx;
   trg=trg(find(trg-tt>0 & trg-tt<size(stim,1)));
   SR(:,tt-maxlag(1)+1)=(resp(trg)'*stim(trg-tt,:))' ./ n - mS*mR;
end
fprintf('\n');

% stimulus 2nd order spatial autocorrelation
disp('Computing spatial autocorrelation...');
sSA2=stim(rgoodidx,:)'*stim(rgoodidx,:) - mS*mS';

% first order temporal autocorrelation
disp('Computing temporal autocorrelation...');
tSA=zeros(diff(maxlag)*2+1,1);
for xx=1:spacecount,
   tstim=stim(rgoodidx,xx);
   tSA=tSA + xcorr(tstim,diff(maxlag),'biased') ./ spacecount;
end
tSA=tSA - mS'*mS./spacecount;


%%
%% 4. correct for stimulus bias (first order in time, second order in space):
%%
% select eigenvalue cutoff for pseudo-inverse. this is a nasty step
% in the process. i've preselected nice values here.

if MODELDATA,
   neigs=60;
elseif strcmp(kernfmt,'space'),
   neigs=15;
else
   neigs=30;
end
fprintf('Removing stim bias (pseudoinv includes %d eigenvectors)...\n',...
        neigs);
H=normalize(SR,sSA2,tSA,[],neigs);

%%
%% 5. predict fit data
%%
tstim=stim'-repmat(mS,[1 size(stim,1)]);
predresp=kernpredict(H,tstim,1,0);  % 1 phase, no rect
vidx=find(~isnan(resp) & ~isnan(predresp));
cc=corrcoef(predresp(vidx),resp(vidx));
Hcc=cc(1,2);

predresp=kernpredict(SR,tstim,1,0);  % 1 phase, no rect
vidx=find(~isnan(resp) & ~isnan(predresp));
cc=corrcoef(predresp(vidx),resp(vidx));
SRcc=cc(1,2);


%%
%% 6. display results
%%
if MODELDATA & strcmp(kernfmt,'space'),
   % load actual model cell STRF:
   load('/auto/k5/david/data/model/model.nr.simp.resp.mat','STRF');
   Hact=STRF.Phase0(:,:,16:-1:6);  % flip in time
   Hact=reshape(Hact,256,11);
   tkern=cat(3,Hact,SR,H);   % concatenate STRFs into a single 3D matrix
   titles={sprintf('actual strf (%s domain)',kernfmt),...
           sprintf('raw STA cc=%.2f (%s domain)',SRcc,kernfmt),...
           sprintf('bias-corrected STRF cc=%.2f (%s domain)',Hcc,kernfmt)};
else
   tkern=cat(3,SR,H);   % concatenate STRFs into a single 3D matrix
   titles={sprintf('cell %s raw STA cc=%.2f (%s domain)',...
                   cellid,SRcc,kernfmt),...
           sprintf('cell %s bias-corrected STRF cc=%.2f (%s domain)',...
                   cellid,Hcc,kernfmt)};
end

% display results
figure;
showkern(tkern,kernfmt,[movx movx],titles);



