% xccorefet.m
%
% script to compute unbiased kernel from preprocessed stimulus and
% response matrices using FET's 2001 algorithm
%
% created SVD 8/24/06 - ripped off of xccore and strfpak
%

disp('xccorefet.m:');

%
% define resampling regimes here
%
% option to branch according to resampfmt
if params.resampcount==1 | length(resp)==0,
   % bootstrap resampling
   rstartidx=1;
   rendidx=size(resp,1);
elseif params.resampfmt==1,
   [rstartidx,rendidx]=resampsegs(resp,params.resampcount);
   
elseif params.resampfmt==2,
   % shuffled response resampling
   rstartidx=starttimes(fidx,attidx);
   rendidx=stoptimes(fidx,attidx);
   resp=resampnoise(resp,params.resampcount);
   rsize=[rsize,size(resp,3)];
   resp=resp(:,:); % reshape to 2D matrix
   
elseif params.resampfmt==3,
   % gamma/poisson spiking model resampling   
   rstartidx=1;
   rendidx=size(resp,1);
   respcount=1;
   rsize(2)=1;
   [mu,alpha,beta]=reversepoisson(resp);
end


if firstseg,
   disp('initializing matrices');
   if 0 & (diff(params.maxlag)==0 | rsize(2)==1),
      corrmtxcount=1;
      singlesSA=1;
   else
      corrmtxcount=rsize(2);
      singlesSA=0;
   end
   
   % first time, set kernel & ac matrices to zero
   %SR=zeros(spacecount,diff(params.maxlag)+1,respcount,params.resampcount);
   n=zeros(respcount,params.resampcount);
   mS=zeros(spacecount,respcount,params.resampcount);
   mR=zeros(respcount,params.resampcount);
   %tSA=zeros(diff(params.maxlag)*2+1,corrmtxcount,params.resampcount);
   %sSA2=zeros(spacecount,spacecount,corrmtxcount,params.resampcount);
   %if params.decorrspace==4,
   %   sSAfull=zeros(spacecount*(diff(params.maxlag)+1),spacecount*(diff(params.maxlag)+1),...
   %       params.resampcount);
   %end
   
   % global stimulus temporal autocorr. works better than resamped?
   %tSA1=zeros(diff(params.maxlag)*2+1,corrmtxcount,params.resampcount);
   %nSA1=0;
   firstseg=0;

   params.maxlag(1)=-params.maxlag(2);
   twindow=params.maxlag;
   nband=spacecount;
   
   % temporal axis range
   tot_corr = diff(twindow) + 1;
   
   % spatial axis range
   spa_corr = (nband * (nband - 1))/2 + nband;
   
   % initialize CS and CSJN variable
   CS = zeros(spa_corr, tot_corr, params.resampcount);
   CS_ns = zeros(1, tot_corr, params.resampcount);
   
   % initialize autoCorr and autoCorrJN variable
   CSR = zeros(nband, tot_corr, params.resampcount,respcount);
   CSR_ns = zeros(1,tot_corr, params.resampcount,respcount);
end

fprintf('resamp=');
for resampidx=1:params.resampcount,
   fprintf('%d ',resampidx);
   for respidx=1:respcount,
      if params.resampfmt~=3,
         tr=resp(:,respidx);
         %tr([1:rstartidx(resampidx,respidx)-1 ...
         %    rendidx(resampidx,respidx)+1:end],:)=nan;
         tr(rstartidx(resampidx,respidx):rendidx(resampidx,respidx),:)=nan;
      else
         tr=reversepoisson(resp(:,2:end),alpha,beta);
      end
      
      rgoodidx=find(~isnan(tr));
      rbadidx=find(isnan(tr));
      
      nlen=length(tr);
      
      stimval=stim';
      mS(:,respidx,resampidx)=mean(stimval(:,rgoodidx),2);
      stimval=stimval-repmat(mS(:,respidx,resampidx),[1,nlen]);
      stimval(:,rbadidx)=0;
      
      % ========================================================
      % AC calculation
      % The algorithm is based on FET's dcp_stim.c
      % ========================================================
      % mod svd 2010-02-13, only calc for respidx==1
      if respidx==1,
         xb=1;
         for ib1 = 1:nband
            for ib2 = ib1:nband
               
               % NEW version of algorithm by using xcorr 
               CS(xb, :, resampidx) = CS(xb, :, resampidx)+ ...
                   xcorr(stimval(ib1,:),stimval(ib2, :),twindow(2));
               xb = xb +1;
            end              % END of ib2
         end                 % END of ib1
         
         % Count the total trials for later normalization 
         lengthVec = ones(1, length(rgoodidx));
         CS_ns(1, :, resampidx) = CS_ns(1, :, resampidx)+ ...
             xcorr(lengthVec, lengthVec, twindow(2));
      end
      
      % ========================================================
      % cross-corr. The algorithm is from FET's dcp_stim_spike
      % ========================================================
      
      mR(respidx,resampidx)=nanmean(tr);
      psthval = tr-mR(respidx,resampidx);
      psthval(isnan(psthval))=0;
      
      % New version of algorithm for computing cross-correlation 
      for ib1 = 1:nband
         CSR(ib1, :,resampidx,respidx) = CSR(ib1, :,resampidx,respidx)+ ...
             xcorr(stimval(ib1, :), psthval, twindow(2));
      end
      
      % For normalization and assign the count_ns
      CSR_ns(:,:,resampidx,respidx) = CS_ns(1, :, resampidx,1);
   end
   
   % update queue if active
   %dbsetqueue;
end
fprintf('\n');

clear stimval


% match H to H in xccore
H=zeros(nband,tot_corr,params.sfscount,respcount,params.resampcount);

for resampidx=1:params.resampcount
   for respidx=1:respcount,
      % take temporal FFT of stim autocorrelation and stim-resp
      % cross-correlation
      % [] is to skip jackknife
      nstd_val = 0.5;
      [fstim, fstim_spike, stim_spike_JNf] = ...
          fft_AutoCrossCorr(CS(:,:,resampidx),CSR(:,:,resampidx,respidx),...
                            [],twindow(2), nband, nstd_val);
      
      nb = nband;
      nt = 2*twindow(2) +1;
      nJN = 1;
      stim_size = size(fstim);
      stim_spike_size = size(fstim_spike);
      stim_spike_JNsize = size(stim_spike_JNf);
      
      % set tolerance values
      Tol_val=10.^(-linspace(1,params.sfsstep+1,params.sfscount));
      ntols = params.sfscount;
      lambda=Tol_val;
      
      for itol=1:ntols
         tol=Tol_val(itol);
         
         % ======================================= 
         % Calculate strf for each tol val.        
         % ======================================= 
         [forward, forwardJN, forwardJN_std] = ...
             cal_Strf(fstim,fstim_spike, stim_spike_JNf,stim_size, ...
                      stim_spike_size,stim_spike_JNsize, nb, nt, nJN, tol);
         
         H(:,:,itol,respidx,resampidx) = forwardJN;
      end
      fprintf('Done calculation of STRF for resampidx=%d\n', resampidx);
      
   end
end

% hacked in from xccore:

if 0,
   irS=squeeze(std(SR,1,1));
   irH=squeeze(std(H,1,1));
   
   allt=std(reshape(permute(H,[1 2 4 5 3]),spacecount*(diff(params.maxlag)+1)* ...
                    respcount*params.resampcount,params.sfscount));
   alls=std(reshape(permute(H,[1 3 4 5 2]),spacecount*params.sfscount* ...
                    respcount*params.resampcount,(diff(params.maxlag)+1)));
   tH=squeeze(std(reshape(permute(H,[1 4 5 2 3]),spacecount* ...
                          respcount*params.resampcount,...
                       (diff(params.maxlag)+1),params.sfscount)));
end

meanCS=mean(CS,3);
clear CS
%pack

mH=mean(H,5);
eH=std(H,1,5) .* sqrt(params.resampcount-1);

disp('created mH');

mSR=mean(CSR,3);
eSR=std(CSR,1,3).* sqrt(params.resampcount-1);

disp('created mSR');

clear CSR

% powunbiased not defined.  can define it somehow?
%   powunbiased=normalizereg(powbiased,mSA2,[],params.sfscount,...
%                            params.sfsstep,1);

% take mean stim/response across resamples
mSall=mean(mS,3);
mRall=mean(mR,2);
ntot=n;
disp('xccorefet.m done');









