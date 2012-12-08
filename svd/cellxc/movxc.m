% function [SR,n,mS,mR,tSA,sSA2]=movxc(stim,resp,maxlag,boundary,
%                                           singlesSA,meansub)
%
% calculates spike-triggered average and stimulus autocorrelation
% matrices given a stimulus and response
%
% Inputs:
% stim - T x N  matrix where space is unwrapped into a single
%        dimension of size N. T is the number of time bins
% resp - T x M  matrix where M are different responses to the same
%        stimulus (eg, resamplings of a neuron, multiple neurons,
%        voxels, etc)
%        Nan values in resp are considered invalid bins for
%        the given response vector and will be skipped in the
%        various calculations.
% maxlag - range of latencies to include in kernel. 
%          maxlag=[t0 t1] calculate SR over latancies t0:t1
%          if maxlag is scalar, use t0=-maxlag, t1=maxlag
% boundary - 'circular' - assume stimulus is circular
%            'zero' - assume zero values outside of valid times (default)
% singlesSA - 0 [default] - separate sSA for each response channel
%             1 - single sSA for all channels (faster, use for PFTH)
% meansub - subtract mean from stim and resp before calculating the
%           various correlation matrices (default 1)
%
% Outputs:
% SR   - stimulus-response cross correlation (N x maxlag*2+1 x M)
% n    - number of samples in each response set (M x 1)
% tSA  - temporal stimulus autocorrelation (maxlag*2+1 x maxlag*2+1 x M)
% sSA2 - 2nd order spatial stimulus autocorrelation (N x N x M)
% mS   - stimulus mean for each response channel (N x M)... may
%        differ along m dimension because some response channels
%        may be invalid at different times
% mR   - mean response for each channel (M x 1)
%
% created SVD 3/20/02 - hacked from movcorrs.m
% modified SVD 7/12/02 - removed sSA1, added mS,mR
%
function [SR,n,mS,mR,tSA,sSA2]=movxc(stim,resp,maxlag,boundary,...
                                     singlesSA,meansub)

if nargin<2,
   disp('SYNTAX: movxc(stim,resp,maxlag,boundary,singlesSA,meansub)');
   disp('First two parameters required');
   return
end

if not(exist('maxlag')),
   maxlag=[-10 10];
end
if length(maxlag)==1,
   minlag=-maxlag;
else
   minlag=maxlag(1);
   maxlag=maxlag(2);
end
if not(exist('boundary')),
   boundary='zero';
end
if not(exist('singlesSA')),
   singlesSA=0;
end
if not(exist('meansub')),
   meansub=1;
end

% figure out spatial and temporal sizes
N=size(stim,2);
M=size(resp,2);
T=size(stim,1);
if size(resp,1)~=T,
   disp('Time dimension must match in stim and resp!');
   return
end
timebincount=(maxlag-minlag)+1; % total number of latency bins (0 +/- maxlag)

% if boundary is 'zero', add on some padding to stim and resp
%if strcmp(boundary,'zero'),
%   stim=[zeros(timebincount,N); stim];
%   resp=[nan*ones(timebincount,M); resp];
%   T=T+timebincount;
%end

if nargout>4,
   dosSA2=1;
else
   dosSA2=0;
end

% create output matrices
SR=zeros(N,timebincount,M);
if singlesSA,
   R=1;
else
   R=M;
end
if dosSA2,
   sSA2=zeros(N,N,R);
end
tSA=zeros(timebincount*2-1,R);
n=zeros(M,1);

if meansub,
   % subtract mean from stim
   mS=mean(stim,1); % mean stim for whole sequence
   for nn=1:N,
      stim(:,nn)=stim(:,nn)-mS(nn);
   end
end

mR=zeros(M,1);     % mean response, averaged over valid time bins
mS=zeros(N,M);     % mean stim, averaged over valid time bins for a
                   % given response channel

tzero=1-minlag;
ttzero=maxlag-minlag+1;

% calculate outputs for each response set
for respidx=1:M,
   
   if tzero>1,
      fprintf('%d/%d:',respidx,M);
   end
   
   % valid time bins for this set
   if strcmp(boundary,'zero'),
      tgood=(find(~isnan(resp((maxlag+1):(end+minlag),respidx)))+maxlag)';
   else
      tgood=find(~isnan(resp(:,respidx)))';
   end
   n(respidx)=length(tgood);
   
   % compute increments to mean stim and resp
   mR(respidx)=sum(resp(tgood,respidx));
   mS(:,respidx)=sum(stim(tgood,:),1)';  % should be approximately 0
   
   if length(tgood)>0 & meansub,
      % subtract mean response from current response set
      resp(tgood,respidx)=resp(tgood,respidx)-mR(respidx)./length(tgood);
      mR(respidx)=sum(resp(tgood,respidx));  % should be 0 now
   end
   
   % step through lags and calculate SR and tSA as needed
   %for tt=minlag-maxlag:maxlag,
   for tt=minlag:maxlag,
      tgoodt=mod(tgood-tt-1,T)+1; % offset for time lag and take modulus
      
      %if tt>=minlag & tt<=maxlag,
         SR(:,tzero+tt,respidx)=stim(tgoodt,:)'*resp(tgood,respidx);
      %end
      
      % moved tSA calc to separate loop to avoid edge effects.
      %if (~singlesSA | respidx==1) & tt<=0,
      %   tSA(ttzero+tt,respidx)=sum(sum(stim(tgoodt,:).*stim(tgood,:)))./N;
      %   tSA(ttzero-tt,respidx)=tSA(ttzero+tt,respidx);
      %end
      
      if tt<0
         fprintf('-');
      elseif tt==0
         fprintf('0');
      else
         fprintf('+');
      end
   end
   
   if ~singlesSA | respidx==1,
      
      fprintf('T');
      for tt=1:N,
         tSA(:,respidx)=tSA(:,respidx)+xcorr(stim(:,tt),ttzero-1,...
                                             'unbiased')./N .* size(stim,1);
      end
      
      % calculate 2nd order spatial autocorrelation matricies
      % and sSA2 (if selected for output)
      % for output
      fprintf('S');
      if dosSA2,
         sSA2(:,:,respidx)=conj(stim(tgood,:))' * stim(tgood,:);
      end
      
   end
   if tzero>1 | respidx==M,
      fprintf('\n');
   end
end






