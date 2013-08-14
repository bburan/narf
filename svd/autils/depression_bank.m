% function dstim=depression_bank(stim,u,tau,fixu)
%
% on each time step,
%  dstim(t) = stim(t) * d(t)
% where d(t) is the amount of depression from previous inputs,
%  d(t)=d(t-1) - stim(t-1)*d(t-1)*u + (1-d(t-1))/tau
% (depression accumulates/recovers on each spectral channel separately)
%
% inputs:
%  stim -spectrogram freq X time
%  u - depression strength as fraction of max of stim
%  tau - recovery time constant in units of time bins
%  fixu - (default 0) if 1, don't normalize u to max of stim, just
%         use absolute value
% 
% created svd 2010-07-08
%
function dstim=depression_bank(stim,u,tau,fixu)

global ESTIMATIONPHASE
global DEPSTIMMAX

if ~exist('fixu','var'),
   fixu=0;
end
if ~exist('verbose','var'),
   verbose=0;
end

dcount=length(tau);
dstim=zeros(dcount*size(stim,1),size(stim,2));
internalfs=100;
istimthresh=0;

if verbose,
   fprintf('Applying depression_bank\n');
end

% don't update if we KNOW we're in the validation phase of cellxcnodb.m
if isempty(ESTIMATIONPHASE) || ESTIMATIONPHASE || isempty(DEPSTIMMAX),
   DEPSTIMMAX=max(stim,[],2);
end

for jj=1:dcount,
   
   if ~fixu,
       ui=u(jj)./DEPSTIMMAX(:);
   else
       ui=u(:,jj);
   end
   
   taui=tau(jj);
   %taui=floor(taui);
   if taui==0,
       ui=0;
   end
   
   taui=abs(taui);
   ui=abs(ui);
   
   if verbose,
      fprintf(' efficacy = %.3f (max(stim)=%.3f)\n',ui,DEPSTIMMAX);
      fprintf(' time constant = %.1f bins (%.3f sec?)\n',...
              taui,taui./internalfs);
   end
   
   if any(ui),
      % apply threshold to stim going into adaptation mechanism
      tstim=(stim>istimthresh).*stim - istimthresh;
      tstim=max(tstim,0);
      di=ones(size(stim));
      if size(tstim,2)<10,
         keyboard
      end
      tstim(:,1:10)=min(tstim(:));
      for ii=2:size(stim,2),
         delta=(1-di(:,ii-1))./taui - ui.*di(:,ii-1).*tstim(:,ii-1);
         di(:,ii)=di(:,ii-1)+delta;
         di(di(:,ii)<0,ii)=0;
         %di(di(:,ii)>1,ii)=1;
      end
      
      dstim((1:size(stim,1))+(jj-1)*size(stim,1),:)=di.*stim;
   else
      dstim((1:size(stim,1))+(jj-1)*size(stim,1),:)=stim;
   end
end

%if size(stim,1)>1,
%    dstim=reshape(dstim,size(stim,1),dcount,size(stim,2));
%    dstim=permute(dstim,[2 1 3]);
%    dstim=reshape(dstim,dcount*size(stim,1),size(stim,2));
%end

