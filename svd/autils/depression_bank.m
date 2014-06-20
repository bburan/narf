% function dstim=depression_bank(stim,u,tau,fixu[=0],crosstalk[=0])
%
% Apply the Tsodyks & Markram synaptic depression model to a
% spectrogram or some other collection of stimulus envelopes.
%
% On each time step,
%  dstim(t) = stim(t) * d(t)
% where d(t) is the amount of depression from previous inputs,
%  d(t)=d(t-1) - stim(t-1)*d(t-1)*u + (1-d(t-1))/tau
% (depression accumulates/recovers on each spectral channel separately)
%
% inputs:
%  stim -spectrogram freq X time
%  u - vector of depression strengths as fraction of max of stim
%  tau - vector of recovery time constants in units of time bins
%  fixu - (default 0) if 1, don't normalize u to max of stim, just
%         use absolute value
%  crosstalk - (default 0) calculate depression using <crosstalk>% of
%              neighboring spectral channels rather than just the 
%              current channel
%
% returns:
%  dstim - depressed stimulus (freq*length(u)) X time matrix. ie, if
%          multiple depression synapses specified, stack the output
%          matrices for the different synapses on top of each other.
% 
% created svd 2010-07-08
%
function dstim=depression_bank(stim,u,tau,fixu,crosstalk)

global ESTIMATIONPHASE
global DEPSTIMMAX

if ~exist('fixu','var'),
   fixu=0;
end
if ~exist('crosstalk','var'),
   crosstalk=0;
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
      if length(crosstalk)>1,
          tstim=crosstalk;
      elseif crosstalk>0,
          % 1 channel do nothing
          if size(stim,1)==2,
              sfilt=[1-crosstalk./100 crosstalk./100;
                     crosstalk./100 1-crosstalk./100];
              tstim=sfilt*tstim;
          elseif size(stim,1)>2,
              sfilt=[crosstalk./100; 1-crosstalk.*2./100; crosstalk./100];
              tstim=rconv2(tstim,sfilt);
          end
          %keyboard
      end
      %tstim(:,1:10)=min(tstim(:));
      
      if taui(1)./ui(1)./10>1,
          subsample=1;
      elseif taui(1)./ui(1)./10>.1,
          subsample=10;
      elseif taui(1)./ui(1)./10>.01,
          subsample=100;
      elseif ui(1)>10,
          % tau too small, effectively no depression.
          % otherwise this will create unstable oscillations
          ui(:)=taui(:)./10*1000;
          subsample=100;
      else
          % tau too small, effectively no depression.
          % otherwise this will create unstable oscillations
          subsample=1;
          ui(:)=0;taui(:)=10;
      end
      taui_ss=taui.*subsample;
      ui_ss=ui./subsample;
      for ii=2:size(stim,2),
          if 1 % subsample as much as is necessary to avoid oscillations
              td=di(:,ii-1);
              for dd=1:subsample,
                  delta=(1-td)./taui_ss - ui_ss.*td.*tstim(:,ii-1);
                  td=td+delta;
                  td(td<0)=0;
              end
              di(:,ii)=td;
          else
              delta=(1-di(:,ii-1))./taui - ui.*di(:,ii-1).*tstim(:,ii-1);
              di(:,ii)=di(:,ii-1)+delta;
              di(di(:,ii)<0,ii)=0;
              %di(di(:,ii)>1,ii)=1;
          end
          
      end
      
      dstim((1:size(stim,1))+(jj-1)*size(stim,1),:)=di.*stim;
   else
      dstim((1:size(stim,1))+(jj-1)*size(stim,1),:)=stim;
   end
end


