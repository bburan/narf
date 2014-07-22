% function dstim=adp_bank(stim,u,tau,fixu[=0],crosstalk[=0])
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
%     or a matrix of process_count X synapse_count strengths
%  tau - vector of recovery time constants in units of time bins
%     or a matrix of process_count X synapse_count time constants
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
function dstim=adp_bank(stim,u,tau,fixu,crosstalk)

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

proc_count=size(tau,1);
syn_count=size(tau,2);
dstim=zeros(syn_count*size(stim,1),size(stim,2));

if verbose,
    fprintf('Applying depression_bank\n');
end

% don't update if we KNOW we're in the validation phase of cellxcnodb.m
if isempty(ESTIMATIONPHASE) || ESTIMATIONPHASE || isempty(DEPSTIMMAX),
    DEPSTIMMAX=max(stim,[],2);
end

for jj=1:syn_count,
    
    if ~fixu,
        ui=(1./DEPSTIMMAX(:))*u(:,jj);
    else
        ui=u(:,jj);
    end
    
    taui=tau(:,jj);
    ui(taui==0)=0;
    
    taui=abs(taui); % no negative time constants
    
    if verbose,
        fprintf(' efficacy = %.3f (max(stim)=%.3f)\n',ui,DEPSTIMMAX);
        fprintf(' time constant = %.1f bins (%.3f sec?)\n',...
                taui,taui./100);
    end
    
    if any(ui),
        % positively rectify stim going into adaptation mechanism
        tstim=(stim>0).*stim;
        di=ones(size(stim));
        if size(tstim,2)<10,
            disp('tstim too short');
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
        end
        %tstim(:,1:10)=min(tstim(:));
        
        % tamp down ui if large to prevent oscillations
        for ii=1:length(taui),
            if taui(ii)/abs(ui(ii))<10,
                ui(ii)=sign(ui(ii))*taui(ii)/10;
            end
        end
        
        for pp=1:proc_count,
            tdi=ones(size(stim));
            for ii=2:size(stim,2),
                td=tdi(:,ii-1);
                if ui(pp)>0
                    delta=(1-td)./taui(pp) - ui(pp).*td.*tstim(:,ii-1);
                    td=td+delta;
                    td(td<0)=0;
                else
                    delta=(1-td)./taui(pp) - ui(pp).*td.*tstim(:,ii-1);
                    td=td+delta;
                    td(td<1)=1;
                end
                tdi(:,ii)=td;
            end
            % accumulate effects
            di=di.*tdi;
        end
        
        dstim((1:size(stim,1))+(jj-1)*size(stim,1),:)=di.*stim;
    else
        % ui==0, no dep/fac, just pass through input
        dstim((1:size(stim,1))+(jj-1)*size(stim,1),:)=stim;
    end
end


