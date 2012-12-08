% function [r]=kernpredict(sfIR,mov,phasecount[=4],brect[=1],convtime[=1],zerobin[=1])
%
% sfIR      - spacecount X tbincount : array of kernels
%             (forward time--gets flipped here before running)
% mov       - spacecount X time :  pre-loaded movie matrix
% phasecount- number of separately rectified phase channels (assume
%             they're all equal size and stacked in the first
%             dimension for compatibility) (default=4)
% brect     - rectify indiv phase channels 0=no, 1=yes (default=1).
% convtime  - convolve in time (1) or just do spatial inner prod
%             (0) (default=1)
%
function [r]=kernpredict_2nd(H,stim,T1,T2)

% debug 2nd order prediction code 
stimlen=size(stim,1);
r=zeros(stimlen,1);

shiftstim=zeros(stimlen,T2-T1+1);
for os=0:(T2-T1),
   % t2 is always >= t1
   shiftstim((1+os):end,os+1)=stim(1:(end-os),:).*...
       stim((1+os):end,:);
end

for t1=T1:T2,
   gidx=max(t1+1,1):min(t1+stimlen,stimlen);
   r(gidx)=r(gidx)+stim(gidx-t1,:)*H(1,t1-T1+1,end);
   
   for t2=t1:T2,
      os=t2-t1;
      if os>0
         r(gidx)=r(gidx)+shiftstim(gidx-t1,os+1)* 2 *...
            H(t2-T1+2,t1-T1+1,end);
      else
         r(gidx)=r(gidx)+shiftstim(gidx-t1,os+1)*...
            H(t2-T1+2,t1-T1+1,end);
      end
   end
end

