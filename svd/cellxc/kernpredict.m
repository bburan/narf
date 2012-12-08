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
function [r]=kernpredict(sfIR,mov,phasecount,brect,convtime,zerobin)

if not(exist('phasecount','var')),
   phasecount=4;
end
if not(exist('brect','var')),
   brect=1;
end
if not(exist('convtime','var')),
   convtime=1;
end
if ~exist('zerobin','var'),
   zerobin=1;
end

sfIRtot=sum(conj(sfIR(:)).*sfIR(:));
movlen=size(mov,2);
upc=size(sfIR,1)/phasecount;

if convtime,
   % flip kernel to reverse time to speed convolution
   
   binsize=size(sfIR,2);
   
   r=zeros(movlen-binsize+1,phasecount);
   if ~sfIRtot,
      if brect,
         r=zeros(movlen-binsize+1,1);
      end
      r=[ones(binsize-1,size(r,2))*nan; r];
      return
   end
   
   
   if upc==1,
      for phaseidx=1:phasecount,
         pr=conv(mov(phaseidx,:),sfIR(phaseidx,:));
         r(:,phaseidx)=pr(binsize:movlen)';
      end
      r=real(r); % same as ( mov X sfIR + conj(mov) X conj(sfIR))./2
      if brect,
         r(find(r<0))=0;
      end
   else
      for phaseidx=1:phasecount,
         
         % figure out relevant stim channels for this phase set
         spacerange=upc*(phaseidx-1)+1:upc*phaseidx;
         
         % compute inner product between stim and STRF for all latencies
         tr=sfIR(spacerange,:)'*mov(spacerange,:);
         
         % sum over latencies
         pr=tr(1,binsize:movlen);
         for tt=2:binsize,
            pr=pr+tr(tt,binsize-tt+1:movlen-tt+1);
         end
         
         % take real part (necessary??)
         pr=real(pr);
         
         % rectify
         if brect,
            pr(find(pr<0))=0;
         end
         
         % save output for this phase channel
         r(:,phaseidx)=pr';
      end
   end
   
   
else
   % do stuff for fixation-triggered pred...
   segcount=size(sfIR,2);
   binsize=1;
   r=zeros(movlen,phasecount,segcount);
   if ~sfIRtot,
      return
   end
   for phaseidx=1:phasecount,
      spacerange=upc*(phaseidx-1)+1:upc*phaseidx;
      
      tr=mov(spacerange,:)' * sfIR(spacerange,:);
      r(:,phaseidx,:)=reshape(tr,movlen,1,segcount);
   end
   
   % take real part???
   r=real(r);
   
   % rectify:
   if brect,
      r(find(r<0))=0;
   end
end

if brect,
   % square before summing over phase channels if brect=1.
   %r=squeeze(sum(r,2).^2);
   r=squeeze(sum(r,2));
end

%r=[-ones(binsize-1,size(r,2)); r];
%r=[zeros(binsize-1,size(r,2)); r];

r=[ones(binsize-zerobin,size(r,2))*nan; r; ones(zerobin-1,size(r,2))*nan];


