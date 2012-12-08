% function [H,lambda]=normalizesmooth(SR,sSA,tSA,nfactors,mintol,dim[=1])
% 
% created SVD 4/12/02 - hacked from normalize.m
% modified SVD 4/15/02 - scale lambdas according to eigenvalues
% 
% Inputs (basically the output from movxc.m or summed output from
%        multiple runs of it):
%      N - number of spatial dimensions
%      U - total number of time lags (eg, length([maxlag(1):maxlag(2)]))
%      M - number of response sets (or response dimensions)
% SR - raw spike-triggered average (correlation between stim and
%      resp) (N x U x M)
% sSA - spatial autocorrelation matrix (N x N x M or N x N and use
%       the same one for each response set--speeds things up)
% tSA - temporal autocorrelation matrix (U x U x M or leave empty
%       [], to skip temporal decorrelation)
% nfactors - number of regularization parameters to test (ranging
%            (10^3 through 10^-3)
% 
% Outputs:
% H - decorrelated kernel estimate(s), N x U x nfactors x M matrix
% lambda - regulatization parameters corresponding to each estimate
%          of H
% 
function [H,tolval]=normalizesmooth(SR,SA,nfactors,mintolval);

global BATQUEUEID

TMINSFS=0.0001;  % hard coded minimum eigenvector for temporal
                 % decorr. hacky but it seems to work.

s=size(SR);
N=s(1);
U=s(2);
M=prod(s(3:end));

if size(SA,1)==N & size(SA,2)==N,
   disp('normalizesmooth.m:  space');
   dim=1;
   stat=0;
elseif size(SA,1)==U & size(SA,2)==U,
   disp('normalizesmooth.m:  time');
   dim=2;
   stat=0;
elseif size(SA,1)==U*2-1 & size(SA,2)==1,
   disp('normalizesmooth.m:  time-stat');
   dim=2;
   stat=1;
else
   disp('unsupported SA/SR size combo');
   return
end
if dim==1,
   dimlen=N;
elseif dim==2,
   dimlen=U;
end

if length(nfactors)>1 | nfactors<=1,
   tolval=nfactors;
   nfactors=length(tolval);
elseif mintolval<=1,
   tolval=[1 10.^(linspace(-1,log10(mintolval),nfactors-1))];
else
   tolval=[1 (10.^(linspace(-1,-mintolval,nfactors-1)))];
end

H=zeros([N,U,nfactors,s(3:end)]);

% as in XC, each response set is essentially independent here
fprintf('respidx(/%d)=',M);
for respidx=1:M,
   
   fprintf('%d ',respidx);
   
   % decorr time if tSA included
   if stat,
      
      % create at toeplitz tSA matrix to get full temporal
      % decorr. seems kind of hokey but avoids edge effects for
      % short kernels
      tSA=zeros(dimlen);
      
      for uu=1:dimlen,
         tSA(:,uu)=SA((dimlen+1-uu):(dimlen*2-uu),1);
      end
   end
   
   % do the decorr
   % first precompute some intermediate matrices to speed things up
   % do some preliminary stuff to speed up decorrelation algebra
   
   % do SVD and parse out zero and non-zero eigenvalues
   if respidx==1 | size(SA,3)>1 | size(SA,4)>1,
      [sU sS sV] = svd(SA(:,:,respidx));
      sD=diag(sS);
      sDnz=sD(find(sD>0));
      sDz=sD(find(sD<=0));
      %sDpct=(sqrt(sD)./sum(sqrt(sD)));
      sDpct=((sD)./sum((sD)));
      %keyboard
      sDcumpct=1-cumsum(sDpct);
      sDi=1./sDnz;
   end
   
   % convert sta into svd domain.
   if dim==1,
      tSR=SR(:,:,respidx);
   else
      tSR=SR(:,:,respidx)';
   end

   tSR=sU'*tSR;

   for sfsidx=1:nfactors,
      
      % hackish thing: first guess is STA w/o spatial correlations removed
      if sfsidx==1,
         H(:,:,sfsidx,respidx)=sV*tSR;
      else
         % scale by regularization paramater here.
         tolscale=(sDcumpct>=tolval(sfsidx));
         borderidx=max([0; find(tolscale)])+1;
         
         pctbelow=sum(sDpct(1:borderidx-1));
         pctminusbelow=1-tolval(sfsidx)-pctbelow;
         
         tolscale(borderidx)=(sDpct(borderidx)-pctminusbelow)/sDpct(borderidx);
         if borderidx>1 & borderidx<dimlen,
            pctsave=tolscale(borderidx)*sDpct(borderidx+1)/2;
            tolscale(borderidx+1)=(pctsave)/sDpct(borderidx+1); 
            tolscale(borderidx-1)=...
                (sDpct(borderidx-1)-pctsave)/sDpct(borderidx-1);
         end
         
         SAspaceinv=sV*diag([(sDi.*tolscale);sDz]);
         H(:,:,sfsidx,respidx)=SAspaceinv*tSR;
      end
      
      if mod(sfsidx,10)==0,
         %fprintf('%d . ',sfsidx);
         dbsetqueue;
      end
   end
   
   %fprintf('\n');
end
fprintf('\n');




