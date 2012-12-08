% function [H,lambda,pctvar]=normalizereg(SR,sSA,tSA,nfactors,lambdaomag,topSR,smoothtime)
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
function [H,lambda,pctvar]=normalizereg(SR,sSA,tSA,nfactors,lambdaomag,topSR,smoothtime);

disp('normalizereg.m:');
global BATQUEUEID

TMINSFS=0.0001;  % hard coded minimum eigenvector for temporal
                 % decorr. hacky but it seems to work.

s=size(SR);
N=s(1);
U=s(2);
M=prod(s(3:end));
if ~exist('topSR','var'),
   topSR=1;
end

if size(sSA,1)~=N,
   disp('Spatial dimension must match in SR and sSA!');
   return
end

if exist('tSA','var') & size(tSA,1)>1,
   decorrtime=1;
   if size(tSA,1)~=U*2-1,
      disp('Temporal dimension must match in SR and tSA!');
      return
   end
   if ~exist('smoothtime','var'),
      smoothtime=1;
   end
else
   decorrtime=0;
   smoothtime=0;
end
if ~exist('lambdaomag','var') | isempty(lambdaomag),
   lambdaomag=4;
end
if ~exist('nfactors','var'),
   lambdaset=1;
   lambda=0;
   nfactors=1;
elseif exist('lambdaomag','var') & length(lambdaomag)>1,
   lambdaset=1;
   lambda=lambdaomag;
else
   lambdaset=0;  %generate lambda from SA matrix
end

if smoothtime>0,
   fprintf('(smoothtime=%d)',smoothtime');
end

H=zeros([N,U,nfactors,s(3:end)]);
ttSR=zeros([N,U,M]);
ttsD=zeros([N,M]);
ttsiD=zeros([N,nfactors,M]);

% as in XC, each response set is essentially independent here
fprintf('respidx(/%d)=',M);
for respidx=1:M,
   
   fprintf('%d ',respidx);
   
   % decorr time if tSA included
   if decorrtime,
      
      % normalize sSA to avoid double stimulus power
      % over-compensation. this is derivable from the idea that
      % tsSA(x1,t1,x2,t2)=sSA(x1,x2) * tSA(t1,t2)
      m=max(tSA(:,respidx));
      
      % create at toeplitz tSA matrix to get full temporal
      % decorr. seems kind of hokey but avoids edge effects for
      % short kernels
      ttSA=zeros(U);
      for uu=1:U,
         ttSA(:,uu)=tSA((U+1-uu):(U*2-uu),1);
         %ttSA(:,uu)=mean(tSA((U+1-uu):(U*2-uu),:),2);
      end
      
      % compute SVD on tSA and then pseudo-inverse
      [tU tS tV] = svd(ttSA);
      teigs=diag(tS)./sum(diag(tS));
      tSAinv=tV*diag(1./diag(tS).*(teigs>TMINSFS))*tU';
      
      tSR=SR(:,:,respidx)*tSAinv';
      if smoothtime==1,
         tSR=conv2(tSR,[0.2 0.6 0.2],'same');
      end
   else
      tSR=SR(:,:,respidx);
      m=1;
   end
   
   ttSR(:,:,respidx)=SR(:,:,respidx);
   
   % spatial decorr
   % first precompute some intermediate matrices to speed things up
   % do some preliminary stuff to speed up decorrelation algebra
   
   % do SVD and parse out zero and non-zero eigenvalues
   if respidx==1 | size(sSA,3)>1 | size(sSA,4)>1,
      %[sU sS sV] = svd(sSA(:,:,respidx)./m);
      [sU sS sV] = svd(sSA(:,:,respidx));
      sD=diag(sS);
      sDnz=sD(find(sD>0));
      sDz=sD(find(sD<=0));
   end
   
   % convert sta into svd domain.
   tSR=sU'*tSR;
   
   % save in mid-processing for display
   ttSR(:,:,respidx)=tSR;
   ttsD(:,respidx)=sD;
   
   if ~lambdaset,
      
      % scale reg parm to match stats of AC power
      % this case: first term scales so that first ac is equiv
      % to third???
      if topSR & length(sD)>=3,
         lambda=sD(1).*sD(3).*10.^(-[-1 linspace(0,lambdaomag,nfactors-1)]);
      elseif topSR,
         lambda=sD(1).*sD(2)./2.*10.^(-[-1 linspace(0,lambdaomag,nfactors-1)]);
      else
         % this case: first term scales so that first ac is equiv
         % to third???
         lambda=sD(1).*sD(3).*10.^(-linspace(0,lambdaomag,nfactors));
      end
      
      %fprintf(' (lambda(2)=%.2f) sfsidx=',lambda(2));
   end
   
   pctvar=zeros(nfactors,1);
   
   for sfsidx=1:nfactors,
      
      % hackish thing: first guess is STA w/o spatial correlations removed
      if lambda(sfsidx)>=lambda(1) & topSR,
         H(:,:,sfsidx,respidx)=sV*tSR;
         %fprintf('(topSR) ');
         if sum(sDnz)>0,
            pctvar(sfsidx)=(sDnz./sum(sDnz))' * sDnz./sum(sDnz);
         end
      else
         % scale by regularization parameter here.
         sDi=1./(sDnz+lambda(sfsidx)./sDnz);
         
         % option: force max sDi >= sDi(1)
         sDi(find(sDi<sDi(1)))=sDi(1);
         
         pctvar(sfsidx)=(sDnz.*sDi)' * sDnz./sum(sDnz);
         
         ttsiD(:,sfsidx,respidx)=[sDi; sDz];
         
         SAspaceinv=sV*diag([sDi;sDz]);
         H(:,:,sfsidx,respidx)=SAspaceinv*tSR;
         if smoothtime==2,
            H(:,:,sfsidx,respidx)=conv2(H(:,:,sfsidx,respidx),...
                                        [0.2 0.6 0.2],'same');
         end
      end
      
   end
   if 0,
      %fprintf('%d . ',sfsidx);
      dbsetqueue;
   end
   
end
fprintf('\n');

%keyboard

if 0,
   setmap(redblue);
   
   %dump some junk out about RF
   PCUSE=100;
   tU=sU(:,1:PCUSE);
   
   mtSR=mean(ttSR,3);
   ttH=ttSR./repmat(reshape(ttsD,[N 1 M]),[1 U]);
   mtH=mean(ttH,3);
   etH=std(ttH,1,3) .* sqrt(M);
   
   msD=mean(ttsD,2);
   g=sqrt(msD.^2./sum(msD.^2));
   
   for pp=1:size(tU,2),
      tU(:,pp)=tU(:,pp)./max(abs((tU(:,pp))));
      tU(:,pp)=tU(:,pp) .* sign(mtH(pp,1));
   end
   showstim(tU,'pfft');
   
   figure
   subplot(3,3,1);
   semilogy(g,'r')
   hold on
   semilogy(abs(mtSR),'b');
   hold off
   a=axis;
   axis([0 PCUSE+1 a(3:4)]);
   title('stimpow/sta');
   
   subplot(3,3,2);
   plot([0 length(mtSR)+1],[0 0],'k--');
   hold on
   plot(mtSR,'b');
   plot(mtH,'k');
   hold off
   a=axis;
   axis([0 PCUSE+1 a(3:4)]);
   
   subplot(3,3,3);
   plot([0 PCUSE+1],[0 0],'k--');
   hold on
   plot(ttSR(1:PCUSE),'b');
   errorbar(mtH(1:PCUSE),etH(1:PCUSE),'k');
   hold off
   a=axis;
   axis([0 PCUSE+1 a(3:4)]);
   title('hest no reg');
   
   sfsrange=round(linspace(2,nfactors,6));
   
   % plot eig domain kernels with sucessively weaker regularization
   for ssi=1:length(sfsrange),
      sfsidx=sfsrange(ssi);
      subplot(3,3,3+ssi);
      plot([0 PCUSE+1],[0 0],'k--');
      hold on
      plot(ttsD(1:PCUSE,1).*ttsiD(1:PCUSE,sfsidx,1),'r');
      plot(ttSR(1:PCUSE,1).*ttsiD(1:PCUSE,sfsidx,1),'k');
      %errorbar(mtH(1:PCUSE),etH(1:PCUSE));
      hold off
      a=axis;
      axis([0 PCUSE+1 -1 1]);
      title(sprintf('hest %.3f %% var',pctvar(sfsidx).*100));
   end
   
   % Hbad means no regularization
   Hbad=sV*mtH;
   
   figure
   showkern(cat(3,H(:,:,sfsrange,1),Hbad),'pfft')
end




