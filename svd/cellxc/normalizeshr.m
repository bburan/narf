% function [H,lambda,pctvar]=normalizeshr(SR,sSA,tSA,nfactors,lambda,topSR,smoothtime)
% 
% created SVD 11/14/02 - hacked from normalizereg.m
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
function [H,lambda,pctvar]=normalizeshr(SR,sSA,tSA,nfactors,lambda,topSR,smoothtime);

disp('normalizeshr.m:');
global BATQUEUEID

TMINSFS=0.0001;  % hard coded minimum eigenvector for temporal
                 % decorr. hacky but it seems to work.

s=[size(SR) 1 1 ];
N=s(1);
U=s(2);
M=s(3);
R=s(4); % number of resamples
if ~exist('topSR','var'),
   topSR=1;
end

if size(sSA,1)~=N,
   disp('Spatial dimension must match in SR and sSA!');
   return
end

if exist('tSA','var') &size(tSA,1)>1, ,
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
   if ~exist('smoothtime','var'),
      smoothtime=0;
   end
end
if smoothtime>0,
   fprintf('smoothtime=%d\n',smoothtime');
end

if ~exist('nfactors','var'),
   nfactors=1;
   lambda=1;
else
   lambda=linspace(sqrt(5),sqrt(0.5),nfactors).^2;
end

H=zeros([N,U,nfactors,M,R]);

% as in XC, each response set is essentially independent here
for respidx=1:M,
   
   fprintf('respidx=%d/%d',respidx,M);
   
   tSR=zeros(N,U,R);
   for resampidx=1:R,
      
      % decorr time if tSA included
      if decorrtime,
         
         % normalize sSA to avoid double stimulus power
         % over-compensation. this is derivable from the idea that
         % tsSA(x1,t1,x2,t2)=sSA(x1,x2) * tSA(t1,t2)
         m=max(tSA(:,respidx,resampidx));
         
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
         
         tSR(:,:,resampidx)=SR(:,:,respidx,resampidx)*tSAinv';
         if smoothtime==1,
            tSR(:,:,resampidx)=conv2(tSR(:,:,resampidx),[0.2 0.6 0.2],'same');
         end
      else
         tSR(:,:,resampidx)=SR(:,:,respidx,resampidx);
         m=1;
      end
   end

   % spatial decorr
   tH=zeros(N,U,R);
   
   % figure out global eig space for shrinkage
   if respidx==1 | size(sSA,3)>1,
      fprintf('-->eigspace');
      SA=mean(sSA(:,:,respidx,:),4);
      [gU gS gV] = svd(SA);
      
      g=diag(gS);
      g=sqrt(g.^2./sum(g.^2));
      gSt=1./diag(gS);
      
      %gSt(find(g<0.0001))=0;
   end
   
   for resampidx=1:R,
      %if (resampidx==1 | R>1),
      %   %[tU tS tV] = svd(sSA(:,:,respidx,resampidx));
      %   %sD=gU'*tV*diag(1./diag(tS))*tU';
      %   st=diag(gU' * sSA(:,:,respidx,resampidx) * gU);
      %   g=sqrt(st.^2./sum(st.^2));
      %   
      %%   sst=[1./st(find(st>0)); ...
      %        zeros(length(find(st<=0)),1)];
      %   %sst(find(g<0.0001))=0;
      %   
      %   sD=diag(sst)*gU';
      %elseif respidx==1 | M>1,
      %   st=diag(gU' * sSA(:,:,respidx,resampidx) * gU);
      %   sst=[1./st(find(st>0)); ...
      %        zeros(length(find(st<=0)),1)];
      %   sD=diag(sst)*gU';
      %end
      
      st=diag(gU' * sSA(:,:,respidx,resampidx) * gU);
      g=sqrt(st.^2./sum(st.^2));
      sst=[1./st(find(st>0)); zeros(size(find(st<=0)))];
      sst(find(g<0.0005))=0;
      sD=diag(sst)*gU';
      
      tH(:,:,resampidx)=sD*tSR(:,:,resampidx);
      
   end
   dbsetqueue;
   
   Hm=mean(tH,3);
   Hs=std(tH,1,3) .* sqrt(R);
   Hs(find(Hs==0))=1;
   pctvar=zeros(nfactors,1);
   
   fprintf('-->shrink');
   
   for sfsidx=1:nfactors,
      Hdamp=(abs(Hm)./(Hs.*lambda(sfsidx)));
      pctvar(sfsidx)=mean(Hs(:)./(Hs(:).*lambda(sfsidx)));
      
      Hdamp=(1-Hdamp.^(-2));
      Hdamp=Hdamp.*(Hdamp>0);
      Hdamp(find(isnan(Hdamp)))=0;
      
      for resampidx=1:R,
         H(:,:,sfsidx,respidx,resampidx)=gU * (tH(:,:,resampidx).*Hdamp);
      end
      
   end
   dbsetqueue;
   
   fprintf('\n');
   
   if 0,
      %dump some junk out about RF
      PCUSE=40;
      tU=gU(:,1:PCUSE);
      for pp=1:size(tU,2),
         tU(:,pp)=tU(:,pp)./max(abs((tU(:,pp))));
         tU(:,pp)=tU(:,pp) .* sign(Hm(pp,1));
      end
      showstim(tU,'pfft',[18 18],4,10);
      
      figure
      subplot(5,1,1);
      semilogy(g)
      hold on
      semilogy(abs(gU'*tSR(:,:,1)),'r');
      hold off
      
      ttSR=mean((gU(:,:)'*tSR(:,:)),2);
      tmH=1./st.*mean((gU(:,:)'*tSR(:,:)),2);
      teH=1./st.*std((gU(:,:)'*tSR(:,:)),1,2) .* sqrt(R);
      
      subplot(5,1,2);
      plot([0 length(ttSR)+1],[0 0],'k--');
      hold on
      plot(ttSR,'r');
      plot(tmH);
      hold off
      
      subplot(5,1,3);
      plot([0 PCUSE+1],[0 0],'k--');
      hold on
      plot(ttSR(1:PCUSE),'r');
      errorbar(tmH(1:PCUSE),teH(1:PCUSE));
      hold off
      
      figure
      showkern(H(:,:,round(linspace(1,nfactors,6)),1),'pfft')
   end
   
end




