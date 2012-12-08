% function [mu,alpha,beta]=reversegamma(x,[alpha],[beta]);
%
% ryan's single trial gamma noise model
%
% x- vector of single trial data (integral number of spikes per bin)
% mu - randomly generated max likelihood modulator function
%
% use bayes rule to generate pdf for modulator mu, given observed
% response, x:
%
% p(mu|x)=p(x|mu) * p(mu) / p(x)
% p(x|mu) = poisson
% x is given
% p(mu) is a gamma distribution
%
% hyperparameters: use alpha and beta test ranges if specified
% otherwise estimate max likely alpha and beta from x
%
% created SVD 2004-04-26
%
function [mu,alpha,beta]=reversegamma(x,alpha,beta);

VERBOSE=0;

xlen=size(x,1);
nx=sum(~isnan(x),2);

gidx=find(nx>0);
x=x(gidx,:);
nx=nx(gidx);

if size(x,2)>1,
   sx=nansum(x')';
   mx=sx./nx;
else
   sx=x;
   mx=x;
end

ax=x(find(~isnan(x)));

if ~exist('alpha','var') | ~exist('beta','var'),
   arange=linspace(0.001,8,25);
   brange=linspace(0.001,2,25);
   
   fprintf('probing alpha=%.2f-%.2f (%d), beta=%.2f-%.2f (%d)\n',...
           arange([1 end]),length(arange),...
           brange([1 end]),length(brange));
   
   ee=zeros(length(arange),length(brange));
   for aa=1:length(arange),
      if mod(aa,5)==0 & VERBOSE,
         fprintf('%.3f ',arange(aa));
      end
      for bb=1:length(brange),
         
         % cost function for single trial
         data=(gamma(ax+arange(aa)) .* ...
               (brange(bb)./(1+brange(bb))).^(ax+arange(aa)) ...
               ./ (gamma(ax+1).*gamma(arange(aa)).*brange(bb).^arange(aa)));
         
         ee(aa,bb)=sum(log(data));
      end
   end
   
   mm=find(ee==max(ee(:)));
   [amax,bmax]=ind2sub([length(arange) length(brange)],mm);
   alpha=arange(amax);
   beta=brange(bmax);
   
   mu=gamrnd(alpha+sx,beta./(1+nx.*beta));
   
   if VERBOSE,
      fprintf('(alpha,beta)=(%.3f,%.3f)\n',alpha,beta);
      
      figure
      subplot(2,1,1);
      imagesc(arange,brange,ee);
      hold on
      plot(alpha,beta,'x');
      hold off
      colormap hot
   
      subplot(2,1,2);
      plot(sortrows([mx mu]));
   end
   
else
   mu=gamrnd(alpha+sx,beta./(1+nx.*beta));
   
end

mu0=mu;
mu=ones(xlen,1).*nan;
mu(gidx)=mu0;

