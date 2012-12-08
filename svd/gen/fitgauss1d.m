% function [beta,beta0]=fitgauss1d(x,y,e,beta0,bcircular,batstr,maxeval);
%
% beta=[m s a d]
%
function [beta,beta0]=fitgauss1d(x,y,e,beta0,bcircular,batstr,maxeval);

if ~exist('e','var'),
   e=[];
end
if ~exist('beta0','var'),
   beta0=[];
end
if ~exist('bcircular','var'),
   bcircular=0;
end
if ~exist('batstr','var'),
   batstr='';
end
if ~exist('maxeval','var'),
   maxeval=5000;
end

VERBOSE=0;
opt=optimset('Display','off','MaxFunEvals',maxeval);
%opt=optimset('Display','iter');
%opt=[];
% fit spatial frequency tuning
scurve=y(:)';
smag=abs(scurve(min(find(abs(scurve)==max(abs(scurve)))))/100);
scurve=scurve./(smag+(smag==0));
scmin=min(scurve);
s0=mean(x(find(scurve-scmin>max(scurve-scmin)/2)));
if isnan(s0),
   s0=x(min(find(scurve==max(scurve)/2)));
end
beta0=[s0 ...
       mean(diff(x))*0.5*length(find( ...
          (scurve-min(scurve))>(max(scurve)-min(scurve))/2))...
       max(scurve)-scmin ...
       scmin];
%if beta0(2)<2,
%   beta0(2)=2;
%end
if bcircular,
   beta0(4)=0;
   betasflb=[min(x) min(diff(x))./2 0 0];
   betasfub=[max(x) inf max(scurve)*1.15 0+max(abs(y(:)))./100];
else   
   betasflb=[min(x) min(diff(x))./2 0 -Inf];
   betasfub=[max(x) inf max(scurve)*1.15 Inf];
end

if betasflb(3)==betasfub(3),
   beta=beta0;
   warning('no variance for fit');
   return
end

beta=lsqcurvefit('gauss1',beta0,x(:)',scurve,betasflb,betasfub,opt,bcircular);
beta0(3:4)=beta0(3:4).*smag;
beta(3:4)=beta(3:4).*smag;

beta=beta(:);
beta0=beta0(:);

if VERBOSE,
   fprintf('%s beta0: m=%.1f s=%.1f a=%.3f d=%.3f bcirc=%d\n',...
           batstr,beta0(1),beta0(2),beta0(3),beta0(4),bcircular);
   fprintf('%s beta : m=%.1f s=%.1f a=%.3f d=%.3f\n',...
           batstr,beta(1),beta(2),beta(3),beta(4));

   % display fit
   figure(1)
   clf
   plot(x,y);
   hold on
   plot(x,gauss1(beta0,x),'r:');
   plot(x,gauss1(beta,x),'r--');
   hold off
end
