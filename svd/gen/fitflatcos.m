% function [beta,beta0]=fitflatcos(x,y,e,beta0,batstr);
%
function [beta,beta0]=fitflatcos(x,y,e,beta0,batstr);

if ~exist('e'),
   e=[];
end
if ~exist('beta0'),
   beta0=[];
end
if ~exist('batstr'),
   batstr='';
end

opt=optimset('Display','off');      

ocurve=y(:)';
e=e(:)';
%omag=ocurve(min(find(abs(ocurve)==max(abs(ocurve)))))/100;
omag=abs(ocurve(min(find(abs(ocurve)==max(abs(ocurve)))))/100);
ocurve=ocurve./(omag+(omag==0));
ocmin=min(ocurve);
%keyboard
m0=mean(x(find(ocurve-ocmin>=max(ocurve-ocmin)/2)));
if isnan(m0),
   m0=x(min(find(ocurve==max(ocurve)/2)));
end

if length(beta0)<4,
   beta0=[m0 ...
          mean(diff(x))*0.5*(length(x)-length(find( ...
             (ocurve-min(ocurve))>(max(ocurve)-min(ocurve))/2)))...
          max(ocurve)-ocmin ...
          ocmin];
else
   beta0=[beta0(1) ...
          beta0(2) ...
          max(ocurve)-ocmin ...
          ocmin];
end
betaorlb=[-Inf 0 0 -Inf];
betaorub=[Inf 120 Inf Inf];
if isempty(e),
   beta=lsqcurvefit('flatcos',beta0,x,ocurve,betaorlb,betaorub,opt);
else
   beta=lsqcurvefit('flatcos',beta0,x,ocurve./e,betaorlb,betaorub,opt,e);
end

beta(1)=mod(beta(1),180);

beta0(3:4)=beta0(3:4).*omag;
beta(3:4)=beta(3:4).*omag;

fprintf('%s   or0: m=%.1f s=%.1f a=%.3f d=%.3f\n',...
        batstr,beta0(1),beta0(2),beta0(3),beta0(4));
fprintf('%s   or : m=%.1f s=%.1f a=%.3f d=%.3f\n',...
        batstr,beta(1),beta(2),beta(3),beta(4));

% force params to be returned in column format
beta=beta(:);
beta0=beta0(:);

return

% display fit
figure(1)
clf
errorbar(x,y,e);
hold on
plot(x,flatcos(beta,x),'r--');
hold off
