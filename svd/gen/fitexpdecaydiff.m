% function [beta,beta0]=fitexpdecaydiff(x,y,e,beta0,batstr);
%
function [beta,beta0]=fitexpdecaydiff(x,y,e,beta0,batstr);

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

% fit positive temporal IR
tcurve=y(:)';
tmag=abs(tcurve(min(find(abs(tcurve)==max(abs(tcurve)))))/100);
tcurve=tcurve./(tmag+(tmag==0));
tcmin=min(tcurve);
betat2lb=[x(2)/2 1 0 x(2) 1 0];
betat2ub=[x(8) 50 Inf x(8) 50 Inf];
beta0=[x(min(find(tcurve==max(tcurve))))-x(3) ...
       x(2) ...
       max(tcurve) ...
       x(min(find(tcurve==max(tcurve)))) ...
       x(2) ...
       max(tcurve)-tcmin];
if beta0(1)<=betat2lb(1),
   beta0(1)=betat2lb(1);
elseif beta0(1)>=x(5),
   beta0(1)=x(5);
end
beta=lsqcurvefit('expdecaydiff',...
                 beta0,x,tcurve,betat2lb,betat2ub,opt);
beta0([3 6])=beta0([3 6]).*tmag;
beta([3 6])=beta([3 6]).*tmag;
fprintf(['%s time0: t1=%.1f tau1=%.1f a1=%.3f',...
         ' t2=%.1f tau2=%.1f a2=%.3f\n'],...
        [batstr],beta0(1),beta0(2),beta0(3),...
        beta0(4),beta0(5),beta0(6));
fprintf(['%s time : t1=%.1f tau1=%.1f a1=%.3f',...
         ' t2=%.1f tau2=%.1f a2=%.3f\n'],...
        [batstr],beta(1),beta(2),beta(3)*100,...
        beta(4),beta(5),beta(6)*100);

beta=beta(:);
beta0=beta0(:);

return

% display fit
figure(1)
clf
plot(x,y);
hold on
plot(x,expdecaydiff(beta,x),'r--');
hold off
