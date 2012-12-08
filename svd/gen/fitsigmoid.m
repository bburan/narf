% function [sigparms,sigerr]=fitsigmoid(stim,resp,showfit[=0],forcevalues);
%
% stim is range of inputs.
% resp is range of outputs.
%
% sigparms=[m:meanx slope a1:minamp a2:maxamp];
%    fed to sigmoid_svd.m such that:
%       s=(a2-a1)./(exp(-(x-m)*slope)+1) + a1;
%
% created SVD 12/02
%  use sigmoid_svd rather than "sigmoid" to avoid conflict with nsltoolbox
%
function [sigparms]=fitsigmoid(stim,resp,showfit,forcevalues,njack);

if ~exist('showfit','var'),
   showfit=0;
end
% currently not used
if ~exist('njack','var') | njack<1,
   njack=1;
end

%beta0=[median(stim) std(resp)/std(stim) ...
%       min(resp)+std(resp(:)) 3.*std(resp(:))];
%lb=[min(stim)  0      min(resp)  min(resp)];
%ub=[max(stim)  inf    max(resp)  max(resp)];

%       std(resp)/std(stim)         ...
%       log(9)./(std(stim).*2)     ...


if mean(abs(stim))==0,
   sigparms=[0 0 0 0];
   return
end

beta0=[median(stim)-1.5*std(stim)  ...
       std(resp)/std(stim)         ...
       mean(resp)-1.5*std(resp(:)) ...
       3.*std(resp(:))                 ];
lb=[min(stim)     eps     -inf                   mean(resp) ];
ub=[max(stim)     inf     mean(resp)             inf         ];

if exist('forcevalues','var'),
   useidx=find(~isnan(forcevalues));
   beta0(useidx)=forcevalues(useidx);
   lb(useidx)=forcevalues(useidx)-eps.*(10.*forcevalues(useidx));
   ub(useidx)=forcevalues(useidx)+eps.*(10.*forcevalues(useidx));
end

if find(abs(ub-lb)==0),
   ub-lb;
   sigparms=[0 0 0 0];
   return
end

fitopt=optimset('Display','off');
sigparms=lsqcurvefit('sigmoid_svd',beta0,stim,resp,lb,ub,fitopt);

% plot fit over data
if showfit,
   hold off
   ptshow=min([length(stim) 500]);
   
   scatter(stim(1:ptshow),resp(1:ptshow),'.');
   xx=linspace(min(stim),max(stim),50);
   hold on
   plot(xx,sigmoid_svd(beta0,xx),'r--');
   plot(xx,sigmoid_svd(sigparms,xx));
   hold off
   drawnow
   %title(sprintf('sigmoid fit: m=%.2f s=%.2f a=%.2f d=%.2f',sigparms));
   title(sprintf('sigmoid fit: x10=%.2f s=%.2f a=%.2f d=%.2f',sigparms));
   legend('init','fit');
end
