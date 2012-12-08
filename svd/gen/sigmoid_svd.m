% function s=sigmoid(parms,x);
%
% s=a./(exp(-(x-m)*slope)+1) + d;
%
% parms=[x10:10%ofmax slope d:minamp a:(maxamp-minamp)];
% m=log(9)./slope + x10;
%
% OLD BAD
% m=(x10+x90)/2;
% s=log(9)./ (0.5*x10 - 0.5*x90);
%
% x is range of inputs.
%
% example:
% x=1:100;
% y=sigmoid([28 0.1 0 1],x);
% plot(x,y);
%
% created SVD 12/02
% modified SVD 1/04 : in parms: made a difference, changed m to x10
%
function s=sigmoid(parms,x);

if sum(parms)==0,
   s=x;
   return
end

x10=parms(1);
slope=parms(2);
%x90=parms(2);
d=parms(3);
a=parms(4);

%m=(x10+x90)/2;
%slope=log(9)./ 0.5 ./ (x90-x10 + eps);  % force non-zero denom
m=log(9)./slope + x10;

s=a./(exp(-(x-m)*slope)+1) + d;





