% function y=tadapt(parms,x)
%
% apply temporal adaptation to vector x using parameters y.
%
% created SVD 2006-01-15
%
function y=tadapt(parms,x)

adaptstr=parms(1);
tdecay=parms(2);
thresh=parms(3);

tpower=x(:)-thresh;
tpower=tpower.*(tpower>0)./max(tpower);

ff=[0 0 0 0 0 0 0 exp(-[0 1 2 3 4 5].*tdecay)]';
ff=ff./sum(ff);
tpower2=rconv2(tpower,ff);

% remove the effects of nans
tpower2(isnan(tpower2))=tpower(isnan(tpower2));

y=tpower.*tpower./((1-adaptstr).*tpower+adaptstr.*tpower2 + ...
                   (tpower+tpower2==0));

