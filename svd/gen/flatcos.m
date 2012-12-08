% function f=flatcos(beta,x,e)
%
% f=d + a cos (x + (60-s)*sin(x-m) - m)
%
% where beta=[m s a d]'
% 
% NOTE: angles should range from 0 to 180 deg!
% (minus 60 so that bw looks like reasonable bw)
% phi
% ie, flat-topped cos with mean m, width s, peak height a (above d),
% offset d, evaluated over the values in x
%
% SVD 10/01
%
function f=flatcos(beta,x,e)

if length(beta)<1,
   m=0;
else
   m=beta(1);
end
if length(beta)<2,
   s=1;
else
   s=beta(2);
end
if length(beta)<3,
   a=1;
else
   a=beta(3);
end
if length(beta)<4,
   d=0;
else
   d=beta(4);
end

% convert to radians if necessary
%if abs(max(x)-min(x))>3.2,
%   x=x./(max(x)-min(x)) .* 2*pi;
%end

% convert from half-circle deg to full-circle radians:
x=x./180*2*pi;
m=m./180*2*pi;
s=(60-s)./180*2*pi;

f=d + a .* cos(x + s.* sin(x-m) - m);

% divide by std err if passed as a parameter (tweak to get fitting
% to work)
if exist('e'),
   f=f./e;
end

