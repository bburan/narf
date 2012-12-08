% function f=gabor1d(beta,x)
%
% f=exp(-(x-m).^2/(2*s^2)) .* a/(sqrt(2*pi*s^2))
%
% where beta=[m s a]
%
% ie, 1d gabor function with mean m, std s, integral a, evaluated
% over the values in x
%
function f=gabor1d(beta,x)

m=beta(1);
s=beta(2);
a=beta(3);

f=exp(-(x-m).^2./(2.*s^2))./(sqrt(2*pi)*s);
f=f.*a;



