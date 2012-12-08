% function f=gaussfp(beta,x)
%
% f=a/(sqrt(2*pi*s^2)) .* exp(-(x-m).^2/(2*s^2)) + d
%
% where beta=[or sf orw sfw a d]
%
% ie, 2D gaussian with mean (or,sf), peak height a (above d),
% offset d, evaluated over the values in x
%
% CURRENTLY FORCING d=0
%
function f=gaussfpN(beta,x)

%f=exp(-(x-m).^2./(2.*s^2))./(sqrt(2*pi)*s);
%f=f./max(f).*a + d;

N=length(beta)/5;

if ~exist('x','var'),
   [xx,yy]=meshgrid(-8:7,-8:7);
   x=[xx(:) yy(:)];
elseif length(x)==1,
   [xx,yy]=meshgrid((-x/2):(x/2-1),(-x/2):(x/2-1));
   x=[xx(:) yy(:)];
end

f=zeros(length(x),1);
for nn=1:N,
   tbeta=beta((nn-1)*5+(1:5));
   f=f+gaussfp(tbeta,x);
end

