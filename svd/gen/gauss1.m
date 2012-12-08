% function f=gauss1(beta,x,bcircular)
%
% f=a/(sqrt(2*pi*s^2)) .* exp(-(x-m).^2/(2*s^2)) + d
%
% where beta=[m s a d]
%
% ie, 1D gaussian with mean m, std s, peak height a (above d),
% offset d, evaluated over the values in x
%
function f=gauss1(beta,x,bcircular)

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
if s==0,
   s=1;
end

if ~exist('bcircular','var'),
   bcircular=0;
end
if bcircular,
   range=max(x)-min(x);
   step=mean(diff(x));
   if m<min(x) | m>max(x),
      m=mod(m-min(x),(range+step))+min(x);
   end
   locount=length(find(x<m));
   hicount=length(find(x>m));
   if locount>hicount,
      x(1:round((locount-hicount)/2))=x(1:round((locount-hicount)/2)) ...
          + range + step;
   elseif hicount>locount,
      xlen=length(x);
      
      x((xlen-round((hicount-locount)/2)+1):end)=...
          x((xlen-round((hicount-locount)/2)+1):end)-range-step;
   end
end

%f=exp(-(x-m).^2./(2.*s^2))./(sqrt(2*pi)*s);
%f=f./max(f).*a + d;
f=exp(-(x-m).^2./(2.*s^2));
f=f.*a + d;


