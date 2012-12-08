% function y=expdecaydiff(beta,x)
% 
%
%
function y=expdecaydiff(beta,x)

lat1=beta(1);
tau1=beta(2);
A1=beta(3);
lat2=beta(4);
tau2=beta(5);
A2=beta(6);

y1=A1*(1-exp(-(x-lat1)/tau1)) .* (x>=lat1);
y2=A2*(1-exp(-(x-lat2)/tau2)) .* (x>=lat2);

y=y1-y2;


return

lat=beta(1);
tau=beta(2);
a=beta(3);


y = a * (x-lat) .* exp(-(x-lat)./tau) .* (x>=lat);



return

