% f=polfp(beta,r,mS,bgpix,stimfilterparms)
function f=polfp(beta,r,mS,bgpix,stimfilterparms)

r=r(1);

z = pol_grat2(r,beta(1),beta(2),beta(3),beta(4),1);
if std(z(:))>0,
   z=(z-mean(z(:)))./std(z(:));
end
fspace=z.*bgpix(2)+bgpix(1);

f=(movpower(fspace,stimfilterparms{:})-mS).*beta(5);
