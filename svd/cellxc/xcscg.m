function h=xcscg(s,r,h0);

addpath /auto/k1/david/code/toolbox/netlab
addpath /auto/k1/david/code/demo

options=zeros(1,18);
options(1)=0;
options(2)=1e-10;
options(3)=1e-10;
%options(4)=0.01;
options(14)=20;
%options(18)=((step_mag/(net.nwts))*(T))/(T_init);

alpha(1)=1;
alpha(2)= 1 ./ std([h0(:); -h0(:)]);
%alpha
h = scg('testerr', h0, options, 'testgrad', s, r, alpha);

return

for itidx=1:10,
   
   % ARD terms. adjust alpha to be 1/var(h)
   %alpha(1)=length(r) ./ sum((h*s'-r').^2);
   alpha(1)=1;
   alpha(2)= 1 ./ var([h(:); -h(:)]);
   
   [itidx alpha testerr(h,s,r,alpha)]
   
   % update h with new alpha
   h = scg('testerr', h, options, 'testgrad', s, r, alpha);
   
end

