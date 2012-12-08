function x=gnoise(m,s,N);

dimcount=length(m);

x=randn(N,dimcount);
x=x*s^(1/2);

for dimidx=1:dimcount,
   x(:,dimidx)=x(:,dimidx)+m(dimidx);
end

