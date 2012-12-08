

T=1000;
N=10;
h=[1 2 zeros(1,N-2)];

stim=randn(T,N);
stim(:,3)=stim(:,3)+stim(:,1)./10;
stim(:,2)=stim(:,2)*1.5;
stim(:,4)=stim(:,4)*30+stim(:,1);
resp=(stim*h');
resp=resp.*abs(resp);
resp1=resp;

htest=zeros(size(h));
for ii=1:10,
   r1=stim;
   d1=sum(stim.^2,1);
   gn=resp1'*r1./d1;
   maxh=find(abs(gn)==max(abs(gn)))
   rpred=gn(maxh).*stim(:,maxh);
   
   [maxh gn(maxh)]
   [norm(resp1) norm(rpred) norm(resp1-rpred)]
   resp1=resp1-rpred;
   
   htest(maxh)=htest(maxh)+gn(maxh);
end

