

m=[0 0];
s=[1 0; 0 1];
N=100;


stim=gnoise(m,s,N);
h=[1 0]';
resp=stim*h;
resp(resp<0)=0;

h0=[0.5 0.5]';
nl0=[0.5];
%nl0=[0.5;1];

[h,nloutparms]=nlfitkern(stim,resp,h0,nl0);




h=zeros(spacecount,params.resampcount);
sfsmax=35;
for respidx=1:params.resampcount,
   [sU sS sV] = svd(sSA2(:,:,respidx));
   
   sS=diag(sS);
   sS=sS(1:sfsmax);
   
   eigstim=stim*sU(:,1:sfsmax);
   eigH0=sU(:,1:sfsmax)'*SR(:,:,respidx) ./ sS;
   nl0=0.1;
   
   resp0=resp;
   resp0(rstartidx(respidx):rendidx(respidx),:)=nan;
   
   [eigH,nloutparms]=nlfitkern(eigstim,resp0,eigH0,nl0);
   
   h(:,respidx)=sU(:,1:sfsmax)*eigH;
   
end

mH=mean(h,2);
eH=std(h,0,2).*sqrt(size(h,2));
hfinal=shrinkage(mH,eH,0.2);

