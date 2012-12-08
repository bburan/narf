function [h,nloutparms]=nlfitkern(stim,resp,h0,nl0);

gidx=find(~isnan(resp(:,1)));
resp=resp(gidx,:);
stim=stim(gidx,:);

r0=stim*h0;
r0=r0-nl0;
r0(find(r0<0))=0;
r0=r0-mean(r0);

r1=resp-mean(resp);
d1=sum(r0.^2);
if d1>0,
   scf=sum(r0.*r1)./d1;
else
   scf=1;
end

h0=h0*scf;

beta0=[h0;nl0];
nlopts=[];
nlopts=optimset('Display','iter');

fprintf('fitting to %d stim/resp samples...\n',length(gidx));

beta=fminsearch('nlkernerr',beta0,nlopts,stim,resp);


h=beta(1:length(h0));
nloutparms=beta(length(h0)+1:end);
