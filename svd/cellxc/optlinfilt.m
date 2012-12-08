

respfile='/auto/k3/david/data/R110A/R110A.gridreview.androsgrid.1.all.d-.psth.neg_flag';

r=respload(respfile);

p=r(:,1);
r=compact_raster_matrix2(r(:,2:end),1);
p=p(size(p,1)-size(r,1)+1:end);
trialcount=size(r,2);

p=[ones(20,1)*nan; p];
r=[zeros(20,trialcount); r];



p=repmat(p,[trialcount 1]);
r=r(:);
minlag=-40;
maxlag=40;

[SR,n,mS,mR,tSA,sSA2]=movxc(r,p,[minlag maxlag]);
[SR,mS,mR,tSA,sSA2]=normalize0(SR,n,mS,mR,tSA,sSA2,1);
[H,nfactors]=normalizetime(SR,tSA);

tt=minlag:maxlag;

figure(1);
clf
plot(tt,H(:,:,5));
