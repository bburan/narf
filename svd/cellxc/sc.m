
% set parms to get xccore to run correctly
params.maxlag=[-2 2];
params.maxlag=[0 0];
framecount=diff(params.maxlag)+1;
resp0idx=1;
params.resampfmt=1;
params.resampcount=10;
spacecount=size(stim,2);
firstseg=1;
xccore;

[H,neigs]=normalize(SR,sSA2,[],10.^(linspace(-1,-params.sfsstep,params.sfscount)));

thf=mean(H,5);
rpred=zeros(resplen,framecount,params.sfscount);
vpred=zeros(vallen,framecount,params.sfscount,nlcount);

predpower=zeros(params.sfscount,1);
predpower2=zeros(params.sfscount,params.resampcount);
vpredpower2=zeros(params.sfscount,params.resampcount);

for ii=1:params.sfscount,
   fprintf('%d ',ii);
   for resampidx=1:params.resampcount
      
      resamplen=rendidx(resampidx)-rstartidx(resampidx)+1;
      
      stim0=stim(rstartidx(resampidx):rendidx(resampidx),:)-repmat(mSall',[resamplen 1]);
      vstim0=cdata.stim-repmat(mSall',[size(cdata.stim,1) 1]);
      
      rpred(rstartidx(resampidx):rendidx(resampidx),:,ii)=stim0*H(:,:,ii,:,resampidx);
      vpred(:,:,ii)=cstim0*H(:,:,ii,:,resampidx);
      
      nnidx0=find(~isnan(aresp(rstartidx(resampidx):rendidx(resampidx))));
      predpower2(ii,resampidx)=std(rpred(nnidx0,resp0idx,ii)-resp(nnidx0));
      
      vpredpower2(ii,resampidx)=std(vpred(vnnidx,resp0idx,ii)-cdata.resp(vnnidx));
      
   end
   
   for respidx=1:framecount,
      ridxdiff=resp0idx-respidx;
      
      if respidx<resp0idx,
         r1=[rpred((ridxdiff+1):end,respidx,ii); ones(ridxdiff,1).*nan];
      else
         r1=[ones(-ridxdiff,1).*nan; rpred(1:(end+ridxdiff),respidx,ii)];
      end
      rpred(:,respidx,ii)=r1;
      
      if respidx<resp0idx,
         r1=[vpred((ridxdiff+1):end,respidx,ii); ones(ridxdiff,1).*nan];
      else
         r1=[ones(-ridxdiff,1).*nan; vpred(1:(end+ridxdiff),respidx,ii)];
      end
      vpred(:,respidx,ii)=r1;
   end
   
   predpower(ii)=std(rpred(nnidx,resp0idx,ii)-resp(nnidx));
   
end
disp('');

figure(1);
clf

plot(neigs,predpower);
hold on
plot(neigs,ones(size(neigs)).*std(aresp(nnidx)),'k');
plot(neigs,predpower2,'r');
hold off

figure(2);
clf
%showidx1=max(find(predpower<std(aresp(nnidx))));
%showidx2=max(find(sum(predpower2<std(aresp(nnidx)),2)));
showidx1=min(find(predpower==min(predpower)));
showidx2=min(find(mean(predpower2,2)==min(mean(predpower2,2))));

showkern(H(:,:,[showidx1 showidx2]),params.kernfmt);


sr2=zeros(128,3);

for aa=1:3
   g=find(~isnan(resp(:,aa)));
   ms=mean(stim(g,:),1);
   sr2(:,aa)=(stim(g,:)-repmat(ms,[length(g) 1]))'*resp(g,aa)./length(g);
end
