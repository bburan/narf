% function res=kvastimsim(cellid,batch)
%
% load results from kernfile in sRunData and display attentional
% modulation info ... kernfile generated from kerncomp4
%
% created SVD 10/18/02 - hacked from kerncomp4res.m
%
function res=kvastimsim(runidx,batch)

dbopen;
if isnumeric(runidx),
   sql=['SELECT * from sRunData WHERE id=',num2str(runidx)];
   rundata=mysql(sql);
   batchcount=length(rundata);
   cellid=rundata.cellid;
else
   cellid=runidx;
   goodbatch=zeros(1,length(batch));
   batchcount=0;
   for ii=1:length(batch),
      sql=['SELECT * from sRunData WHERE cellid="',cellid,'"',...
           ' AND batch=',num2str(batch(ii))];
      trd=mysql(sql);
      if ~isempty(trd),
         goodbatch(ii)=1;
         batchcount=batchcount+1;
         rundata(batchcount)=trd;
      end
   end
   batch=batch(find(goodbatch));
end

if batchcount==0,
   disp('no entries found in db!');
   if nargout>0,
      r=0;
   end
   return
end

global GCOLORMAP
GCOLORMAP=redblue;

rcsetstrings;

if length(rundata)>1,
   disp('cellres-kerncomp: only displaying first batch');
end

fprintf('%s.m: loading %s\n',mfilename,...
        [rundata(1).respath,rundata(1).kernfile,'.gz']);
zload([rundata(1).respath,rundata(1).kernfile,'.gz']);

% show spatial kernels for each attention condition and jackknife

kernfmt=params.kernfmt;
if strcmp(kernfmt,'pfftgr'),
   kernfmt='pfft';
   iconside=[sqrt(spacecount*2) sqrt(spacecount*2)];
end

H=zeros([size(vstrf(1).h,1),params.bootcount,attcount]);
mS=zeros([size(vstrf(1).mS,1),params.bootcount,attcount]);
H(:,:,1)=cat(2,vstrf(1,1,1,:).h);
mS(:,:,1)=cat(2,vstrf(1,1,1,:).mS);

if isfield(params,'oddout') & params.oddout,
   for attidx=1:attcount-1,
      H(:,:,attidx+1)=cat(2,vstrf(4+attidx,attidx,1,:).h);
      mS(:,:,attidx+1)=cat(2,vstrf(4+attidx,attidx,1,:).mS);
   end
else
   for attidx=2:attcount,
      H(:,:,attidx)=cat(2,vstrf(end,attidx-1,1,:).h);
      mS(:,:,attidx)=cat(2,vstrf(end,attidx-1,1,:).mS);
   end
end
stimcount=size(bstim,1);

mH=H(:,:,1);
kernsim=(bstim-repmat(mS(:,1)',[stimcount,1]))*mH;
kernpref=(mpatches(:,2:end)'-repmat(mS(:,1)',[size(mpatches,2)-1,1]))*mH;

mks=mean(kernsim,2);




respsim=zeros(paircount,params.bootcount);
for pidx=1:paircount,
   p1=pairidx(pidx,1);
   p2=pairidx(pidx,2);
   
   for bootidx=1:params.bootcount,
      if std(H(:,bootidx,p1+1))>0 & std(H(:,bootidx,p2+1))>0,
         respsim(pidx,bootidx)=xcov(H(:,bootidx,p1+1),...
                                    H(:,bootidx,p2+1),0,'coeff');
      end
   end
end

mrs=mean(respsim,2);
ers=std(respsim,0,2).*sqrt(params.bootcount);
rs=shrinkage(mrs,ers,0.5);


tp=(mpatches(:,2:end)-repmat(mS(:,1),[1 attcount-1]));
targsim=zeros(paircount,params.bootcount);
prefsim=zeros(paircount,params.bootcount);
dcsim=zeros(paircount,params.bootcount);

for pidx=1:paircount,
   p1=pairidx(pidx,1);
   p2=pairidx(pidx,2);
   
   targsim(pidx,:)=tp(:,p1)'*tp(:,p2) ./...
       (norm(tp(:,p1)).*norm(tp(:,p2)));
   prefsim(pidx,:)=1-abs(kernpref(p1,:)-kernpref(p2,:))./...
       abs(kernpref(p1,:)+kernpref(p2,:));
   
   dcsim(pidx,:)=abs(globaldc(p1,:)-globaldc(p2,:));
end

sum(~isnan(bresp(:,:)))
[pairidx targsim(:,1)]

[xcov(targsim(:),respsim(:),0,'coeff') ...
 xcov(prefsim(:),respsim(:),0,'coeff') ...
 xcov(targsim(:),dcsim(:),0,'coeff') ...
 xcov(prefsim(:),dcsim(:),0,'coeff') ...
]

figure(1);
clf
subplot(2,1,1);
scatter(targsim(:),respsim(:));

subplot(2,1,2);
scatter(prefsim(:),respsim(:));

res.ncount=sum(~isnan(bresp(:,:)));
res.prefsim=prefsim;
res.respsim=respsim;
res.targsim=targsim;
res.dcsim=dcsim;
res.targvresp=xcov(targsim(:),respsim(:),0,'coeff');
res.prefvresp=xcov(prefsim(:),respsim(:),0,'coeff');
res.targvdc=xcov(targsim(:),dcsim(:),0,'coeff');
res.prefvdc=xcov(prefsim(:),dcsim(:),0,'coeff');
res.predxc=predxc;


return
keyboard





mH=mean(H(:,:,1),2);
kernsim=(bstim-repmat(mS(:,1)',[stimcount,1]))*mH;
kernpref=(mpatches(:,2:end)'-repmat(mS(:,1)',[size(mpatches,2)-1,1]))*mH;





bincount=10;
[crap,sortidx]=sort(kernsim);
ll=round(linspace(1,stimcount+1,bincount+1));
kedges=[crap(ll(1:end-1)); inf];

rbinned=zeros(bincount,attcount-1);
for attidx=1:attcount-1,
   aokidx=find(~isnan(bresp(:,attidx+1)));
   
   for binidx=1:bincount,
      iidx=aokidx(find(kernsim(aokidx) >=kedges(binidx) & ...
                       kernsim(aokidx) < kedges(binidx+1)));
      if length(iidx)==0,
         iidx=(find(kernsim >=kedges(binidx) & ...
                    kernsim < kedges(binidx+1)));
      end
      rbinned(binidx,attidx)=nanmean(bresp(iidx,1));
   end
end

tp=(mpatches(:,2:end)-repmat(mS(:,1),[1 attcount-1]));
targsim=zeros(paircount,1);
prefsim=zeros(paircount,1);
respsim=zeros(paircount,1);
respdiff=zeros(paircount,1);

for pidx=1:paircount,
   p1=pairidx(pidx,1);
   p2=pairidx(pidx,2);
   
   targsim(pidx)=tp(:,p1)'*tp(:,p2) ./...
       (norm(tp(:,p1)).*norm(tp(:,p2)));
   prefsim(pidx)=1-abs(kernpref(p1)-kernpref(p2))./...
       abs(kernpref(p1)+kernpref(p2));
   
   respsim(pidx)=xcov(...
      (rbinned(:,p1)-mean(rbinned,2)),...
      (rbinned(:,p2)-mean(rbinned,2)),...
      0,'coeff');
   %respsim(pidx)=xcov(...
   %   gsmooth(rbinned(:,p1)-mean(rbinned,2),bincount/10),...
   %   gsmooth(rbinned(:,p2)-mean(rbinned,2),bincount/10),...
   %   0,'coeff');
   respdiff(pidx)=sum(abs(rbinned(:,p1)-rbinned(:,p2)));
   
end

plot(gsmooth(rbinned,bincount/10),'LineWidth',2)
legend('1','2','3','4')

sum(~isnan(bresp(:,:)))

[pairidx targsim prefsim respsim]

[xcov(targsim,respsim,0,'coeff') xcov(prefsim,respsim,0,'coeff')]

return

keyboard



return



% figure out number of valid fixations in each attentional condition
ncount=squeeze(sum(~isnan(bresp(:,end,:)),1));
cumncount=[0;cumsum(ncount(2:end))];
anyokidx=find(sum(~isnan(bresp(:,end,2:end)),3));

anycount=length(anyokidx);
stimcount=size(bstim,1);
spacecount=size(bstim,2);

tstim=bstim;
mS=mean(tstim,1);
tstim=tstim-repmat(mS,[stimcount 1]);
tsp=sqrt(sqrt(sum(tstim.^2,1)));
tstim=tstim./repmat(tsp,[stimcount 1]);
tsp=sqrt(sum(tstim.^2,2));
tstim=tstim./repmat(tsp,[1,spacecount]);

realattidx=zeros(size(bresp,1),1);
for attidx=2:attcount,
   tokidx=find(~isnan(bresp(:,1,attidx)));
   realattidx(tokidx)=attidx-1;
end

for ii=1:anycount,
   stimidx=anyokidx(ii);
   
   simvec=tstim(stimidx,:)*tstim(anyokidx,:)';
   
   NN=5;
   attcount=4;
   mm=zeros(1,attcount+1);
   rr=zeros(1,attcount+1);
   ss=zeros(1,attcount+1);
   mm(1)=realattidx(stimidx);
   rr(1)=bresp(stimidx,1).*1000;
   
   for attidx=1:4,
      
      ts=simvec(find(realattidx(anyokidx)==attidx &...
                     anyokidx~=stimidx));
      tr=bresp(find(realattidx(anyokidx)==attidx &...
                    anyokidx~=stimidx),1);
      
      [crap,simidx]=sort(-ts);
      
      mm(attidx+1)=mean(ts(simidx(1:NN))).*100;
      rr(attidx+1)=mean(tr(simidx(1:NN))).*1000;
      ss(attidx+1)=std(tr(simidx(1:NN))).*1000./sqrt(NN);
   end
   
   [mm; rr;ss]
   
   %[realattidx(simidx(1:NN)) -crap(1:NN)'.*100 bresp(simidx(1:NN),1).*1000]
   
   %showkern(tstim(simidx(1:10),:)','pfft');
   
   pause
end
