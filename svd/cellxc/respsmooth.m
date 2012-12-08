
cellid='m0118';

[cellfiledata,cellids,cellfileids]=dbgetscellfile('cellid',cellid,'runclassid',3,...
                                                  'respfmtcode',1,'fmt','downsamp');

pfth=[];
for ii=1:length(cellfiledata),
   fprintf('loading %s\n',cellfiledata(ii).respfile);
   r=load([cellfiledata(ii).path,cellfiledata(ii).respfile]);
   pfth=cat(1,pfth,r.pfth);
end
   

fixcount=size(pfth,1);
msperbin=6;
pfth=pfth(:,51:250);
lags=r.movformat.tstart:r.movformat.tstop;
lags=lags(51:250);



bincount=size(pfth,2);
pfth(:,bincount:(msperbin*ceil(bincount/msperbin)))=nan;
bincount=size(pfth,2);
bincount2=bincount/msperbin;
lags2=linspace(lags(1)+msperbin/2,lags(1)+bincount-msperbin/2,bincount2);

pfth2=squeeze(sum(reshape(pfth,fixcount,msperbin,bincount2),2));


%pfth(:,end)=nan;
%h=diff(find(pfth'==1));

mm=nanmean(pfth2);
pfth2=pfth2(:,find(~isnan(mm)));
lags2=lags2(find(~isnan(mm)));
mm=mm(find(~isnan(mm)));
for ii=1:size(pfth2,2),
   pfth2(find(isnan(pfth2(:,ii))),ii)=mm(ii);
   pfth2(:,ii)=pfth2(:,ii)-mm(ii);
   
end

sa=pfth2'*pfth2 ./ fixcount;
[u,s,v]=svd(sa);


figure(1);
clf
mm=mm./sqrt(sum(mm.^2));
for ii=1:25,
   subplot(5,5,ii);
   plot(lags2,mm,'r--');
   hold on
   plot(lags2,u(:,ii));
   hold off
   axis([lags2(1) lags2(end)  -1 1]);
   title(sprintf('%s eig %d',cellid,ii));
end

legend('mean','eig');
