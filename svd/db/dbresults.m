% function dbresults(batchid,pccount,latidx)

function dbresults(batchid,pccount,latidx)

dbopen;

if ~exist('batchid'),
   batchid=11
end

sql=['SELECT * FROM sResults WHERE batch=',num2str(batchid)];
resdata=mysql(sql);

for ii=1:length(resdata),
   eval(char(resdata(ii).matstr));
end

okidx=zeros(length(preddata),1);
for ii=1:length(preddata),
   if ~isempty(preddata(ii).p),
      okidx(ii)=1;
   end
end
ttt=find(okidx)';

if 0,
   keyboard
end
if 0,
   for ii=ttt,
      size(preddata(ii).p)
   end
end


cellid={preddata(ttt).cellid};
cellcount=length(cellid);
cellcount=cellcount+1;
cellid{cellcount}='mean';

runidx=find(okidx);
runidx(cellcount)=0;

predxc=cat(5,preddata(ttt).predxc);
tpredxc=predxc;
tpredxc(find(isnan(tpredxc)))=0;
predxc=cat(5,predxc,mean(tpredxc,5));
%attsimxc=cat(4,preddata(ttt).attsimxc);
%attsimxc=cat(4,attsimxc,mean(attsimxc,4));

TSON=0;
if TSON,
   targsim=cat(4,preddata(ttt).targsim);
   targsim=cat(4,targsim,mean(targsim,4));
end

% make the last entry in each summary statistic the mean across
% cells
if length(preddata(ttt(1)).p)>1,
   p=repmat(zeros(size(preddata(ttt(1)).p)),[1 1 1 length(ttt)]);
   for ii=1:length(ttt),
      pcount0=min([size(preddata(ttt(ii)).p,1) size(p,1)]);
      pcount=min([size(preddata(ttt(ii)).p,2) size(p,2)]);
      pcount2=min([size(preddata(ttt(ii)).p,3) size(p,3)]);
      p(1:pcount0,1:pcount,1:pcount2,ii)=...
          preddata(ttt(ii)).p(1:pcount0,1:pcount,1:pcount2);
   end
   p=cat(4,p,mean(p,4));
else
   p=zeros(size(predxc));
end

% display prediction scores for this latency bin (4 means 0-100 ms after
% fixation onset)
latcount=size(predxc,3);
latbase=min([latcount 1])-1;
if ~exist('latidx'),
   latidx=min([latcount 6]);
end
PCLOW=10;
if ~exist('pccount'),
   pccount=18
end
nlidx=1;

% PTHRESH=0.012; % p < 1-(1-0.05)^(1/4)
PTHRESH=0.008; % p < 1-(1-0.05)^(1/6
%PTHRESH=0.05
figure(1);
clf
minx=min(min(min(predxc(:,:,latidx,:))));
maxx=max(max(max(predxc(:,:,latidx,:))));
if minx<-maxx & maxx>0,
   minx=-maxx;
end
for ii=1:cellcount,
   subplot(ceil(cellcount/4),4,ii);
   imagesc(predxc(:,:,latidx,nlidx,ii),[minx maxx]);
   axis image
   axis off
   title(sprintf('%s (%d)',cellid{ii},runidx(ii)));
end
colorbar
colormap(hot);

attcount=size(predxc,1);
attpcount=size(p,1);

disp('cell  SMT (pred in/out) (att comp 1-2 1-3 1-4 2-3 2-4 3-4)');
pxc=zeros(cellcount,attcount);
pxcout=zeros(cellcount,attcount);
pxcatt=[];
sigcount=[0 0];
attpairs=[1 2; 1 3; 1 4; 2 3; 2 4; 3 4];

xcallin=[];
xcallout=[];
pall=[];

attbetter=zeros(cellcount-1,attcount-1);
for ii=1:cellcount,
   fprintf('%-5s',cellid{ii});
   for attidx=1:attcount,
      %tpredxc=predxc(:,:,9,nlidx,ii);
      tpredxc=predxc(:,:,latidx,nlidx,ii);
      %tpredxc(find(isnan(tpredxc)))=0;
      pxc(ii,attidx)=tpredxc(attidx,attidx);
      pxcout(ii,attidx)=mean(tpredxc(attidx,[2:attidx-1 (attidx+1):end]));
      
      if attidx>1 & ii<cellcount,
         xcallin=[xcallin tpredxc(attidx,attidx).*ones(1,attcount-2)];
         xcallout=[xcallout tpredxc(attidx,[2:attidx-1 (attidx+1):end])];
      end
      
      pxcout(ii,attidx)=mean(tpredxc(attidx,[2:attidx-1 (attidx+1):end]));
      fprintf(' %5.2f/%5.2f',pxc(ii,attidx),pxcout(ii,attidx));
      
      if ii<cellcount & attidx>1 & pxc(ii,attidx)>pxcout(ii,attidx),
         attbetter(ii,attidx-1)=1;
      end
   end
   
   % print out in vs out pred xcs. transfer over to pairwise
   % comparison if it looks like it's gonna help.
   
   if preddata(ttt(1)).p~=0,
      for attidx=1:attpcount,
         pall=[pall; p(attidx,pccount,latidx-latbase,ii)];
         if p(attidx,pccount,latidx-latbase,ii)<=PTHRESH,
            sextra='*';
            %pxcatt[pxcatt; pxc(ii,attidx) pxcout(ii,attidx)];
            pxcatt=[pxcatt; pxc(ii,attpairs(attidx,1)+1) ...
                    pxc(ii,attpairs(attidx,2)+1)];
         else
            sextra=' ';
         end
         fprintf('%s%5.3f',sextra,p(attidx,pccount,latidx-latbase,ii));
      end
      
      % if qualifies, mark as cell with significant attentional effect
      if min(p(1:attpcount,PCLOW,latidx-latbase,ii))<=PTHRESH,
         fprintf('*%d',PCLOW);
         sigcount(2)=sigcount(2)+1;
      end
      if min(p(1:attpcount,pccount,latidx-latbase,ii))<=PTHRESH,
         fprintf('*%d',pccount);
         sigcount(1)=sigcount(1)+1;
      end
   else
      pxcatt=[0 0];
   end
   fprintf('\n');
end
fprintf(['sig att: %d/%d %.2f%% (%dPC) %d/%d %.2f%% (%dPC)',...
         ' p<%.3f respidx=%d\n'],...
        sigcount(2),cellcount-1,sigcount(2)./(cellcount-1).*100,...
        PCLOW,...
        sigcount(1),cellcount-1,sigcount(1)./(cellcount-1).*100,...
        pccount,PTHRESH,latidx);
xcallin=xcallin(find(~isnan(xcallout)));
xcallout=xcallout(find(~isnan(xcallout)));
[pp,mm,ss]=randpairtest(xcallin,xcallout,2000);

fprintf('mean xc in/out: %.3f/%.3f (p<%.3f): %d/%d (%.1f%%) attstates\n',...
        mean(xcallin),mean(xcallout),pp,...
        sum(xcallin>xcallout),length(xcallout),...
        sum(xcallin>xcallout)/length(xcallout)*100);

figure(2);
clf
subplot(2,3,1);
scatter(xcallout(:),xcallin(:));
hold on
plot([0 1.0],[0 1.0],'k--');
plot(nanmean(xcallout),nanmean(xcallin),'r+','LineWidth',2);
hold off
axis equal
xlabel('pcx out');
ylabel('pcx in');
title('attention preds in vs. out');

subplot(2,3,2);
scatter(mean(pxcout(:,2:end),2),pxc(:,1));
hold on
plot([0 1.0],[0 1.0],'k--');
plot(median(mean(pxcout(:,2:end),2)),median(pxc(:,1)),...
     'r+','LineWidth',2);
hold off
axis equal
xlabel('pcx att in');
ylabel('pcx all att');
title('data lim: all vs. atten preds');

subplot(2,3,3);
mxcin=zeros(latcount,1);
mxcout=zeros(latcount,1);
for attidx=2:attcount,
   mxcin=mxcin+squeeze(predxc(attidx,attidx,:,nlidx,end)) ./ (attcount-1);
   mxcout=mxcout+squeeze(sum(predxc(attidx,[2:attidx-1 (attidx+1):end],...
                                     :,nlidx,end),2)) ./ ...
          ((attcount-2)*(attcount-1));
end
plot([mxcout mxcin]);
legend('out','in');
%scatter(pxcatt(:,2),pxcatt(:,1));
%hold on
%plot([0 1.0],[0 1.0],'k--');
%plot(median(pxcatt(:,2)),median(pxcatt(:,1)),...
%     'r+','LineWidth',2);
%hold off
%axis equal
%xlabel('xc att out');
%ylabel('xc att in');
title(sprintf('att mod preds (latidx=%d,pccount=%d)',latidx,pccount));

subplot(2,3,4);
%plot(squeeze(sum((min(p(2:end,:,:,:),[],1)<=PTHRESH),4)),'LineWidth',2);
plot(squeeze(sum((min(p(1:end,:,:,:),[],1)<=PTHRESH),4)),'LineWidth',2);
hold on
%plot(squeeze(sum((min(p(2:end,:,latidx-latbase,:),[],1)<=PTHRESH),4)),...
%     'kx','LineWidth',2);
plot(squeeze(sum((min(p(1:end,:,latidx-latbase,:),[],1)<=PTHRESH),4)),...
     'kx','LineWidth',2);
hold off
xlabel('pc count');
ylabel('N signif atten');
legend('lat1','lat2');
title('sig att cells vs. pc count');

subplot(2,3,5);
hist(pall,15);
title('pair p values')

subplot(2,3,6);
hist(xcallin(:)-xcallout(:),15);
title('pair pred diff (in - out)');

if TSON,
   subplot(2,3,5);
   ts1=targsim(:,pccount,1,:);
   ts2=targsim(:,pccount,2,:);
   tp=p(:,pccount,latidx-latbase,:);
   sigsim=[ts1(:) ts2(:) tp(:)];
   scatter(sigsim(:,1),sigsim(:,3));
   title('att diff vs targ diff (unnorm)');
   xlabel('targ1.*targ2');
   ylabel('att p');
   
   subplot(2,3,6);
   scatter(sigsim(find(tp<PTHRESH),2),sigsim(find(tp<PTHRESH),3));
   title('att diff vs targ diff (norm)');
   xlabel('angle(targ1,targ2)');
   ylabel('att p');
end

%keyboard
return

figure(3);
clf
atxc=squeeze(attsimxc(:,:,pccount,:));
minx=min(atxc(:));
maxx=max(atxc(:));
if minx<-maxx & maxx>0,
   minx=-maxx;
end
for ii=1:cellcount-1,
   subplot(ceil(cellcount/4),4,ii);
   imagesc(atxc(:,:,ii),[minx maxx]);
   axis image
   axis off
   title(sprintf('atxc %s (%d)',cellid{ii},runidx(ii)));
end
colorbar
colormap(hot);



disp('cell  T   (pred in/out)');
for ii=1:cellcount,
   fprintf('%-5s',cellid{ii});
   for attidx=1:attcount,
      pxc=predxc(attidx,attidx,2,ii);
      pxcout=mean(predxc(attidx,[2:attidx-1 (attidx+1):end],latidx,ii));
      fprintf(' %5.3f (%4.1f/%4.1f)',p(attidx,1,1,ii),pxc,pxcout);
   end
   
   % if qualifies, mark as cell with significant attentional effect
   if min(p(:,1,1,ii))<PTHRESH,
      fprintf('*');
   end
   fprintf('\n');
end


