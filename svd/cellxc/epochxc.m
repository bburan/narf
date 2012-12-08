% function [epochxc,epocherr,ract,rpred]=epochxc(rpredpsth,resp,stim);
function [epochxc,epocherr,ract,rpred]=epochxc(rpredpsth,resp,stim);

% choose threshold automatically
[firstbins,lastbins]=movgetepochs(stim);
% [firstbins,lastbins]=movgetepochs(stim,60);

binok=zeros(length(firstbins),1);
for ii=1:length(firstbins),
   %if length(find(resp(firstbins(ii):lastbins(ii))==-1))==0 ...
   %       & lastbins(ii)-firstbins(ii) > 5,
   if length(find(~isnan(resp(firstbins(ii):lastbins(ii)))))>6,
      binok(ii)=1;
   end
end
firstbins=firstbins(find(binok))+2;
lastbins=lastbins(find(binok));
fprintf('Found %d fixations.\n',length(firstbins));

if 0,
   figure(gcf);
   clf
   plot(mean(stim));
   hold on
   plot(firstbins,300,'bo')
   plot(lastbins,300,'ro') 
   hold off
end

if length(firstbins)<=2,
   %keyboard
   ract=0;
   rpred=zeros(1,size(rpredpsth,2));
   epochxc=zeros(size(rpredpsth,2),1);
   epocherr=zeros(size(rpredpsth,2),1);
   return
end

ract=zeros(length(firstbins),1);
rpred=zeros(length(firstbins),size(rpredpsth,2));
epochxc=zeros(size(rpredpsth,2),1);
epocherr=zeros(size(rpredpsth,2),1);

%size(epochxc)

for ii=1:length(firstbins),
   
   % sum over non -1 bins, get index for each epoch for each pred
   % column  and for resp
   trespseg=resp(firstbins(ii):lastbins(ii));
   tpredseg=rpredpsth(firstbins(ii):lastbins(ii),:);
   ract(ii)=sum(trespseg(find(not(isnan(trespseg(:,1))))));
   rpred(ii,:)=sum(tpredseg(find(not(isnan(tpredseg(:,1)))),:),1);
   %ract(ii)=mean(trespseg(find(not(isnan(tpredseg(:,1))))));
   %rpred(ii,:)=mean(tpredseg(find(not(isnan(tpredseg(:,1)))),:),1);
end
if max(ract)>0,
   ract=ract./(max(ract));
end
rpredmax=max(rpred,[],1);

for jj=1:size(rpredpsth,2),
   if max(rpred(:,jj))>0,
      rpred(:,jj)=rpred(:,jj)./rpredmax(jj);
   end
   if sum(abs(rpred(:,jj)-mean(rpred(:,jj))))>0 & sum(abs(ract-mean(ract)))>0,
      %epochxc(jj)=xcov(ract,rpred(:,jj),0,'coeff');
      [epochxc(jj),epocherr(jj),tt,p]=randxcov(ract,rpred(:,jj),0,500);
   end
end




