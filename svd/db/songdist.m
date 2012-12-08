% function songmatch(svdcount);
%
% svdcount - by default, chooses svd cutoff that gives min
%            prediction error on prior ratings.
%
function songdist(svdcount);

addpath /auto/k1/david/code/toolbox/kmeans

SKIPSINGLES=1;

%
% load all ratings from celldb
%
dbopen;
sql=['SELECT DISTINCT songid FROM dEval,dMusic',...
     ' WHERE dEval.songid=dMusic.id',...
     ' AND not(dup) AND not(filemissing) ORDER BY songid'];
evaldata=mysql(sql);
songid=cat(1,evaldata.songid);
songcount=length(songid);

sql=['SELECT DISTINCT user FROM dEval order by user'];
userdata=mysql(sql);
usercount=length(userdata);

%
% fill up a song X user matrix of ratings. empty values are Nan
%
evalmtx=nan.*zeros(length(songid),usercount);
users={};
for ii=1:usercount,
   users{ii}=userdata(ii).user;
   sql=['SELECT DISTINCT songid,rating FROM dEval,dMusic',...
        ' WHERE dEval.user="',userdata(ii).user,'"',...
        ' AND dEval.songid=dMusic.id',...
        ' AND not(dup) AND not(filemissing)',...
        ' ORDER BY songid'];
   evaldata=mysql(sql);
   
   usongid=cat(1,evaldata.songid);
   urating=cat(1,evaldata.rating);
   
   ff=find(ismember(songid,usongid));
   evalmtx(ff,ii)=urating;
end

% normalize rating matrix
mm=nanmean(evalmtx);
ss=nanstd(evalmtx);
nn=sum(~isnan(evalmtx));

evalmtx=(evalmtx-repmat(mm,[songcount,1])) ./ ...
        repmat(ss,[songcount,1]);
%evalmtx(find(isnan(evalmtx)))=0;

% remove songs only rated by minsongs people. keep this?
minsongs=2;

singleids=find(sum(~isnan(evalmtx),2)<minsongs);
if SKIPSINGLES,
   ss=find(sum(~isnan(evalmtx),2)>=minsongs);
   tevalmtx=evalmtx(ss,:);
   tsongid=songid(singleids);
else
   tevalmtx=evalmtx;
   tsongid=songid;
end

fprintf('%d/%d songs rated by >=%d users, %d users participating.\n',...
        size(tevalmtx,1),songcount,minsongs,usercount);

% 
% compute covaraiance matrix and PCs
%
cc=zeros(usercount);
n=zeros(usercount);
for u1=1:usercount,
   for u2=u1:usercount,
      n(u1,u2)=sum(~isnan(tevalmtx(:,u1).*tevalmtx(:,u2)));
      n(u2,u2)=n(u1,u1);
      cc(u1,u2)=nanmean(tevalmtx(:,u1).*tevalmtx(:,u2));
      cc(u2,u1)=cc(u1,u2);
   end
end

cc(isnan(cc))=0;
[u,s,v]=svd(cc);

if ~exist('svdcount','var'),
   svdcount=5;
end

% fill in unranked songs with pcs

zmtx=evalmtx;
zmtx(find(isnan(zmtx)))=0;
emtx=zmtx*u(:,1:svdcount);
pmtx=emtx*u(:,1:svdcount)';
%pmtx=zmtx;

nc=15;

[z, c] = kmeans(pmtx',nc);

cmembers=histc(c,1:nc);
[cc,orderidx]=sort(-cmembers);



pcol={'kx-','ro-','bs-','g*-','c--','kx-','ro-','bs-','g*-','c--'};
for ii=1:min([nc length(pcol)]),
   
   if ii==1,
      figure(1);
      clf
      leglabel={};
   elseif ii==min([nc length(pcol)])/2+1;
      hold off
      xticks(1:usercount,users);
      legend(leglabel);
      
      figure(2);
      clf
      leglabel={};
   end
   
   tu=z(:,orderidx(ii));
   
   %tu=tu./norm(tu);
   plot(tu,pcol{ii});
   hold on;
   
   leglabel={leglabel{:},...
             ['c',num2str(ii),' (',num2str(cmembers(orderidx(ii))),')']};
   
   % find closest songs to center of this cluster
   tu=z(:,orderidx(ii));
   dd=std(zmtx'-repmat(z(:,orderidx(ii)),[1 songcount]));
   %dd(singleids)=100;
   [dd0,closeidx]=sort(dd);
   fprintf('cluster %d:\n',ii);
   for jj=1:10,
      sql=['SELECT * FROM dMusic where id=',num2str(songid(closeidx(jj)))];
      songdata=mysql(sql);
      
      fprintf('%2d. %.4f (%d):',jj,dd0(jj),...
              sum(~isnan(evalmtx(closeidx(jj),:))));
      fprintf(' %s -- %s\n',songdata.artist,songdata.title);
   end
end

hold off
xticks(1:usercount,users);
legend(leglabel);


keyboard
return



testidx=1;


matchmtx=mean((repmat(pmtx(testidx,:),length(songid),1)-pmtx).^2,2);

[xx,matchids]=sort(matchmtx);

nsongs=15;
nearsongids=songid(matchids(1:nsongs));

sql=['SELECT * FROM dMusic where id=',num2str(songid(testidx))];
songdata=mysql(sql);

fprintf('%s-%s nearest songs:\n',songdata.artist,songdata.title);

for ii=1:nsongs,
   
   sql=['SELECT * FROM dMusic where id=',num2str(nearsongids(ii))];
   songdata=mysql(sql);
   
   fprintf('%5.2f: %s - %s\n',matchmtx(matchids(ii)),...
           songdata.artist,songdata.title);
end

evalmtx([testidx; matchids(1:nsongs)],:);

if 0,
   nrand=100;
   
   disp('randomizing...');
   tcc=eye(usercount);
   reigs=zeros(usercount,nrand);
   for ii=1:nrand,
      tem=reshape(shuffle(evalmtx(:)),size(evalmtx));
      for u1=1:usercount,
         for u2=(u1+1):usercount,
            tcc(u1,u2)=nanmean(evalmtx(:,u1).*shuffle(evalmtx(:,u2)));
            tcc(u2,u1)=tcc(u1,u2);
         end
      end
      tcc(isnan(tcc))=0;
      [ru,rs,rv]=svd(tcc);
      reigs(:,ii)=diag(rs);
   end
   
   emean=mean(reigs,2);
   estd=std(reigs,0,2);
   emax=emean+estd;
   emin=emean-estd;
   %reigs=sort(reigs')';
   %emax=reigs(:,round(nrand*0.95));
   %emean=median(reigs,2);
   %emin=reigs(:,round(nrand*0.05));

else
   emean=ones(usercount,1);
   emax=emean;
   emin=emean;
end

ee=zeros(usercount,1);
testids=find(~isnan(evalmtx(:,useridx)));
for svdidx=1:usercount,
   zmtx=evalmtx;
   zmtx(find(isnan(zmtx)))=0;
   zmtx(:,useridx)=0;
   emtx=zmtx*u(:,1:svdidx);
   tmtx=emtx*u(:,1:svdidx)';
   if SKIPSINGLES,
      % discount predictions where only one person has rated the song
      tmtx(singleids,:)=tmtx(singleids,:).*0.75;
   end
   
   %ee(svdidx)=xcov(tmtx(testids,useridx),evalmtx(testids,useridx),0,'coeff');
   ee(svdidx)=mean((tmtx(testids,useridx)-evalmtx(testids,useridx)).^2);
end

if ~exist('svdcount','var') | svdcount>usercount | svdcount<1,
   svdcount=min(find(ee==min(ee)));
   %svdcount=min(find(ee==max(ee)));
end
fprintf('predicting with %d principal components\n',svdcount);

zmtx=evalmtx;
zmtx(find(isnan(zmtx)))=0;
emtx=zmtx*u(:,1:svdcount);
pmtx=emtx*u(:,1:svdcount)';
if SKIPSINGLES,
   % discount predictions where only one person has rated the song
   pmtx(singleids,:)=pmtx(singleids,:).*0.75;
end

%
% FIND AND LIST BEST/WORST PREDICTIONS
%

predids=find(isnan(evalmtx(:,useridx)));

predratings=[pmtx(predids,useridx) songid(predids)];
[predratings,refids]=sortrows(predratings);
refids=predids(refids);

nsongs=15;
for ii=[length(predratings):-1:(length(predratings)-nsongs+1) 1:nsongs],
   
   if ii==length(predratings),
      fprintf('\n%s %d predicted best:\n',users{useridx},nsongs);
   elseif ii==1,
      fprintf('\n%s %d predicted worst:\n',users{useridx},nsongs);
   end
   
   sql=['SELECT * FROM dMusic where id=',...
        num2str(predratings(ii,2))];
   songdata=mysql(sql);
   sql=['SELECT * FROM dMusic where id=',...
        num2str(predratings(ii,2))];
   songdata=mysql(sql);
   
   fprintf('%5.1f:',predratings(ii,1)*ss(useridx)+mm(useridx));
   for jj=1:min([svdcount 5]),
      fprintf('%3.0f',emtx(refids(ii),jj).*10);
   end
   
   fprintf(' %s -- %s\n',songdata.artist,songdata.title);
end

%
% PLOT SUMMARY STATS
%

figure(1);
clf

subplot(2,1,1);
pcol={'kx-','ro-','bs-','g*-','c--'};
leglabel={};
for ii=1:min([svdcount length(pcol)]),
   tu=u(:,ii);
   %if sum(tu)<0,
   %   tu=-tu;
   %end
   plot(tu,pcol{ii});
   hold on;
   
   leglabel{ii}=['pc ',num2str(ii)];
end
hold off

xticks(1:usercount,users);
legend(leglabel);

subplot(2,1,2);
ev=diag(s);
ev=ev./sum(ev).*100;

eescf=max(ev); %100
plot(ee*eescf,'r-');
hold on
plot(ev);
plot(svdcount,ee(svdcount)*eescf,'o');
%plot([1 usercount],[1 1]*eescf,'r--');
hold off

legend('pred cc','% eig variance');
xlabel('eigenvector rank');


emtx=zmtx*u(:,:);
goodknown=shuffle(find(evalmtx(:,useridx)>1));
goodknown=goodknown(1:min([100 length(goodknown)]));
badknown=shuffle(find(evalmtx(:,useridx)<-1));
badknown=badknown(1:min([100 length(badknown)]));
otherknown=shuffle(find(evalmtx(:,useridx)<=1 & evalmtx(:,useridx)>=-1));
otherknown=otherknown(1:min([100 length(otherknown)]));
goodpred=refids(length(predratings):-1:(length(predratings)-nsongs+1));
badpred=refids(1:nsongs);

dd=u(useridx,:);
dd=mean(abs(emtx(goodpred,:)),1);
[tt drank]=sort(-abs(dd(1:min([5 svdcount]))));
d1=drank(1);
d2=drank(2);

figure(2);

clf
scatter(emtx(otherknown,d1),emtx(otherknown,d2),'k.');
hold on
scatter(emtx(goodknown,d1),emtx(goodknown,d2),'g.');
scatter(emtx(badknown,d1),emtx(badknown,d2),'r.');
scatter(emtx(goodpred,d1),emtx(goodpred,d2),'go');
scatter(emtx(badpred,d1),emtx(badpred,d2),'ro');

a=axis;
plot([a(1) a(2)],[0 0],'k--');
plot([0 0],[a(3) a(4)],'k--');

hold off
xlabel(sprintf('d1=%d',d1));
ylabel(sprintf('d2=%d',d2));
title('PC contributions to predictions');
