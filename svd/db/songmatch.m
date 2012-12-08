% function songmatch(svdcount,user);
%
% svdcount - by default, chooses svd cutoff that gives min
%            prediction error on prior ratings.
%
function songmatch(svdcount,user);

if ~exist('user','var');
   user=getenv('USER');
end

dbopen;

if 0,
   sql=['SELECT DISTINCT artist,title,file,songid FROM dEval,dMusic',...
        ' WHERE dEval.songid=dMusic.id',...
        ' AND not(dup) AND not(filemissing) ORDER BY songid'];
   evaldata=mysql(sql);
   songid=cat(1,evaldata.songid);
   songcount=length(songid);
   titles={};
   artists={};
   files={};
   for ii=1:length(evaldata),
      titles{ii}=evaldata(ii).title;
      artists{ii}=evaldata(ii).artist;
      files{ii}=evaldata(ii).file;
   end
   
   keyboard
   save /auto/k5/david/ratings evalmtx users titles artists files songid
end

% get rating matrix
minevals=8;
[evalmtx,users,songid]=songevalmtx(minevals);

if 0,
   disp('testing with top 100 only!');
   load /auto/k1/hayden/songids.mat
   dropidx=find(~ismember(songid,top100_songid));
   evalmtx(dropidx,useridx)=nan;
end


[cc,n,mm,ss]=songdbcorrmtx(evalmtx);

minsongs=2;
singids=find(sum(~isnan(evalmtx),2)>=minsongs);
evalmtx=evalmtx(singids,:);
songid=songid(singids);
songcount=length(songid);

evalmtx=(evalmtx-repmat(mm,[songcount,1])) ./ ...
        repmat(ss,[songcount,1]);

useridx=find(strcmp(user,users));
usercount=length(users);
if isempty(useridx),
   fprintf('user %s: no match\n',user);
   return
end

% find eigenvectors
[u,s,v]=svd(cc);

disp('finding optimal number of eigenvectors...');
ee=zeros(usercount,1);
testids=find(~isnan(evalmtx(:,useridx)));
for svdidx=1:usercount,
   zmtx=evalmtx;
   zmtx(find(isnan(zmtx)))=0;
   zmtx(:,useridx)=0;
   emtx=zmtx*u(:,1:svdidx);
   tmtx=emtx*u(:,1:svdidx)';
   
   ee(svdidx)=xcov(tmtx(testids,useridx),evalmtx(testids,useridx),0,'coeff');
end

if ~exist('svdcount','var') | svdcount>usercount | svdcount<1,
   svdcount=min(find(ee==max(ee)));
end
fprintf('predicting with %d eigenvectors\n',svdcount);

sd=diag(s);
sfrac=sum(sd(1:svdcount))./sum(sd);

fprintf('%d eigenvectors (%.1f%% of variance)\n',svdcount,sfrac*100);
uadj=u./sqrt(sfrac);

zmtx=evalmtx;
zmtx(find(isnan(zmtx)))=0;
emtx=zmtx*uadj(:,1:svdcount);
pmtx=emtx*uadj(:,1:svdcount)';

%
% FIND "MOST INFORMATIVE" SONGS
%
if 1,
   vv=nanstd(evalmtx')';
   ct=sum(~isnan(evalmtx),2);
   vv=vv.*sqrt(ct);
   [vv,topsongs]=sort(-vv);
   vv=-vv;
   disp('Top 10 most informative songs [=std*sqrt(ratingcount)]');
   for sidx=1:10,
      sql=['SELECT * FROM dMusic where id=',num2str(songid(topsongs(sidx)))];
      songdata=mysql(sql);
      fprintf('songid: %5d std=%.1f ct=%d %s - %s\n',...
              songid(topsongs(sidx)),vv(sidx),ct(topsongs(sidx)),...
              songdata.artist,songdata.title);
   end
   
else
   % OLD: "TYPICAL" songs for each eig
   for svdidx=1:min([10 svdcount]),
      [zz,tt]=sort(emtx(:,svdidx));
      for sidx=[1 2 length(songid)+(-1:0)],
         sql=['SELECT * FROM dMusic where id=',num2str(songid(tt(sidx)))];
         songdata=mysql(sql);
         fprintf('eig: %d proj: %.4f songid: %d %s - %s\n',...
                 svdidx,zz(sidx),songid(tt(sidx)),...
                 songdata.artist,songdata.title);
      end
   end
end

%
% FIND AND LIST BEST/WORST PREDICTIONS
%

predids=find(isnan(evalmtx(:,useridx)));

predratings=[pmtx(predids,useridx) songid(predids)];
[predratings,refids]=sortrows(predratings);
refids=predids(refids);

for uidx=1:usercount,
   fprintf('vs. %s: %.3f\n',users{uidx},cc(useridx,uidx));
end
meaneval=nanmean(evalmtx(:,[1:useridx-1 useridx+1:end])')';
fidx=find(~isnan(meaneval) & ~isnan(evalmtx(:,useridx)));
fprintf('vs. all: %.3f\n',...
        xcov(meaneval(fidx),evalmtx(fidx,useridx),0,'coeff'));

goodsongs=15;
badsongs=5;
for ii=[length(predratings):-1:(length(predratings)-goodsongs+1) 1:badsongs],
   
   if ii==length(predratings),
      fprintf('\n%s %d predicted best:\n',users{useridx},goodsongs);
   elseif ii==1,
      fprintf('\n%s %d predicted worst:\n',users{useridx},badsongs);
   end
   
   sql=['SELECT * FROM dMusic where id=',num2str(predratings(ii,2))];
   songdata=mysql(sql);
   %sql=['SELECT * FROM dMusic where id=',num2str(predratings(ii,2))];
   %songdata=mysql(sql);
   
   fprintf('%5.1f:',predratings(ii,1)*ss(useridx)+mm(useridx));
   %for jj=1:min([svdcount 5]),
   %   fprintf('%3.0f',emtx(refids(ii),jj).*10);
   %end
   
   fprintf('(%d) %s - %s\n',...
           predratings(ii,2),songdata.artist,songdata.title);
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
plot(ee,'r-');
hold on
plot(ev./sum(ev));
plot(svdcount,ee(svdcount),'o');
%plot([1 usercount],[1 1],'r--');
hold off

legend('pred cc','% eig variance');
xlabel('eigenvector rank');


%emtx=zmtx*u(:,1:svdcount);
goodknown=shuffle(find(evalmtx(:,useridx)>1));
goodknown=goodknown(1:min([100 length(goodknown)]));
badknown=shuffle(find(evalmtx(:,useridx)<-1));
badknown=badknown(1:min([100 length(badknown)]));
otherknown=shuffle(find(evalmtx(:,useridx)<=1 & evalmtx(:,useridx)>=-1));
otherknown=otherknown(1:min([100 length(otherknown)]));
goodpred=refids(length(predratings):-1:(length(predratings)-goodsongs+1));
badpred=refids(1:badsongs);

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

keyboard
