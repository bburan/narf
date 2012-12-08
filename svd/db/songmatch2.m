% function songmatch2(svdcount,user);
%
% svdcount - by default chooses svd cutoff that gives min
%            prediction error on prior ratings.
%
function songmatch2(svdcount,user);

if ~exist('user','var');
   user=getenv('USER');
end

dbopen;

% get rating matrix
minevals=8;
[evalmtx,users,songid]=songevalmtx(minevals);

if 0,
   disp('testing with top 100 only!');
   load /auto/k1/hayden/songids.mat
   dropidx=find(~ismember(songid,top100_songid));
   evalmtx(dropidx,useridx)=nan;
end

minsongs=2;
singids=find(sum(~isnan(evalmtx),2)>=minsongs);
evalmtx=evalmtx(singids,:);
songid=songid(singids);
songcount=length(songid);

mm=nanmean(evalmtx);
ss=nanstd(evalmtx);
evalmtx=(evalmtx-repmat(mm,[songcount,1])) ./ ...
        repmat(ss,[songcount,1]);

minuserevals=100;
ecount=sum(~isnan(evalmtx));
tevalmtx=evalmtx(:,find(ecount>minuserevals));
ucount=size(tevalmtx,2);
fprintf('minuserevals %d leaves %d users with enough ratings for stats\n',...
        minuserevals,ucount);
if 0,
   jackcount=10;
   svdmax=10;
   tpredmtx=zeros(songcount,ucount,svdmax);
   for jj=1:jackcount,
      valbins=round(songcount*(jj-1)/jackcount+1):round(songcount*jj/jackcount);
      usebins=[1:round(songcount*(jj-1)/jackcount) ...
               round(songcount*jj/jackcount+1):songcount];
      
      cc=zeros(ucount);
      n=zeros(ucount);
      for u1=1:ucount,
         for u2=u1:ucount,
            n(u1,u2)=sum(~isnan(tevalmtx(usebins,u1).*tevalmtx(usebins,u2)));
            n(u2,u2)=n(u1,u1);
            cc(u1,u2)=nanmean(tevalmtx(usebins,u1).*tevalmtx(usebins,u2));
            cc(u2,u1)=cc(u1,u2);
         end
      end
      
      [u,s,v]=svd(cc);
      
      for uidx=1:ucount,
         ttevalmtx=tevalmtx(valbins,:);
         ttevalmtx(find(isnan(ttevalmtx)))=0;
         ttevalmtx(:,uidx)=0;
         eigmtx=ttevalmtx*u;
         for svdidx=1:svdmax,
            tpredmtx(valbins,uidx,svdidx)=...
                eigmtx(:,1:svdidx)*u(uidx,1:svdidx)';
         end
      end
   end
   
   valxc=zeros(svdmax,1);
   zzidx=find(~isnan(tevalmtx));
   a=tevalmtx(zzidx);
   for svdidx=1:svdmax,
      b=tpredmtx(:,:,svdidx);
      valxc(svdidx)=xcov(a,b(zzidx),0,'coeff');
   end
   
   svdcount=max(find(valxc==max(valxc)));
elseif ~exist('svdcount','var') | svdcount==0 | svdcount>ucount,
   svdcount=6;
end
fprintf('using %d principal components\n',svdcount);

cc=zeros(ucount);
n=zeros(ucount);
for u1=1:ucount,
   for u2=u1:ucount,
      n(u1,u2)=sum(~isnan(tevalmtx(:,u1).*tevalmtx(:,u2)));
      n(u2,u2)=n(u1,u1);
      cc(u1,u2)=nanmean(tevalmtx(:,u1).*tevalmtx(:,u2));
      cc(u2,u1)=cc(u1,u2);
   end
end

[u,s,v]=svd(cc);
a=tevalmtx;
a(find(isnan(a)))=0;
evaleigs=a*u(:,1:svdcount);

% now need to find optimal linear mapping from evaleigs to evalmtx
useridx=find(strcmp(user,users));
usercount=length(users);
if isempty(useridx),
   fprintf('user %s: no match\n',user);
   return
end

userproj=zeros(svdcount,usercount);

for uidx=1:usercount,
   sidx=find(~isnan(evalmtx(:,uidx)));
   teval=evalmtx(sidx,uidx);
   eeac=evaleigs(sidx,:)'*evaleigs(sidx,:)./length(sidx);
   if length(sidx)<svdcount,
      fprintf('user %s: only %d ratings. can''t do inv\n',...
              users{uidx},length(sidx));
   else
      userproj(:,uidx)=eeac^-1 * (evaleigs(sidx,:)'*teval)./length(sidx);
   end
end

% predict songs for useridx

notevalidx=1:size(evalmtx,1);
testevals=evaleigs(notevalidx,:)*userproj(:,useridx);

notevalidx=find(isnan(evalmtx(:,useridx)));
predevals=evaleigs(notevalidx,:)*userproj(:,useridx);

% note: sorting from worst to best!
[predratings,xx]=sortrows([predevals,songid(notevalidx)]);

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
   
   fprintf('%5.1f:',predratings(ii,1)*ss(useridx)+mm(useridx));
   %for jj=1:min([svdcount 5]),
   %   fprintf('%3.0f',emtx(refids(ii),jj).*10);
   %end
   
   fprintf('(%d) %s - %s\n',...
           predratings(ii,2),songdata.artist,songdata.title);
end


%dump best songs
outfile='/auto/music0/svdout/cppred';
fprintf('creating 100 best song dump file: %s\n',outfile);
fid=fopen(outfile,'w');
outpath='/auto/music0/svdout/best100';
for ii=[length(predratings):-1:(length(predratings)-100+1)],
   
   sql=['SELECT * FROM dMusic where id=',num2str(predratings(ii,2))];
   songdata=mysql(sql);
   title=songdata.title;
   title(find(title=='/'))='_';
   title(find(title=='"'))='_';
   artist=songdata.artist;
   artist(find(artist=='/'))='_';
   artist(find(artist=='"'))='_';

   fprintf(fid,'cp "%s" "%s/%0.3d %s - %s.mp3"\n',songdata.file,outpath,...
           length(predratings)-ii+1,artist,title);
end
fclose(fid);
unix(['chmod a+x ',outfile]);


keyboard

% fill recon vector in gUserPrefs;
disp('saving reconstruction functions to musicdb');
for ii=1:usercount,
   sql='UPDATE gUserPrefs set ';
   for jj=1:svdcount,
      sql=[sql,sprintf('p%d=%.3f,',jj,userproj(jj,ii))];
   end
   sql=[sql(1:end-1) , ' WHERE userid="',users{ii},'"'];
   mysql(sql);
end

disp('saving eigspace projected song ratings');

mysql('TRUNCATE dEigProj');
for ii=1:length(songid),
   if mod(ii,500)==1,
      eigsql=['INSERT INTO dEigProj (songid',sprintf(',e%d',1:svdcount),') VALUES '];
   end
   eigsql=[eigsql,'(',num2str(songid(ii)),sprintf(',%.5f',evaleigs(ii,:)),'),'];
   if mod(ii,500)==0 | ii==length(songid),
      eigsql=eigsql(1:end-1);
      mysql(eigsql);
   end
end

   
   
return


[cc,n,mm,ss]=songdbcorrmtx(evalmtx);

minsongs=2;
singids=find(sum(~isnan(evalmtx),2)>=minsongs);
evalmtx=evalmtx(singids,:);
songid=songid(singids);
songcount=length(songid);

evalmtx=(evalmtx-repmat(mm,[songcount,1])) ./ ...
        repmat(ss,[songcount,1]);



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


