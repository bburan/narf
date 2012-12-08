function [evalmtx,users,songid]=songevalmtx(minevals);

if ~exist('minevals','var'),
   minevals=8;
end

genreset={'"Rock","Classic Rock","Hard Rock","Punk","Alternative","Indie"',...
          '"Ambient","Easy Listening","Acoustic","Meditative"',...
          '"Christmas or winter holiday"',...
          '"Electronic","Techno","Electronic Rock","Trance","Goa","Drum & Bass"',...
          '"Jazz","NuJazz","Latin Jazz","Modern Jazz","Acid Jazz","Blues"',...
          '"Hard Rock","Punk","Hardcore","Metal"',...
         };
genrecount=length(genreset);

dbopen;

sql=['SELECT user,avg(rating)+std(rating)/4 as defrat FROM dEval group by user order by user'];
userdata=mysql(sql);
usercount=length(userdata);

%sql=['SELECT DISTINCT songid FROM dEval,dMusic',...
%     ' WHERE dEval.songid=dMusic.id',...
%     ' AND not(dup) AND not(filemissing) AND not(crap) ORDER BY songid'];
sql=['SELECT DISTINCT dMusic.id as songid FROM dMusic',...
     ' WHERE not(dup) AND not(filemissing) AND not(crap) ORDER BY dMusic.id'];
evaldata=mysql(sql);
songid=cat(1,evaldata.songid);
songcount=length(songid);

evalmtx=nan.*zeros(length(songid),usercount+genrecount);
users={};
for ii=1:usercount,
   users{ii}=userdata(ii).user;
   
   if 0,
   defrat=userdata(ii).defrat;
   
   sql=['SELECT DISTINCT dMusic.id as songid,dup,1 as rating ',...
        ' FROM dMusic LEFT JOIN dEval ',...
        ' ON (dEval.user="',userdata(ii).user,'"',...
        ' AND dMusic.id=dEval.songid)',...
        ' WHERE dMusic.source="',userdata(ii).user,'"',...
        ' AND isnull(dEval.songid)',...
        ' AND not(filemissing)',...
        ' ORDER BY songid'];
   srcdata=mysql(sql);
   ssongid=cat(1,srcdata.songid);
   sdupid=cat(1,srcdata.dup);
   ssongid(sdupid>0)=sdupid(sdupid>0);
   srating=cat(1,srcdata.rating);
   
   ff=find(ismember(songid,ssongid));
   evalmtx(ff,ii)=defrat;
   end
   
   sql=['SELECT DISTINCT songid,rating FROM dEval,dMusic',...
        ' WHERE dEval.user="',userdata(ii).user,'"',...
        ' AND dEval.songid=dMusic.id',...
        ' AND not(dup) AND not(crap) AND not(filemissing)',...
        ' ORDER BY songid'];
   evaldata=mysql(sql);
   
   usongid=cat(1,evaldata.songid);
   urating=cat(1,evaldata.rating);
   
   ff=find(ismember(songid,usongid));
   evalmtx(ff,ii)=urating;
   
end

for ii=1:length(genreset),
    sql=['SELECT DISTINCT songid FROM dGenre',...
        ' WHERE genre in (',genreset{ii},')',...
        ' ORDER BY songid'];
    evaldata=mysql(sql);
    
    usongid=cat(1,evaldata.songid);
    
    ff=find(ismember(songid,usongid));
    evalmtx(:,ii+usercount)=0;
    evalmtx(ff,ii+usercount)=1;
    
    gg=strsep(genreset{ii},',');
    gg=gg{1}(2:(end-1));
    users{ii+usercount}=gg;
end


if 0,
   disp('testing with top 100 only!');
   load /auto/k1/hayden/songids.mat
   dropidx=find(~ismember(songid,top100_songid));
   evalmtx(dropidx,useridx)=nan;
end

evalcount=sum(~isnan(evalmtx));
guidx=find(evalcount>=minevals);
evalmtx=evalmtx(:,guidx);
totusers=length(users); % before excisig those who haven't rated enough
users={users{guidx}};
usercount=length(users);

fprintf('%d songs, %d evals req''d. %d/%d users participating.\n',...
        songcount,minevals,usercount,totusers);



