% function songdbsync(svdcount)
%
% svdcount - by default, chooses svd cutoff that gives min
%            prediction error on prior ratings.
%
function songdbsync(svdcount)

oo=0;
quitnow=0;

SYNCJLG=0;

while ~quitnow,
   dbopen('polka.isr.umd.edu','david','nine1997','music');
   
   % figure out date of last check
   sql='SELECT mdbsyncdate+0 as mdbsyncdate FROM dJuke WHERE playlist="nsl"';
   syncdata=mysql(sql);
   
   lastdate=num2str(syncdata.mdbsyncdate);
   
   if SYNCJLG,
      disp('checking for modified songs');
      
      % override to get all changes since branching
      %lastdate='20050729184300'
      
      cmd=['ssh david@fusiform.bic.berkeley.edu "mysql --host=sql' ...
           ' --user=david --password=1997nine --raw' ...
           ' --database=cell' ...
           ' --exec=''select * from dMusic' ...
           ' where lastmod>',lastdate,';''"'];
      [w,s]=unix(cmd);
      
      % populate a songdata structure. NOTE: all values are
      % currently saved as strings
      s=strsep(s,char(10),1);
      if length(s)>0,
         titles=strsep(strtrim(s{1}),char(9));
         fprintf('%d mod songs\n',length(s)-1);
      else
         disp('  no modified dMusic entries...');
      end
      
      for ii=2:length(s),
         if mod(ii,50)==0,
            fprintf('.');
         end
         ss=strsep(s{ii},char(9),1);
         if length(ss)>=length(titles),
            
            songdata=[];
            for jj=1:length(titles),
               if length(ss{jj})==0,
                  songdata=setfield(songdata,titles{jj},'');
               elseif strcmp(ss{jj},'NULL'),
                  songdata=setfield(songdata,titles{jj},'');
               else
                  songdata=setfield(songdata,titles{jj},ss{jj});
               end
            end
            sql=['SELECT * FROM dMusic WHERE id=',songdata.id];
            localsongdata=mysql(sql);
            
            if ~isempty(localsongdata),
               songdata.track
               if ~strcmp(localsongdata.title,songdata.title) | ...
                     ~strcmp(localsongdata.artist,songdata.artist) | ...
                     ~strcmp(localsongdata.album,songdata.album) | ...
                     localsongdata.track~=str2num(songdata.track),
                  fprintf('%d\n%s --> %s\n%s --> %s\n%s --> %s\n%d --> %d\n',...
                          ii,localsongdata.title,songdata.title, ...
                          localsongdata.artist,songdata.artist, ...
                          localsongdata.album,songdata.album, ...
                          localsongdata.track, ...
                          str2num(songdata.track));
                  yn=input('change [n]?','s');
                  if ~isempty(yn) & yn(1)=='y',
                     sql=['UPDATE dMusic SET',...
                          ' artist="',songdata.artist,'",',...
                          ' title="',songdata.title,'",',...
                          ' album="',songdata.album,'",',...
                          ' track=',songdata.track,...
                          ' WHERE id=',songdata.id]
                     mysql(sql);
                  end
               end
            else
               songdata
               fprintf('local song not found. pause before add\n');
               keyboard
               sqlinsert('dMusic','id',songdata.id,...
                         'title',songdata.title,...
                         'track',songdata.track,...
                         'artist',songdata.artist,...
                         'album',songdata.album,...
                         'file',songdata.file,...
                         'playsec',songdata.playsec,...
                         'genre',songdata.genre,...
                         'comment',songdata.comment,...
                         'note',songdata.note,...
                         'user',songdata.user,...
                         'id3bad',songdata.id3bad,...
                         'source',songdata.source);
            end
         end
      end
   else
      disp('skipping syncjlg');
   end
   
   evaldata=[]; % set evaldata to empty in case s is empty
   if SYNCJLG,
      disp('downloading evals');
      
      oo=oo+1;
      
      cmd=['ssh david@fusiform.bic.berkeley.edu "mysql --host=sql' ...
           ' --user=david --password=1997nine --raw --skip-column-names' ...
           ' --database=cell' ...
           ' --exec=''select now();''"'];
      [w,checkdate]=unix(cmd);
      checkdate=checkdate(1:end-1);
      
      cmd=['ssh david@fusiform.bic.berkeley.edu "mysql --host=sql' ...
           ' --user=david --password=1997nine --raw --skip-column-names' ...
           ' --database=cell' ...
           ' --exec=''select dEval.id,dEval.songid,dEval.user,dEval.rating,dEval.uid,dMusic.file from dEval,dMusic' ...
           ' where dEval.songid=dMusic.id AND dEval.lastmod>',lastdate,';''"'];
      [w,s]=unix(cmd);
      
      if w,
         disp('error connecting to jlg.berkeley.edu');
      elseif isempty(s),
         disp('  no new dEval entries');
      elseif s(1)=='<',
         disp(['error connecting to mdb:',s]);
      else
         s=strsep(s,char(10));
         for ii=2:length(s),
            ss=strsep(s{ii},char(9));
            if length(ss)>=6,
               evaldata(ii-1).id=ss{1};
               evaldata(ii-1).songid=ss{2};
               evaldata(ii-1).user=ss{3};
               evaldata(ii-1).rating=ss{4};
               evaldata(ii-1).uid=ss{5};
               evaldata(ii-1).file=ss{6};
            end
         end
         
         % checkdate set by output of neweval
         dbopen;
         mysql(['update dJuke set mdbsyncdate="',checkdate,'" WHERE playlist="nsl"']);
      end
   else
      evaldata=[];
   end
   
   if length(evaldata)>0,
      fprintf('  %d evals to check\n',length(evaldata));
      
      % connect to local sql
      dbopen('polka','david','nine1997','music');
      
      sql=['SELECT * FROM dEval'];
      oldevaldata=mysql(sql);
      idset=cat(1,oldevaldata.songid);
      
      for ii=1:length(evaldata),
         if mod(ii,1000)==0,
            fprintf('%d\n',ii);
         end
         
         match=0;
         mset=find(evaldata(ii).songid==idset);
         if length(mset)>0,
            jj=0;
            while (jj<length(mset) & ~match),
               jj=jj+1;
               if strcmp(oldevaldata(mset(jj)).user,evaldata(ii).user),
                  if oldevaldata(mset(jj)).rating==evaldata(ii).rating,
                     match=1;
                  else
                     match=2;
                  end
               end
            end
         end
         
         if ~match,
            sql=['SELECT * FROM dMusic where id=', ...
                 num2str(evaldata(ii).songid)];
            songdata=mysql(sql);
            
            if length(songdata)==0,
               disp('song not found locally!');
               keyboard
            elseif songdata.filemissing | songdata.crap | songdata.dup,
               match=1;
            end
         end
         
         if ~match
            fprintf('new rating: user %s, songid %d  %d\n',...
                    evaldata(ii).user,evaldata(ii).songid,...
                    evaldata(ii).rating);
            sqlinsert('dEval','user',evaldata(ii).user,...
                      'songid',evaldata(ii).songid,...
                      'rating',evaldata(ii).rating,'note','mdb');
         elseif match==2,
            fprintf('rating changed: user %s, songid %d %d --> %d\n',...
                    evaldata(ii).user,evaldata(ii).songid,...
                    oldevaldata(mset(jj)).rating,evaldata(ii).rating);
            sql=['UPDATE dEval set rating=',num2str(evaldata(ii).rating),...
                 ' WHERE songid=',num2str(evaldata(ii).songid),...
                 ' AND user="',evaldata(ii).user,'"'];
            mysql(sql);
         end
      end
   end
   
   %%
   % done checking. now do stats
   %%
   
   % get rating matrix
   minevals=8;
   [evalmtx,users,songid]=songevalmtx(minevals);
   
   if 0,
      disp('testing with top 100 only!');
      load /auto/k1/hayden/songids.mat
      dropidx=find(~ismember(songid,top100_songid));
      evalmtx(dropidx,useridx)=nan;
   end
   fullsongid=songid;
   fullevalmtx=evalmtx;
   
   % only take songs that have been rated by >=2 people
   minsongs=2;
   singids=find(sum(~isnan(evalmtx),2)>=minsongs);
   evalmtx=evalmtx(singids,:);
   songid=songid(singids);
   songcount=length(songid);
   
   mm=nanmean(evalmtx);
   ss=nanstd(evalmtx);
   evalmtx=(evalmtx-repmat(mm,[songcount,1])) ./ ...
           repmat(ss,[songcount,1]);
   fullevalmtx=(fullevalmtx-repmat(mm,[length(fullevalmtx),1])) ./ ...
       repmat(ss,[length(fullevalmtx),1]);
   
   % only include users that have rated by >=100 songs
   minuserevals=50;
   ecount=sum(~isnan(evalmtx));
   keepuseridx=find(ecount>=minuserevals);
   fullevalmtx=fullevalmtx(:,keepuseridx);
   tevalmtx=evalmtx(:,keepuseridx);
   ucount=size(tevalmtx,2);
   fprintf('minuserevals %d leaves %d users with ratings for stats\n',...
           minuserevals,ucount);
   for uuidx=1:length(keepuseridx),
      fprintf('%s: %d\n',users{keepuseridx(uuidx)},ecount(keepuseridx(uuidx)));
   end
   
   %
   % explore svd threshold options
   %
   
   %[ua,sa,va]=songdbspace(tevalmtx,10);
   
   if ~exist('svdcount','var'),
      svdcount=10;
   end
   fprintf('using %d principal components\n',svdcount);
   
   % compute correlation matrix for big users
   cc=zeros(ucount);
   n=zeros(ucount);
   for u1=1:ucount,
      for u2=u1:ucount,
         gii=find(~isnan(tevalmtx(:,u1)) & ~isnan(tevalmtx(:,u2)));
         if length(gii)>0,
            n(u1,u2)=length(gii);
            n(u2,u1)=n(u1,u2);
            cc(u1,u2)=mean(tevalmtx(gii,u1).*tevalmtx(gii,u2));
            cc(u2,u1)=cc(u1,u2);
         end
      end
   end
   cc(isnan(cc))=0;
   
   % regularize
   lown=find(n<10);
   cc(lown)=cc(lown).*n(lown)./10;
   
   % project songs onto 1...svdcount eigenvectors. this info
   % will be saved in dEigProj
   
   [u,s,v]=svd(cc);
   a=fullevalmtx;
   a(find(isnan(a)))=0;
   fullevaleigs=a*u(:,1:svdcount);
   %keyboard
   %fullevaleigs=fullevaleigs./repmat(sum(fullevaleigs.^2,2),[1 svdcount]);
   
   a=tevalmtx;
   a(find(isnan(a)))=0;
   evaleigs=a*u(:,1:svdcount);
   
   % now need to find optimal linear mapping from evaleigs to
   % evalmtx. this will be saved in gUserPrefs
   rr=1:svdcount;
   usercount=length(users);
   userproj=zeros(svdcount,usercount);
   invusers=[];
   
   for uidx=1:usercount,
      sidx=find(~isnan(evalmtx(:,uidx)));
      teval=evalmtx(sidx,uidx);
      eeac=evaleigs(sidx,rr)'*evaleigs(sidx,rr)./length(sidx);
      if length(sidx)<=15,
         fprintf('user %s: only %d ratings. can''t do inv\n',...
                 users{uidx},length(sidx));
      else
         fprintf('user %s: %d ratings. doing inv\n',...
                 users{uidx},length(sidx));
         userproj(rr,uidx)=eeac^-1 * (evaleigs(sidx,rr)'*teval)./length(sidx);
         invusers=[invusers uidx];
      end
   end
   
   %
   % save to local database
   %
   % save eigenvectors for projecting ratings into eigenvector space:
   disp('saving eignevectors');
   mysql('TRUNCATE dEig');
   for rr=1:length(keepuseridx),
      sql='INSERT INTO dEig (rr,cc,userid,u) VALUES ';
      for cc=1:svdcount,
         sql=[sql,sprintf('(%d,%d,"%s",%.5f),',rr,cc,users{keepuseridx(rr)},u(rr,cc))];
      end
      sql=sql(1:end-1);
      mysql(sql);
   end
   
   % fill recon vector in gUserPrefs;
   disp('saving reconstruction functions to local db');
   for ii=1:usercount,
      sql=['SELECT * FROM gUserPrefs WHERE userid="',users{ii},'"'];
      udata=mysql(sql);
      
      if length(udata)==0,
         fprintf('adding user %s to local db\n',users{ii});
         [xx,uidnum]=sqlinsert('gUserPrefs','userid',users{ii},...
                               'seclevel',0,'lab','mdb');
      else
         uidnum=udata.id;
      end
      
      sql='UPDATE gUserPrefs set ';
      for jj=1:svdcount,
         sql=[sql,sprintf('p%d=%.3f,',jj,userproj(jj,ii))];
      end
      sql=[sql(1:end-1), ...
           sprintf(',avgrating=%.3f,stdrating=%.3f',mm(ii),ss(ii)), ...
           ' WHERE userid="',users{ii},'"'];
      mysql(sql);
      
      % book-keeping: force uid to match local database
      sql=['UPDATE dEval SET uid=',num2str(uidnum),...
           ' WHERE user="',users{ii},'"',...
           ' AND (isnull(uid) OR uid<>',num2str(uidnum),')'];
      mysql(sql);
   end
   
   % a relic of a simpler time
   %disp('saving eigspace-projected song ratings');
   %eigdata=mysql(['SELECT dEigProj.* FROM dEigProj,dMusic',...
   %               ' WHERE dMusic.id=dEigProj.songid',...
   %               ' AND not(dMusic.dup) ORDER BY songid']);
   %esongid=cat(1,eigdata.songid);
   %allsongid=union(esongid,fullsongid);
   
   disp('saving eigspace-projected song ratings');
   
   % don't truncate.  Load in all the projections with negative
   % songids and then swap in one fell swoop.
   %mysql('TRUNCATE dEigProj');
   for ii=1:length(fullsongid),
      if mod(ii,400)==1,
         eigsql=['INSERT INTO dEigProj (songid',...
                 sprintf(',e%d',1:svdcount),...
                 sprintf(',b%d',1:svdcount),') VALUES '];
         %        sprintf(',st%d',1:5),...
      end
      
      %ee=find(esongid==allsongid(ii));
      %if ~isempty(ee),
      %   st=sprintf(',%.6f',eigdata(ee).st1,eigdata(ee).st2,...
      %              eigdata(ee).st3, eigdata(ee).st4,eigdata(ee).st5);
      %else
      %   st=',0,0,0,0,0';
      %end
      
      %ff=find(fullsongid==allsongid(ii));
      %if ~isempty(ff),
      %   eigsql=[eigsql,'(',num2str(allsongid(ii)),...
      %           sprintf(',%.5f',fullevaleigs(ff,:)),st,...
      %           sprintf(',%.5f',fullevaleigs(ii,:)),'),'];
      %else
      eigsql=[eigsql,'(',num2str(-fullsongid(ii)),...
              sprintf(',%.5f',fullevaleigs(ii,:)),...
              sprintf(',%.5f',fullevaleigs(ii,:)),'),'];
      %        st,...
      %end
      
      if mod(ii,400)==0 || ii==length(fullsongid),
         eigsql=eigsql(1:end-1);
         mysql(eigsql);
      end
   end
   
   % swap in new projections (quickly!)
   sql=['DELETE FROM dEigProj WHERE songid>0;'];
   mysql(sql);
   sql=' UPDATE dEigProj SET songid=-songid;';
   mysql(sql);
   
   eigsql=['INSERT INTO dEigProj (songid',sprintf(',e%d',1:svdcount),')' ...
           ' SELECT DISTINCT id as songid', sprintf(',0 e%d',1:svdcount) ...
           ' FROM dMusic LEFT JOIN dEigProj ON dMusic.id=dEigProj.songid',...
           ' WHERE isnull(dEigProj.songid)',...
           ' AND not(dup) AND not(filemissing) AND not(crap) ORDER BY songid'];
   mysql(eigsql);
   
   
   % clean up a couple things to help queries run faster:
   mysql('UPDATE dMusic set source="" where isnull(source)');
   mysql('UPDATE dMusic set genre="" where isnull(genre)');
   mysql('UPDATE dMusic set playsec=0 where isnull(playsec)');
   
   
   if 1,
      disp('not putting projections back up to jlg jukebox');
      
   elseif 1,
      disp('saving eigspace projected song ratings to jlg');
      cmd=['ssh david@fusiform.bic.berkeley.edu "mysql --host=sql' ...
           ' --user=david --password=1997nine --raw --skip-column-names' ...
           ' --database=cell' ...
           ' --exec=''TRUNCATE dEigProj''"'];
      [w,s]=unix(cmd);
      for ii=1:length(fullsongid),
         if mod(ii,50)==1,
            eigsql=['INSERT INTO dEigProj (songid',...
                    sprintf(',e%d',1:svdcount),') VALUES '];
         end
         eigsql=[eigsql,'(',num2str(fullsongid(ii)),...
                 sprintf(',%.3f',fullevaleigs(ii,:)),'),'];
         if mod(ii,50)==0 | ii==length(fullsongid),
            eigsql=eigsql(1:end-1)
            
            cmd=['ssh david@fusiform.bic.berkeley.edu "mysql --host=sql' ...
                 ' --user=david --password=1997nine' ...
                 ' --database=cell' ...
                 ' --exec=''', eigsql, ';''"'];
            [w,s]=unix(cmd);
            %mysql(eigsql);
         end
      end
      cmd=['ssh david@fusiform.bic.berkeley.edu "mysql --host=sql' ...
           ' --user=david --password=1997nine --raw --skip-column-names' ...
           ' --database=cell' ...
           ' --exec=''', eigsql, '''"'];
      [w,s]=unix(cmd);
   else
      %
      % save eigenvectors to remote database
      %
      disp('saving eigs to mdb');
      mysql('open','kamzik.org','kamzik_david','nine1997');
      mysql('use kamzik_music');
      
      % fill recon vector in gUserPrefs;
      disp('saving reconstruction functions to mtdb');
      for ii=1:usercount,
         sql='UPDATE gUserPrefs set ';
         for jj=1:svdcount,
            sql=[sql,sprintf('p%d=%.3f,',jj,userproj(jj,ii))];
         end
         sql=[sql(1:end-1), ...
              sprintf(',avgrating=%.3f,stdrating=%.3f',mm(ii),ss(ii)), ...
              ' WHERE userid="',users{ii},'"'];
         mysql(sql);
      end
      
      disp('saving eigspace projected song ratings');
      mysql('TRUNCATE dEigProj');
      for ii=1:length(songid),
         if mod(ii,250)==1,
            eigsql=['INSERT INTO dEigProj (songid',...
                    sprintf(',e%d',1:svdcount),...
                    sprintf(',b%d',1:svdcount),') VALUES '];
         end
         eigsql=[eigsql,'(',num2str(songid(ii)),...
                 sprintf(',%.5f',evaleigs(ii,:)),...
                 sprintf(',%.5f',evaleigs(ii,:)),'),'];
         if mod(ii,500)==0 | ii==length(songid),
            eigsql=eigsql(1:end-1);
            mysql(eigsql);
         end
      end
      
      sql=['SELECT genre,avg(e1) as e1,avg(e2) as e2,avg(e3) as e3,',...
           'avg(e4) as e4,avg(e5) as e5,avg(e6) as e6',...
           ' FROM dEigProj,dGenre',...
           ' WHERE dEigProj.songid=dGenre.songid GROUP BY genre'];
      genredata=mysql(sql);
      
      
      mysql('close');
      
      if 0
         [s,status]=...
             urlread(['http://www.kamzik.org/mdb/redoeigproj.php3?recalc=1']);
         if ~status,
            disp('error connecting to mdb');
         elseif s(1)=='<',
            disp(['error connecting to mdb:',s]);
         else
            disp(s);
         end
      end
      
      disp('pausing 0 seconds');
      pause(0);
   end
   
   songdbadjust;
   
   %disp('no new ratings. pausing 30 seconds');
   %pause(30);
   
   
   disp('quitting. infinite loop disabled.');
   quitnow=1;
end

return






