% function mdbstats(svdcount,user);
%
% svdcount - by default, chooses svd cutoff that gives min
%            prediction error on prior ratings.
%
function mdbstats(svdcount,user);

if ~exist('user','var');
   user=getenv('USER');
end
SKIPSINGLES=1;
oo=0;
while 1,
   dbopen;
   
   sql='SELECT mdbsyncdate FROM dJuke';
   syncdata=mysql(sql);
   
   lastdate=syncdata.mdbsyncdate;
   
   mysql close
   
   disp('downloading evals');
   
   if 1
      evaldata=[]; % set evaldata to empty in case s is empty
      oo=oo+1;
      lastdate=[lastdate(1:10),'%20',lastdate(12:19)];
      [s,status]=...
          urlread(['http://www.kamzik.org/mdb/neweval.php3?lastdate=',...
                   lastdate]);
      if ~status,
         disp('error connecting to mdb');
      elseif s(1)=='<',
         disp(['error connecting to mdb:',s]);
      else
         fprintf('oo=%d: %s\n',oo,s);
         
         eval(s);
         
         % checkdate set by output of neweval
         dbopen;
         mysql(['update dJuke set mdbsyncdate="',checkdate,'"']);
      end
   else
      disp('connecting to mdb...');
      mysql('open','kamzik.org','kamzik_david','nine1997');
      mysql('use kamzik_music');
      
      sql=['SELECT * FROM dEval WHERE date_sub(lastmod, INTERVAL 2 HOUR)>="',lastdate,'"'];
      evaldata=mysql(sql);
      
      %disconnect from mdb
      mysql('close');
   end
   
   fprintf('%d evals to check\n',length(evaldata));
   
   if length(evaldata)>0,
      % connect to local sql
      dbopen;
      
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
            
            if songdata.filemissing | songdata.crap | songdata.dup,
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
      
      %%
      % done checking. now do stats
      %%
      
      minevals=8;
      [evalmtx,users,songid]=songevalmtx(minevals);
      usercount=length(users);
            
      %minsongs=2;
      %singids=find(sum(~isnan(evalmtx),2)>=minsongs);
      %evalmtx=evalmtx(singids,:);
      %songid=songid(singids);
      songcount=length(songid);
      
      [cc,n,mm,ss]=songdbcorrmtx(evalmtx);
      
      usercount=length(users);
      
      [u,s,v]=svd(cc);
      
      svdcount=8;
      sd=diag(s);
      sfrac=sum(sd(1:svdcount))./sum(sd);
      
      fprintf('%d eigenvectors (%.1f%% of variance)\n',svdcount,sfrac*100);
      uadj=u./sqrt(sfrac);
      
      %
      % test: save to local database
      %
      disp('saving eigs locally');
      mysql('TRUNCATE dEig');
      for rr=1:usercount,
         sql=['select * from gUserPrefs where userid="',users{rr},'"'];
         udata=mysql(sql);
         
         if length(udata)>0,
            uidnum=udata.id;
         else
            sqlinsert('gUserPrefs','userid',users{rr},'seclevel',0,'lab','mdb');
            sql=['select * from gUserPrefs where userid="',users{rr},'"'];
            udata=mysql(sql);
            uidnum=udata.id;
         end
         
         sql=sprintf(['UPDATE gUserPrefs set avgrating=%.3f,stdrating=%.3f',...
                      ' WHERE id=%d'],mm(rr),ss(rr),uidnum);
         mysql(sql);
         
         % book-keeping: force uid to match local database
         sql=['UPDATE dEval SET uid=',num2str(uidnum),...
              ' WHERE user="',users{rr},'"',...
              ' AND (isnull(uid) OR uid<>',num2str(uidnum),')'];
         mysql(sql);
         
         sql='INSERT INTO dEig (rr,cc,uid,u) VALUES ';
         for cc=1:svdcount,
            sql=[sql,sprintf('(%d,%d,%d,%.5f),',rr,cc,uidnum,uadj(rr,cc))];
         end
         sql=sql(1:end-1);
         mysql(sql);
      end
      
      if 0,
         % don't bother with this. it doesn't seem to speed things up
         disp('saving eigproj locally');
         tevalmtx(find(isnan(tevalmtx)))=0;
         umtx=tevalmtx*uadj(:,1:svdcount);
         mysql('TRUNCATE dEigProj');
         for ii=1:length(songid),
            if mod(ii,500)==1,
               eigsql='INSERT INTO dEigProj (songid,cc,uproj) VALUES ';
            end
            for cc=1:svdcount,
               eigsql=[eigsql,sprintf('(%d,%d,%.5f),',songid(ss(ii)),...
                                      cc,umtx(ii,cc))];
            end
            if mod(ii,500)==0 | ii==length(songid),
               eigsql=eigsql(1:end-1);
               mysql(eigsql);
            end
         end
      end
      
      mysql('close');
      
      %
      % save eigenvectors to remote database
      %
      disp('saving eigs to mdb');
      mysql('open','kamzik.org','kamzik_david','nine1997');
      mysql('use kamzik_music');
      
      sql='select * from gUserPrefs';
      udata=mysql(sql);
      
      eigsql='INSERT INTO dEig (rr,cc,uid,u) VALUES ';
      for rr=1:usercount,
         umatch=0;
         for uu=1:length(udata),
            if strcmp(udata(uu).userid,users{rr}),
               umatch=uu;
               %fprintf('umatch=1 for userid=%s uid=%d\n',users{rr},udata(uu).id);
            end
         end
         
         if umatch>0,
            uidnum=udata(umatch).id;
         else
            sqlinsert('gUserPrefs','userid',users{rr},'seclevel',0,'lab','mdb');
            sql=['select * from gUserPrefs where userid="',users{rr},'"'];
            tudata=mysql(sql);
            uidnum=tudata.id;
         end
         
         % book-keeping: make sure means and stds are up to date
         % shouldn't need to do this regularly!
         sql=sprintf(['UPDATE gUserPrefs set avgrating=%.3f,stdrating=%.3f',...
                      ' WHERE id=%d'],mm(rr),ss(rr),uidnum);
         %mysql(sql);
         
         % book-keeping: force uid to match local database
         sql=['UPDATE dEval SET uid=',num2str(uidnum),...
              ' WHERE user="',users{rr},'"',...
              ' AND (isnull(uid) OR uid<>',num2str(uidnum),')'];
         mysql(sql);
         
         %append the eig insert string
         for cc=1:svdcount,
            eigsql=[eigsql,sprintf('(%d,%d,%d,%.5f),',rr,cc,uidnum,uadj(rr,cc))];
         end
      end
      % actually save the eigs
      eigsql=eigsql(1:end-1);
      
      % don't erase til ready to add new eigs immediately!
      mysql('TRUNCATE dEig');
      mysql(eigsql);
      mysql('close');
      
      [s,status]=...
          urlread(['http://www.kamzik.org/mdb/redoeigproj.php3?recalc=1']);
      if ~status,
         disp('error connecting to mdb');
      elseif s(1)=='<',
         disp(['error connecting to mdb:',s]);
      else
         disp(s);
      end
      
      disp('pausing 0 seconds');
      pause(0);
   else
      disp('pausing 30 seconds');
      pause(30);
   end
end

return

% test out eig stuff

dbopen;
keyboard

sid=15;
sql=['SELECT cc,sum(dEig.u * (dEval.rating-gUserPrefs.avgrating)/ ',...
     '              gUserPrefs.stdrating) as uproj',...
     ' FROM (dEig INNER JOIN dEval ON dEig.uid=dEval.uid)', ...
     ' INNER JOIN gUserPrefs ON gUserPrefs.id=dEval.uid',...
     ' WHERE dEval.songid=',num2str(sid),...
     ' GROUP BY cc'];
projdata=mysql(sql);
testspace=cat(1,projdata.uproj);

uid=1; % david locally
useridx=3; % david in alphabetical order

projmat=[];
for cc=1:svdcount,
   sql=['SELECT songid,sum(dEig.u * (dEval.rating-gUserPrefs.avgrating)/ ',...
        '              gUserPrefs.stdrating) as uproj,',...
        ' count(dEval.rating) as ratcount,',...
        ' sum(dEval.uid=',num2str(uid),') as urat',...
        ' FROM (dEig INNER JOIN dEval ON dEig.uid=dEval.uid)', ...
        ' INNER JOIN gUserPrefs ON gUserPrefs.id=dEval.uid',...
        ' WHERE cc=',num2str(cc),...
        ' GROUP BY songid',...
        ' HAVING ratcount>1 ORDER BY songid'];
   projdata=mysql(sql);
   projmat=[projmat cat(1,projdata.uproj)];
   songidin=cat(1,projdata.songid);
   ratcount=cat(1,projdata.ratcount);
   urat=cat(1,projdata.urat);
   
end

predmat=projmat * uadj(:,1:svdcount)';


keyboard




