% function mdbstats(svdcount);
%
% svdcount - by default, chooses svd cutoff that gives min
%            prediction error on prior ratings.
%
function mdbstats2(svdcount);

oo=0;
while 1,
   dbopen;
   
   % figure out date of last check
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
   end
   
   if length(evaldata)>0,
      fprintf('%d evals to check\n',length(evaldata));
      
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
      
      if ~exist('svdcount','var'),
         svdcount=6;
      end
      fprintf('using %d principal components\n',svdcount);
      
      % compute correlation matrix for big users
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
      
      % project songs onto 1...svdcount eigenvectors. this info
      % will be saved in dEigProj
      [u,s,v]=svd(cc);
      a=tevalmtx;
      a(find(isnan(a)))=0;
      evaleigs=a*u(:,1:svdcount);

      % now need to find optimal linear mapping from evaleigs to
      % evalmtx. this will be saved in gUserPrefs
      usercount=length(users);
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
      
      %
      % test: save to local database
      %
      
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
      
      disp('saving eigspace projected song ratings');
      mysql('TRUNCATE dEigProj');
      for ii=1:length(songid),
         if mod(ii,500)==1,
            eigsql=['INSERT INTO dEigProj (songid',...
                    sprintf(',e%d',1:svdcount),') VALUES '];
         end
         eigsql=[eigsql,'(',num2str(songid(ii)),...
                 sprintf(',%.5f',evaleigs(ii,:)),'),'];
         if mod(ii,500)==0 | ii==length(songid),
            eigsql=eigsql(1:end-1);
            mysql(eigsql);
         end
      end
      
      mysql('close');
      
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
         if mod(ii,500)==1,
            eigsql=['INSERT INTO dEigProj (songid',...
                    sprintf(',e%d',1:svdcount),') VALUES '];
         end
         eigsql=[eigsql,'(',num2str(songid(ii)),...
                 sprintf(',%.5f',evaleigs(ii,:)),'),'];
         if mod(ii,500)==0 | ii==length(songid),
            eigsql=eigsql(1:end-1);
            mysql(eigsql);
         end
      end
      
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
   else
      disp('no new ratings. pausing 30 seconds');
      pause(30);
   end
end

return






