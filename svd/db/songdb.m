% function songdb(topdir);
%
% maintain dMusic table
%
% probably will only work if you're logged in to mescaline?
%
% topdir is the path to the files you want to search for mp3s
% (remember to run updatedb on mescaline first!)
%
function songdb(topdir,user,machine);

dbopen;

disp('songdb.m: starting (remember to run updatedb on mescaline first!)');
disp('          (also make sure mp3info is installed on this computer!)');
if ~exist('topdir','var');
   origtopdir='/auto/music0/';
   topdir='/data/music/';
   CHECKNEW=1;
   DOALL=1;
else
   origtopdir=topdir;
   CHECKNEW=1;
   DOALL=0;
   if strcmp(topdir(1:12),'/auto/music/'),
      topdir=['/home/',topdir(13:end)];
      if ~exist('machine','var');
         machine='mescaline';
      end
   elseif strcmp(topdir(1:13),'/auto/music2/'),
      topdir=['/home2/',topdir(14:end)];
      if ~exist('machine','var');
         machine='mescaline';
      end
   elseif strcmp(topdir(1:13),'/auto/music0/'),
      topdir=['/mnt/data/music/',topdir(14:end)];
      if ~exist('machine','var');
         machine='valium';
      end
   end
end

if ~exist('user','var');
   user=getenv('USER');
else
   fprintf('user=%s\n',user);
end

% look for new mp3s in the two major locations
if CHECKNEW,
   disp('checking for new files not in db:');
   
   infile=['/tmp/',getenv('USER'),'.songlist'];
   
   if DOALL,
      disp('scanning music1...');
      %unix(['ssh lsd locate /home/music | grep "[mM][pP]3" > ',infile]);
      unix(['ssh mescaline locate /home | grep ".[mM][pP]3" >> ',infile]);
      disp('scanning music2...');
      unix(['ssh mescaline locate /home2 | grep ".[mM][pP]3" >> ',infile]);
   else
      disp(['scanning :',origtopdir,' ...']);
      
      host=strsep(getenv('HOST'),'.');
      host=host{1};
      if 1 | strcmp(machine,host),
         
         [b,ttd]=basename(origtopdir);
         [s,w]=unix(['locate -U ',ttd,' -o /tmp/songdb.tmp']);
         [s,w]=unix(['locate -d /tmp/songdb.tmp "',origtopdir,...
                     '" | grep ".[mM][pP]3" >> ',infile]);
         [s,w]=unix(['\rm /tmp/songdb.tmp']);
         
      else
         
         [b,ttd]=basename(topdir);
         
         disp(['ssh ',machine,' ''locate "',topdir,...
               '" | grep ".[mM][pP]3" >> ',infile,'''']);
         
         [s,w]=unix(['ssh ',machine,' '' locate -U ',ttd,...
                     ' -o /tmp/songdb.tmp''']);
         [s,w]=unix(['ssh ',machine,...
                     ' ''locate -d /tmp/songdb.tmp "',topdir,...
                     '" | grep ".[mM][pP]3" >> ',infile,'''']);
         [s,w]=unix(['ssh ',machine,' '' \rm /tmp/songdb.tmp''']);
         
         %[s,w]=unix(['ssh ',machine,' ''locate "',topdir,...
         %      '" | grep ".[mM][pP]3" >> ',infile,'''']);
      end
   end
   
   fid=fopen(infile,'r');
   
   [a,count]=fscanf(fid,'%c');
   alen=length(a);
   fclose(fid);
   
   delete(infile);
   
   s0=0;
   while s0<alen,
      
      rangelen=100;
      s1=[];
      while isempty(s1);
         
         searchrange=(s0+1):min([s0+rangelen alen]);
         s1=s0+min(find(a(searchrange)==10));
         
         if s0+rangelen>alen & isempty(s1),
            s1=alen+1;
         end
         
         rangelen=rangelen+100;
      end
      
      ts=a(s0+1:s1-1);
      if strcmp(ts(1:6),'/home/'),
         ts=['/auto/music',ts(6:end)];
      end
      if strcmp(ts(1:6),'/home2'),
         ts=['/auto/music2',ts(7:end)];
      end
      if strcmp(ts(1:15),'/mnt/data/music'),
         ts=['/auto/music0',ts(16:end)];
      end
      
      s0=s1;
      
      if ~strcmp(ts(end-3:end),'.mp3') & ~strcmp(ts(end-3:end),'.MP3'),
         % force skip if ext is not .mp3
         songdata=1;
      else
         sql=['select * from dMusic where file="',ts,'"'];
         songdata=mysql(sql);
      end
      
      if length(songdata)>0,
         fprintf('.');
      else
         % insert backslash before special characters ...
         % eg, replace ! with \! in filename
         extidx=sort(cat(2,findstr(ts,'!'),findstr(ts,'$'),...
                         findstr(ts,' '),findstr(ts,''''),...
                         findstr(ts,'('),findstr(ts,')'),...
                         findstr(ts,'{'),findstr(ts,'}'),...
                         findstr(ts,'['),findstr(ts,']'),...
                         findstr(ts,'&')));
         ts2=ts;
         for jj=1:length(extidx),
            ts2=[ts2(1:(extidx(jj)-1)),'\',ts2(extidx(jj):end)];
            extidx((jj+1):end)=extidx((jj+1):end)+1;
         end
         
         [s,w]=unix(['song2db ',ts2]);
         
         % if filename at begining of result... means there was an error.
         if ~s & w(1)=='/',
            s=1;
         end
         
         if s==0,
            % found id3 info
            sql=['INSERT INTO dMusic (title,track,artist,album,playsec,',...
                 'user,genre,comment,file,source,dateadded) VALUES (',...
                 w(1:end-1),',"',ts,'","',user,'",now())']
            fprintf('\nadding: %s\n',ts);
            [res,aff,songid]=mysql(sql);
            
            sql=['SELECT id,genre FROM dMusic WHERE id=',num2str(songid)];
            gendata=mysql(sql);
            if ~strcmp(gendata.genre,''),
               sqlinsert('dGenre',...
                         'songid',songid,...
                         'genre',gendata.genre);
            end
         else
            fprintf('\nid3 tags bad for:\n%s\n',ts);
            
            if 0
               yn=input('add tags manually ([y]/n)? ','s');
            else
               yn='n';
            end
            if isempty(yn) | strcmp(yn,'y'),
            
               [b,p]=filebasename(ts);
               title=b(1:end-4);
               
               if length(title)>2,
                  if ~isempty(str2num(title(1:2))),
                     track=str2num(title(1:2));
                     title=title(3:end);
                  else
                     track=0;
                  end
               end
               
               p=p(1:end-1);
               [b,p]=filebasename(p);
               album=b;
               
               p=p(1:end-1);
               [b,p]=filebasename(p);
               artist=b;
               
               artist(find(artist=='_'))=' ';
               album(find(album=='_'))=' ';
               title(find(title=='_'))=' ';
               while (title(1)==' ' | title(1)=='-') & length(title)>1,
                  title=title(2:end);
               end
               
               [s,w]=unix(['mp3info -p "%S" ',ts2]);
               if s==0 & w(1)~='/',
                  playsec=str2num(w);
               else
                  [s,w]=unix(['du -sk ',ts2]);
                  nn=min(find(w==' ' | w=='/'));
                  playsec=round(str2num(w(1:(nn-1)))/1024*60);
               end
               source=user;
               
               tt=input(['artist [',artist,']: '],'s');
               if ~isempty(tt),
                  artist=tt;
               end
               tt=input(['album [',album,']: '],'s');
               if ~isempty(tt),
                  album=tt;
               end
               tt=input(['track [',num2str(track),']: '],'s');
               if ~isempty(tt),
                  track=str2num(tt);
                  if isempty(track),
                     track=0;
                  end
               end
               
               tt=input(['title [',title,']: '],'s');
               if ~isempty(tt),
                  title=tt;
               end
               tt=input(['playsec [',num2str(playsec),']: '],'s');
               if ~isempty(tt),
                  playsec=str2num(tt);
                  if isempty(playsec),
                     playsec=0;
                  end
               end
               if strcmp(user,'david'),
                  tt=input(['source [',source,']: '],'s');
                  if ~isempty(tt),
                     source=tt;
                  end
               end
               
               dbopen; % in case connection has timed out while
                       % user was entering things by hand
               sql=['INSERT INTO dMusic (artist,album,track,title,playsec,',...
                    'file,source,user,comment,id3bad,dateadded) VALUES (',...
                    '"',artist,'","',album,'",',num2str(track),',',...
                    '"',title,'",',num2str(playsec),',"',ts,'","',user,...
                    '","',user,'","manual id3 by songdb",1,now())'];
               mysql(sql);
               
            else
               % no valid id3 info
               fprintf('adding (id3bad): %s\n',ts);
               
               dbopen; % in case connection has timed out while
                       % user was entering things by hand
               
               %sqlinsert('dMusic',...
                %         'file',ts,'user',user,'id3bad',1);
               sql=['INSERT INTO dMusic (',...
                    'file,source,user,comment,id3bad,dateadded) VALUES (',...
                    '"',ts,'","',user,...
                    '","',user,'","manual id3 by songdb",1,now())'];
               mysql(sql);
            end
         end
      end
   end
   fprintf('\n');
end




