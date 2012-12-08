% songdbid3.m  fix broken id3's
function songdbid3(topdir);

dbopen;

if exist('topdir','var'),
   sql=['SELECT * FROM dMusic where id3bad and not(crap) and not(dup)',...
        ' AND (isnull(artist) OR artist="")',... %OR artist="Unknown"
        ' AND not(filemissing)',...
        ' AND file like "%',topdir,'%"'];
else
   %sql=['SELECT * FROM dMusic where not(crap) and not(dup)',...
   %     ' AND (artist="Unknown Artist") AND not(filemissing)'];
   sql=['SELECT * FROM dMusic where id3bad and not(crap) and not(dup)',...
        ' AND (isnull(artist) OR artist="") AND not(filemissing)'];
end
filedata=mysql(sql);

for ii=1:length(filedata);
   ts=filedata(ii).file;
   
   fprintf('id3 tags bad for:\n%s\n',ts);
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
   if track==0,
      track=filedata(ii).track;
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
   
   fprintf('currently:\nart: %s\nalb: %s\ntrk: %d\ntit: %s\n',...
           filedata(ii).artist,filedata(ii).album,...
           filedata(ii).track,filedata(ii).title);
   fprintf('file: %s\n',filedata(ii).file);
   yn=input('add tags manually ([y]/n)? ','s');
   
   if isempty(yn) | strcmp(yn,'y'),
      if ~exist(ts,'file');
         disp('Fetching ...');
         tts=strrep(ts,' ','\ ');
         [s,w]=unix(['musicget2 "',tts,'"']);
      end
      [s,w]=unix(['mp3info -p "%S" "',ts,'"']);
      if s==0 & w(1)~='/',
         playsec=str2num(w);
      else
         [s,w]=unix(['du -sk "',ts,'"']);
         nn=min(find(w==' ' | w=='/'));
         playsec=round(str2num(w(1:(nn-1)))/1024*60);
      end
      
      if strcmp(artist(1:min([length(artist) 5])),'music'),
         artist=album;
         album='Unknown';
      end
      if strcmp(artist,'sivak') | strcmp(artist,'ben') | ...
            strcmp(artist,'elder') | strcmp(artist,'voytek'),
         artist='Unknown';
      end
      if ~isempty(findstr(title,' - ')),
         aa=findstr(title,' - ');
         if aa>1,
            artist=title(1:aa-1);
         end
         if aa<length(title)-2,
            title=title((aa+3):end);
         end
         album='';
      end
      while length(title)>0 & title(1)==' ',
         title=title(2:end);
      end
      while length(title)>0 & title(end)==' ',
         title=title(1:(end-1));
      end
      while length(artist)>0 & artist(1)==' ',
         artist=artist(2:end);
      end
      while length(artist)>0 & artist(end)==' ',
         artist=artist(1:(end-1));
      end
     
      yn='n';
      while ~isempty(yn) & yn(1)=='n',
         tt=input(['artist [',artist,']: '],'s');
         if ~isempty(tt),
            artist=tt;
            album='';
         end
         tt=input(['album [',album,']: '],'s');
         if ~isempty(tt),
            album=tt;
         end
         tt=input(['track [',num2str(track),']: '],'s');
         if ~isempty(tt),
            track=str2num(tt);
         end
         tt=input(['title [',title,']: '],'s');
         if ~isempty(tt),
            title=tt;
         end
         tt=input(['playsec [',num2str(playsec),']: '],'s');
         if ~isempty(tt),
            playsec=str2num(tt);
         end
         yn=input('proceed ([y]/n)? ','s');
      end
      
      if isempty(track),
         track='NULL';
      else
         track=num2str(track);
      end
      
      sql=['UPDATE dMusic SET ',...
           ' comment="id3s fixed by songdbid3.m",',...
           ' artist="',artist,'",',...
           ' album="',album,'",',...
           ' track=',track,',',...
           ' title="',title,'",',...
           ' playsec=',num2str(playsec),...
           ' WHERE id=',num2str(filedata(ii).id)];
      mysql(sql);
   end
end

   

return

% check existing files to see if they actually exist on disk
if CHECKOLD,
   
   infile='/auto/k1/david/tmp/songlist';

   disp('scanning lsd...');
   unix(['locate /home/music | grep [mM][pP]3 > ',infile]);
   disp('scanning mdma...');
   unix(['ssh mdma locate /local/MP3s | grep ".[mM][pP]3" >> ',infile]);
   
   fid=fopen(infile,'r');
   
   [a,count]=fscanf(fid,'%c');
   alen=length(a);
   fclose(fid);
   
   
   sql='SELECT id,file from dMusic;';
   filedata=mysql(sql);
   
   for ii=1:length(filedata),
      ts=filedata(ii).file;
      if strcmp(ts(1:10),'/auto/musi'),
         ts(1:5)='/home';
      end
      if strcmp(ts(1:8),'/auto/x0'),
         ts=['/local',ts(9:end)];
      end
      xx=findstr(a,ts);
      if isempty(xx),
         fprintf('not found! - %s\n',filedata(ii).file);
      end
      if mod(ii,1000)==0,
         fprintf('%d ',ii);
      end
   end
   fprintf('\n');
end


% look for new mp3s in the two major locations
if CHECKNEW,

   infile='/auto/k1/david/tmp/songlist';

   disp('scanning lsd...');
   unix(['locate /home/music | grep [mM][pP]3 > ',infile]);
   disp('skipping mdma scan...');
   %disp('scanning mdma...');
   %unix(['ssh mdma locate /local/MP3s | grep ".[mM][pP]3" >> ',infile]);
   
   fid=fopen(infile,'r');
   
   [a,count]=fscanf(fid,'%c');
   alen=length(a);
   fclose(fid);
   
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
      if strcmp(ts(1:5),'/home'),
         ts(1:5)='/auto';
      end
      if strcmp(ts(1:6),'/local'),
         ts=['/auto/x0',ts(7:end)];
      end
      
      s0=s1;
      
      %fprintf('%s\n',ts);
   
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
         [s,w]=unix(['song2db "',ts,'"']);
         
         % if filename at begining of result... means there was an error.
         if ~s & w(1)=='/',
            s=1;
         end
      
         if s==0,
            % found id3 info
            sql=['INSERT INTO dMusic (title,track,artist,album,playsec,',...
                 'user,genre,comment,file) VALUES (',...
                 w(1:end-1),',"',ts,'")'];
            fprintf('adding: %s\n',ts);
         else
            % no valid id3 info
            sql=['INSERT INTO dMusic (file,user,id3bad) VALUES (',...
                 '"',ts,'","',getenv('USER'),'",1)'];
            fprintf('adding (id3bad): %s\n',ts);
            
         end
         
         mysql(sql);
      end
   end
end




