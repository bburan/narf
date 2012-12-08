


dbopen('polka.isr.umd.edu','david','nine1997','music');

sql=['SELECT dMusic.*, dGenre.genre as genre2 from dMusic LEFT JOIN ' ...
     'dGenre ON dMusic.id=dGenre.songid where album="" and ' ...
     'dMusic.genre="Trance" and source="tac"',...
     ' AND not(file like "/auto/music4/tac%") group by dMusic.id'];
d=mysql(sql);
for ii=1:length(d),
    newgenre=d(ii).genre2;
    newsource='';
    if ~isempty(findstr(d(ii).file,'/voytek/')),
        newsource='voytek';
    elseif ~isempty(findstr(d(ii).file,'/stephen_download/')),
        newsource='david';
    elseif ~isempty(findstr(d(ii).file,'/mp3download/')),
        newsource='david';
    elseif ~isempty(findstr(d(ii).file,'/stephen/')),
        newsource='david';
    elseif ~isempty(findstr(d(ii).file,'/stephen_random/')),
        newsource='david';
    elseif ~isempty(findstr(d(ii).file,'/DJ Entropy/')),
        newsource='david';
    elseif ~isempty(findstr(d(ii).file,'/Four Tet/')),
        newsource='david';
    elseif ~isempty(findstr(d(ii).file,'/ryan/')),
        newsource='prenger';
    elseif ~isempty(findstr(d(ii).file,'/kendrick/')),
        newsource='kendrick';
    elseif ~isempty(findstr(d(ii).file,'/kdonald/')),
        newsource='kdonald';
    elseif ~isempty(findstr(d(ii).file,'/hayden/')),
        newsource='hayden';
    elseif ~isempty(findstr(d(ii).file,'/ben/')),
        newsource='hayden';
    elseif ~isempty(findstr(d(ii).file,'/Ben/')),
        newsource='hayden';
    elseif ~isempty(findstr(d(ii).file,'/Suspicious Minds/')),
        newsource='hayden';
    elseif ~isempty(findstr(d(ii).file,'/Disneyland/')),
        newsource='hayden';
    elseif ~isempty(findstr(d(ii).file,'/kathleen/')),
        newsource='kathleen';
    elseif ~isempty(findstr(d(ii).file,'/tac/')),
        newsource='tac';
    elseif ~isempty(findstr(d(ii).file,'/dab/')),
        newsource='dab';
    elseif ~isempty(findstr(d(ii).file,'/elder/')),
        newsource='elder';
    elseif ~isempty(findstr(d(ii).file,'/englitz/')),
        newsource='englitz';
    elseif ~isempty(findstr(d(ii).file,'/knathan/')),
        newsource='knathan';
    elseif ~isempty(findstr(d(ii).file,'/wafting/')),
        newsource='wafting';
    elseif ~isempty(findstr(d(ii).file,'/sivak/')),
        newsource='sivak';
    elseif ~isempty(findstr(d(ii).file,'/kanold/')),
        newsource='kanold';
    elseif ~isempty(findstr(d(ii).file,'/mikesMusic/')),
        newsource='mikeg';
    elseif ~isempty(findstr(d(ii).file,'/mikeg/')),
        newsource='mikeg';
    elseif ~isempty(findstr(d(ii).file,'/mainland/')),
        newsource='mainland';
    elseif ~isempty(findstr(d(ii).file,'/mainland/')),
        newsource='mainland';
    elseif ~isempty(findstr(d(ii).file,'/shinji/')),
        newsource='shinji';
    elseif ~isempty(findstr(d(ii).file,'/mnima/')),
        newsource='mnima';
    elseif ~isempty(findstr(d(ii).file,'/alexzo/')),
        newsource='alexzo';
    elseif ~isempty(findstr(d(ii).file,'/jamie-misc/')),
        newsource='mazer';
    elseif ~isempty(findstr(d(ii).file,'/Skatalites/')),
        newsource='mazer';
    elseif ~isempty(findstr(d(ii).file,'/daniel/')),
        newsource='daniel';
    elseif ~isempty(findstr(d(ii).file,'/shin/')),
        newsource='shin';
    elseif ~isempty(findstr(d(ii).file,'/hazeleyes/')),
        newsource='hazeleyess87';
    elseif ~isempty(findstr(d(ii).file,'/bfen/')),
        newsource='laxwrestler';
    elseif ~isempty(findstr(d(ii).file,'/kghose/')),
        newsource='gghose';
    elseif ~isempty(findstr(d(ii).file,'/avi/')),
        newsource='avi';
    elseif ~isempty(findstr(d(ii).file,'/therick/')),
        newsource='therick';
    elseif ~isempty(findstr(d(ii).file,'/jzemsk/')),
        newsource='jzemsk';
    elseif ~isempty(findstr(d(ii).file,'/majid/')),
        newsource='majid';
    elseif ~isempty(findstr(d(ii).file,'/nima/')),
        newsource='mnima';
    elseif ~isempty(findstr(d(ii).file,'/willmore/')),
        newsource='willmore';
    elseif d(ii).id>=12357 && d(ii).id<26289 &&...
            ~isempty(findstr(d(ii).file,'/music2/')),
        newsource='ekm';
 
    
    end
    
    

    if ~isempty(newsource),
        sql=['UPDATE dMusic set genre="',newgenre,'",'...
             'source="',newsource,'" WHERE id=',num2str(d(ii).id)];
        fprintf('file %s : %s\n',d(ii).file,sql);
        mysql(sql);
    end
end

        
        



dupdata=mysql('select dEval.*,dMusic.dup FROM dEval inner join dMusic on dEval.songid=dMusic.id where dMusic.dup>0;');

for ii=1:length(dupdata),
    sql=['UPDATE dEval set songid=',num2str(dupdata(ii).dup),...
         ' WHERE songid=',num2str(dupdata(ii).songid)]
    mysql(sql);
end

missingdata=mysql('select dMusic.* FROM dEval inner join dMusic on dEval.songid=dMusic.id where dMusic.filemissing;');

for ii=1:length(missingdata),
    fprintf('%s - %s - %s\n',missingdata(ii).artist,...
           missingdata(ii).title,missingdata(ii).file);
    sql=['UPDATE dMusic SET filemissing=0 WHERE id=',...
         num2str(missingdata(ii).id)];
    mysql(sql);
    
end

    
    sql=['SELECT * FROM dMusic WHERE artist="',...
         missingdata(ii).artist,'"',...
         ' AND title="',missingdata(ii).title,'"',...
         ' AND not(filemissing) AND not(dup)'];
    mdata=mysql(sql);
    
    
    sql=['UPDATE dEval set songid=',num2str(dupdata(ii).dup),...
         ' WHERE songid=',num2str(dupdata(ii).songid)]
    mysql(sql);
end



dbopen

spokendata=mysql(['select distinct songid from dGenre where genre in' ...
            ' ("Audio Book","Book","Spoken Word")']);

for ii=1:length(spokendata),
   sql=['UPDATE dMusic set skiprand=1 where id=',...
        num2str(spokendata(ii).songid)];
   mysql(sql);
   
end

dbopen
sql=['select * from dMusic where skiprand=0 AND playsec>1800'];
longdata=mysql(sql);
for ii=1:length(longdata),
   sql=['UPDATE dMusic set skiprand=1 where id=',...
        num2str(longdata(ii).id)];
   mysql(sql);
end




dbopen;

sql=['select * from dMusic where filemissing=2 ORDER BY file'];
missingsongdata=mysql(sql);

for ii=1:length(missingsongdata),
   if exist(missingsongdata(ii).file,'file'),
      fprintf('file now exists: %s\n',missingsongdata(ii).file);
      sql=['UPDATE dMusic SET filemissing=0 WHERE id=',...
           num2str(missingsongdata(ii).id)];
      mysql(sql);
   end
end



sql=['select * from dMusic where source="kdonald" and filemissing=2 ORDER BY file'];
missingsongdata=mysql(sql);

% dump to file
f=fopen('/tmp/missingsongs_kevin.txt','w');
for ii=1:length(missingsongdata),
   fprintf(f,'%s\n',missingsongdata(ii).file);
end



return

dbopen;

sql=['select * from dMusic where file like "/auto/music3/hayden/%"',...
     'and not(dup) and not(filemissing);']

newsongdata=mysql(sql);

for ii=1:length(newsongdata),
   sql=['SELECT * FROM dMusic WHERE artist="',newsongdata(ii).artist,...
        '" AND title="',newsongdata(ii).title,'"',...
        ' AND filemissing'];
   oldsongdata=mysql(sql);
   if length(oldsongdata)==1,
      fprintf('%d: match for %s -- %s\n',ii,...
              newsongdata(ii).artist,newsongdata(ii).title)
      sql=['UPDATE dEval SET songid=',num2str(newsongdata(ii).id),...
           ' WHERE songid=',num2str(oldsongdata.id)];
      %mysql(sql);
      sql=['UPDATE dPlaylist SET songid=',num2str(newsongdata(ii).id),...
           ' WHERE songid=',num2str(oldsongdata.id)];
      %mysql(sql);
   else
      
      fprintf('%d: NO match for %s -- %s\n',ii,...
              newsongdata(ii).artist,newsongdata(ii).title)
      
   end
end

      



return


sql=['SELECT * FROM dMusic',...
     ' WHERE artist like "various%";'];
songdata=mysql(sql);

for ii=1:length(songdata),
   
   ts=songdata(ii).file;
   [b,p]=filebasename(ts);
   title=b(1:end-4);
   
   if length(title)>2,
      if ~isempty(str2num(title(1:2))),
         track=str2num(title(1:2));
         title=title(3:end);
      else
         track=songdata(ii).track;
      end
   end
   
   album=songdata(ii).album;
   artist=songdata(ii).artist;
   
   title(find(title=='_'))=' ';
   while (title(1)==' ' | title(1)=='-') & length(title)>1,
      title=title(2:end);
   end
   
   fprintf('%d: %s\n',songdata(ii).id,ts);
   
   yn='n';
   while (~isempty(yn) & strcmp(yn,'n')),
      
      fprintf('art: %s->%s\nalb: %s->%s\ntrk: %d->%d\ntit: %s->%s\n',...
              songdata(ii).artist,artist,...
              songdata(ii).album,album,...
              songdata(ii).track,track,...
              songdata(ii).title,title);
      
      if ~isempty(yn) & strcmp(yn,'n'),
         
         tt=input(['artist [',artist,']: '],'s');
         if ~isempty(tt),
            artist=tt;
         end
         %tt=input(['album [',album,']: '],'s');
         %if ~isempty(tt),
         %   album=tt;
         %end
         %tt=input(['track [',num2str(track),']: '],'s');
         %if ~isempty(tt),
         %   track=str2num(tt);
         %end
         tt=input(['title [',title,']: '],'s');
         if ~isempty(tt),
            title=tt;
         end
         
         yn=input('procede ([y]/n/s(kip))? ','s');
      end
   end
   if ~strcmp(yn,'s'),
      sql=['UPDATE dMusic SET ',...
           ' comment="id3s fixed by songdbid3.m",',...
           ' artist="',artist,'",',...
           ' album="',album,'",',...
           ' track=',num2str(track),',',...
           ' title="',title,'"',...
           ' WHERE id=',num2str(songdata(ii).id)];
      
      mysql(sql);
   end  
end



return

dbopen;

sql=['SELECT dMusic.id,count(dEval.rating) as ratcount',...
     ' FROM dMusic LEFT JOIN dEval',...
     ' ON dMusic.id=dEval.songid',...
     ' GROUP BY dMusic.id'];
songdata=mysql(sql);

for ii=1:length(songdata),
   if mod(ii,1000)==0, 
      ii
   end
   if songdata(ii).ratcount==0,
      sql=['UPDATE dMusic set temprating=NULL WHERE id=', ...
           num2str(songdata(ii).id)];
      mysql(sql);
   else
      sql=['SELECT round(avg((dEval.rating-gUserPrefs.avgrating)/',...
           'gUserPrefs.stdrating)*sqrt(count(dEval.id)),3) as rating',...
           ' FROM dEval INNER JOIN gUserPrefs',...
           ' ON dEval.user=gUserPrefs.userid',...
           ' WHERE dEval.songid=',num2str(songdata(ii).id)];
      evaldata=mysql(sql);
      
      sql=['UPDATE dMusic set temprating=',num2str(evaldata.rating),...
           ' WHERE id=',num2str(songdata(ii).id)];
      
      %keyboard
      
      mysql(sql);
   end
end

return



dbopen;

sql=['select album,artist,count(id) as count,min(id) as minid from dMusic',...
     ' where track in (0,1) and album<>"" and not(isnull(album))',...
     ' group by artist,album having count(id)>3'];
albumdata=mysql(sql);

for ii=1:length(albumdata);
   sql=['SELECT * FROM dMusic WHERE album="\',albumdata(ii).album,'\""',...
        ' AND artist="',albumdata(ii).artist,'"'];
   songdata=mysql(sql);
   
   for jj=1:length(songdata),
      tf=songdata(jj).file;
      b=basename(tf);
      tt=sscanf(b,'%d');

      if ~isempty(tt),
         if length(tt)>1,
            b
            tt
            keyboard
         end
         if songdata(jj).track~=tt,
            fprintf('%s - %s : track -> %d\n',songdata(jj).artist,...
                    songdata(jj).title,tt);
            
            sql=['UPDATE dMusic set track=',num2str(tt),...
                 ' WHERE id=',num2str(songdata(jj).id)];
            mysql(sql);
         end
      end
   end
end




return

%% TRY TO EXTEND TRUNCATED MP3INFO USING FILE NAMES

dbopen;

sql=['SELECT * FROM dMusic WHERE length(title)=30'];
songdata=mysql(sql);

for ii=1:length(songdata),
   ts=songdata(ii).file;
   [b,p]=filebasename(ts);
   title=b(1:end-4);
   
   if length(title)>2,
      if ~isempty(str2num(title(1:2))),
         track=str2num(title(1:2));
         title=title(3:end);
      else
         track=songdata(ii).track;
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
   
   fprintf('%d: %s\n',songdata(ii).id,ts);
   
   yn='n';
   while (~isempty(yn) & strcmp(yn,'n')),
      
      fprintf('art: %s ->%s\nalb: %s ->%s\ntrk: %d ->%d\ntit: %s ->%s\n',...
              songdata(ii).artist,artist,...
              songdata(ii).album,album,...
              songdata(ii).track,track,...
              songdata(ii).title,title);
      
      yn=input('procede ([y]/n/s(kip))? ','s');
   
      if ~isempty(yn) & strcmp(yn,'n'),
         
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
         end
         tt=input(['title [',title,']: '],'s');
         if ~isempty(tt),
            title=tt;
         end
      end
   end
   if ~strcmp(yn,'s'),
      sql=['UPDATE dMusic SET ',...
           ' comment="id3s fixed by songdbid3.m",',...
           ' artist="',artist,'",',...
           ' album="',album,'",',...
           ' track=',num2str(track),',',...
           ' title="',title,'"',...
           ' WHERE id=',num2str(songdata(ii).id)];
      
      mysql(sql);
   end
end

return

