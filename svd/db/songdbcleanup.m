%% DELETE MISSING FILES AND ASSOCIATED RANKINGS/REQUEST RECORDS/GENRES

CHECKFORMISSING=0;
DLMISSING=0;
DELETEDUPS=0;
FIXDUP1=0;
FIXDUPEVAL=0;
FIXRAP=1;

if FIXRAP,
    rdata=mysql(['select dMusic.id,bitrate,dup,filemissing,',...
                 ' dGenre.genre as cgenre,file',...
                 ' FROM dMusic LEFT JOIN dGenre ON dMusic.id=dGenre.songid',...
                 ' WHERE source="hayden" and dMusic.genre="Rap"',...
                 ' and album="" and not(dup)']);
    for ii=1:length(rdata),
        ff=rdata(ii).file;
        fprintf('%d/%d\n',ii,length(rdata));
        fprintf('%s\n',ff);
        fprintf('genre=%s\n',rdata(ii).cgenre);
        
        if ~isempty(findstr(ff,'/auto/music/stephen_random/')),
            s='david';
        elseif ~isempty(findstr(ff,'/auto/music0/ben/')),
            s='hayden';
        elseif ~isempty(findstr(ff,'/auto/music2/ben/')),
            s='hayden';
        elseif ~isempty(findstr(ff,'/auto/music3/hayden/')),
            s='hayden';
        elseif ~isempty(findstr(ff,'/auto/music4/hayden/')),
            s='hayden';
        elseif ~isempty(findstr(ff,'/auto/music/kendrick')),
            s='kendrick';
        elseif ~isempty(findstr(ff,'/auto/music0/kendrick')),
            s='kendrick';
        elseif ~isempty(findstr(ff,'/auto/music/wafting/')),
            s='wafting';
        else
            s=input('source? ','s');
        end
        sql=['UPDATE dMusic SET source="',s,'",genre="',...
             rdata(ii).cgenre,'"',...
             ' WHERE id=',num2str(rdata(ii).id)]
        mysql(sql);
        
    end
    

end

if FIXDUPEVAL,
    dbopen('polka.isr.umd.edu','david','nine1997','music');
    disp('checking for dup chains');
    sql=['SELECT dMusic.*,dM2.dup dup2,dM2.artist as artist2,dM2.title as title2'...
         ' FROM dMusic,dMusic dM2'...
         ' WHERE dMusic.dup=dM2.id AND dMusic.dup>1 AND dM2.dup>1'];
    dupdups=mysql(sql);
    for ii=1:length(dupdups),
        sql=['UPDATE dMusic SET dup=',num2str(dupdups(ii).dup2),...
             ' WHERE dMusic.id=',num2str(dupdups(ii).id)]
        mysql(sql);
    end
    
    dbopen('polka.isr.umd.edu','david','nine1997','music');
    disp('checking ratings of songs marked dup');
    sql=['SELECT dEval.*,dMusic.dup'...
         ' FROM dEval,dMusic'...
         ' WHERE dEval.songid=dMusic.id AND dMusic.dup>1'];
    dupevals=mysql(sql);
    
    for ii=1:length(dupevals),
        sql=['SELECT * FROM dMusic WHERE id=',num2str(dupevals(ii).songid)];
        odata=mysql(sql);
        sql=['SELECT * FROM dMusic WHERE id=',num2str(dupevals(ii).dup)];
        ndata=mysql(sql);
        fprintf('%d %s - %s - %s --> %d %s - %s - %s\n',...
                odata.id,odata.artist,...
                odata.album,odata.title,...
                ndata.id,ndata.artist,...
                ndata.album,ndata.title);
        
        sql=['SELECT * FROM dEval WHERE songid=',...
             num2str(ndata.id),...
             ' AND user="',dupevals(ii).user,'"'];
        ddata=mysql(sql);
        if ~isempty(ddata),
            fprintf('user %s dup ratings old: %d new: %d\n',...
                    dupevals(ii).user,dupevals(ii).rating,ddata.rating);
            sql=['DELETE FROM dEval',...
                 ' WHERE id=',num2str(dupevals(ii).id)]
        else
            fprintf('user %s rating: %d\n',...
                    dupevals(ii).user,dupevals(ii).rating);
            sql=['UPDATE dEval SET songid=',num2str(ndata.id),...
                 ' WHERE id=',num2str(dupevals(ii).id)]
        end

        %pause;
        mysql(sql);
    end
    
end

if FIXDUP1,
    dbopen('polka.isr.umd.edu','david','nine1997','music');
    disp('checking for the existence of duplicate files');
    filedata=mysql(['SELECT id,file,artist,title,album,dup,bitrate',...
                    ' FROM dMusic',...
                    ' WHERE dup=1;']);
    for ii=1:length(filedata),
        fprintf('%d : %s - %s - %s\n',filedata(ii).id,...
                filedata(ii).artist,filedata(ii).album,...
                filedata(ii).title);
        
        sql=['SELECT * FROM dMusic',...
             ' WHERE artist="',filedata(ii).artist,'"',...
             ' AND album="',filedata(ii).album,'"',...
             ' AND title="',filedata(ii).title,'"',...
             ' AND dup=0'];
        matchdata=mysql(sql);
        if isempty(matchdata),
            sql=['SELECT * FROM dMusic',...
                 ' WHERE artist="',filedata(ii).artist,'"',...
                 ' AND title="',filedata(ii).title,'"',...
                 ' AND dup=0'];
            matchdata=mysql(sql);
        end
        if length(matchdata)==1,
            fprintf('found single match: %d: %s - %s - %s\n',...
                    matchdata.id,matchdata.artist,matchdata.album,...
                    matchdata.title);
            sql=['UPDATE dMusic SET dup=',num2str(matchdata.id),...
                 ' WHERE id=',num2str(filedata(ii).id)];
            mysql(sql);
        else
            fprintf('found %d matches\n',length(matchdata));

            %keyboard
        end
    end
end


if DELETEDUPS,
   dbopen('polka.isr.umd.edu','david','nine1997','music');
   disp('checking for the existence of duplicate files');
   filedata=mysql(['SELECT id,file,dup,bitrate,playsec from dMusic',...
                   ' WHERE dup AND file like "/auto/music%";']);
   
   for ii=1:length(filedata),
      ts=filedata(ii).file;
      delid=filedata(ii).id;
      % reload old data in case it has changed during this loop
      sql=['SELECT * FROM dMusic WHERE id=',num2str(filedata(ii).id)];
      oldfiledata=mysql(sql);
      if exist(ts,'file') && filedata(ii).dup>1,
         fprintf('(id %d) dup found: %s (%d kb/s %d sec)\n',...
                 filedata(ii).id,ts,...
                 filedata(ii).bitrate,filedata(ii).playsec);
         sql=['SELECT * FROM dMusic WHERE id=',num2str(oldfiledata.dup)];
         repfiledata=mysql(sql);
         if ~isempty(repfiledata) && repfiledata.dup>0,
            disp('updating orphaned dup link');
            sql1=['UPDATE dMusic set dup=',...
                  num2str(repfiledata.dup),' WHERE id=',...
                  num2str(oldfiledata.id)]
            %pause
            mysql(sql1);
            sql=['SELECT * FROM dMusic WHERE id=',...
                 num2str(repfiledata.dup)];
            repfiledata=mysql(sql);
            
         end
         if ~isempty(repfiledata) && exist(repfiledata.file,'file'),
            fprintf('(id %d) match:  %s (%d kb/s %d sec)\n',...
                    repfiledata.id,repfiledata.file,...
                    repfiledata.bitrate,repfiledata.playsec);
            keptbr=repfiledata.bitrate;
            oldbr=oldfiledata.bitrate;
            if isempty(oldbr),oldbr=0; end
            if (keptbr==-1 && oldbr>192) || ...
               (keptbr>0 && keptbr<oldbr) || ...
               keptbr==-2,
               fprintf('bitrate higher for file marked as dup?\n');
               if ~isempty(repfiledata.albumid),
                  sql1=['UPDATE dMusic set dup=0,albumid=',...
                        num2str(repfiledata.albumid),' WHERE id=',...
                        num2str(oldfiledata.id)]
               else
                  sql1=['UPDATE dMusic set dup=0 WHERE id=',...
                        num2str(oldfiledata.id)]
               end
               sql2=['UPDATE dMusic set dup=',num2str(oldfiledata.id),...
                     ' WHERE id=', num2str(repfiledata.id),...
                     ' OR dup=', num2str(repfiledata.id)]
               
               %pause
               mysql(sql1);
               mysql(sql2);
               mysql(['UPDATE dEval set songid=',num2str(oldfiledata.id),...
                      ' WHERE songid=',num2str(repfiledata.id)]);
               ts=repfiledata.file;
               delid=repfiledata.id;
            end
            
            disp(['I want to delete: ',ts]);
            sql=['UPDATE dMusic set filemissing=1',...
                 ' WHERE id=',num2str(delid)];
            if isempty(repfiledata.albumid) ||...
                isempty(oldfiledata.albumid) || ...
                    oldfiledata.albumid~=repfiledata.albumid,
                disp('but there is a mismatch in albumid');
                pause(0.01);
          else
                delete(ts);
                mysql(sql);
                disp('DELETING AND PAUSING 0.5 SEC!');
                pause(0.5);
            end
         else
            fprintf('dup file not found %s\nproposed change:\n',...
                    repfiledata.file);
            sql1=['UPDATE dMusic set filemissing=0,dup=0 WHERE id=',...
                 num2str(oldfiledata.id)]
            sql2=['UPDATE dMusic set filemissing=1,dup=',...
                 num2str(oldfiledata.id),...
                 ' WHERE id=',num2str(repfiledata.id)]
            disp(['Nothing to delete in this case since file already ' ...
                  'missing']);
            %pause
            mysql(sql1);mysql(sql2);
            mysql(['UPDATE dEval set songid=',num2str(oldfiledata.id),...
                 ' WHERE songid=',num2str(repfiledata.id)]);
         end
         
      end
      if mod(ii,500)==0,
         fprintf('%d ',ii);
      end
   end
   fprintf('\n');
   disp('need to run songdbsync now!');
   
   return
end

if CHECKFORMISSING,
   dbopen('polka.isr.umd.edu','david','nine1997','music');
   disp('checking for missing files. this takes a while!');
   sql=['SELECT id,file,replace(file,'' '',''\\ '') as file2',...
        ' FROM dMusic WHERE not(crap) AND not(dup) AND not(filemissing);'];
   filedata=mysql(sql);
   
   for ii=1:length(filedata),
      ts=filedata(ii).file;
      if ~exist(ts,'file'),
         fprintf('not found! - %s\n',filedata(ii).file);
         
         ts2=strrep(ts,'/auto/','/data/');
         if exist(ts2,'file'),
             fprintf('but bebop copy does: %s\n',ts2);
             
             ['cp "',ts2,'" "',ts,'"']
             pause
         end
         
         
         %dbopen('polka','david','nine1997','music');
         %if DLMISSING,
         %   [b,p]=basename(filedata(ii).file2); 
         %   disp(['musicget2 "',p,'*"']);

            %keyboard
            %   [s,w]=unix(['musicget2 "',p,'*"']);
            %   disp(w);
            %else
            %   sql=['UPDATE dMusic SET filemissing=1 WHERE id=',num2str(filedata(ii).id)];
            %   mysql(sql);
            %end
      end
      if mod(ii,500)==0,
         fprintf('%d ',ii);
      end
   end
   fprintf('\n');
end

return

if 0,
   sql=['SELECT * FROM dMusic where filemissing'];
   filedata=mysql(sql);
   if length(filedata)>0,
      fprintf('found %d missing files to delete.\n',length(filedata));
      
      if 1,
         disp('skipping delete');
      else
         if length(filedata)>10,
            for ii=1:length(filedata)
               fprintf('%d %s: %s\n',filedata(ii).id,filedata(ii).artist,...
                       filedata(ii).title);
            end
            disp('this seems like a lot. pausing.');
            keyboard
         end
         
         sql=['DELETE FROM dMusic WHERE filemissing'];
         mysql(sql);
      end
   end
end

sql=['SELECT dPlaylist.* FROM dPlaylist LEFT JOIN dMusic',...
     ' ON dPlaylist.songid=dMusic.id',...
     ' WHERE isnull(dMusic.id)'];
baddata=mysql(sql);
if length(baddata)>0,
   fprintf('deleting %d orphaned playlist entries.\n',length(baddata));
   for ii=1:length(baddata),
      sql=['DELETE FROM dPlaylist where id=',num2str(baddata(ii).id)];
      mysql(sql);
   end
end

sql=['SELECT dGenre.* FROM dGenre LEFT JOIN dMusic',...
     ' ON dGenre.songid=dMusic.id',...
     ' WHERE isnull(dMusic.id)'];
baddata=mysql(sql);
if length(baddata)>0,
   fprintf('deleting %d orphaned genre entries.\n',length(baddata));
   for ii=1:length(baddata),
      sql=['DELETE FROM dGenre where id=',num2str(baddata(ii).id)];
      mysql(sql);
   end
end

sql=['SELECT dEval.* FROM dEval LEFT JOIN dMusic',...
     ' ON dEval.songid=dMusic.id',...
     ' WHERE isnull(dMusic.id)'];
baddata=mysql(sql);
if length(baddata)>0,
   fprintf('deleting %d orphaned eval entries.\n',length(baddata));
   for ii=1:length(baddata),
      sql=['DELETE FROM dEval where id=',num2str(baddata(ii).id)];
      mysql(sql);
   end
end

% correct any invalid reqcounts, based on whether request was
% canceled or duplicate songs were merged
sql=['SELECT dMusic.id,dMusic.reqcount,count(dPlaylist.id) as newreqcount',...
     ' FROM dPlaylist INNER JOIN dMusic',...
     ' ON dPlaylist.songid=dMusic.id',...
     ' WHERE dPlaylist.user<>"xmms"',...
     ' GROUP BY dMusic.id;'];
reqdata=mysql(sql);

mismatchcount=0;
for ii=1:length(reqdata),
   if reqdata(ii).reqcount~=reqdata(ii).newreqcount,
      mismatchcount=mismatchcount+1;
      sql=['UPDATE dMusic set reqcount=',num2str(reqdata(ii).newreqcount),...
           ' WHERE id=',num2str(reqdata(ii).id)];
      mysql(sql);
   end
end

if mismatchcount>0,
   fprintf('fixed %d reqcount mismatches.\n',mismatchcount);
end

return

sql=['select dMusic.title,dMusic.artist,dMusic.dup,dEval.* from dEval,dMusic',...
     ' WHERE dEval.songid=dMusic.id and (dup)'];
rdata=mysql(sql)

for ii=1:length(rdata),
   if rdata(ii).dup==1,
      sql=['SELECT * FROM dMusic where artist="',rdata(ii).artist,'"',...
           ' AND title="',rdata(ii).title,'"',...
           ' AND not(dup)'];
      ddata=mysql(sql);
      if length(ddata)==1,
         
         fprintf('%d. fixing dup link for %s/%s (id %d) from 1 to %d\n',...
                 ii,rdata(ii).artist,rdata(ii).title,rdata(ii).songid,...
                 ddata.id);
         sql=['UPDATE dMusic set dup=',num2str(ddata.id),...
              ' WHERE id=',num2str(rdata(ii).songid)];
         mysql(sql);
         goodsongid=ddata.id;
      else
        fprintf('%d. no obvious dup fix for %s/%s (id %d) from songid=1\n',...
                 ii,rdata(ii).artist,rdata(ii).title,rdata(ii).songid);
          
         goodsongid=0;
      end
   else
      
      sql=['SELECT * FROM dMusic where id=',num2str(rdata(ii).dup)];
      ddata=mysql(sql);
      goodsongid=ddata.id;
   end
   
   if goodsongid,
   
      sql=['select * FROM dEval WHERE songid=',num2str(ddata.id),...
           ' AND user="',rdata(ii).user,'"'];
      edata=mysql(sql);
      if isempty(edata),
         fprintf('%d. %s/%s (id %d) dupped to %s/%s (id %d)\n',...
                 ii,rdata(ii).artist,rdata(ii).title,rdata(ii).songid,...
                 ddata.artist,ddata.title,ddata.id);
         sql=['UPDATE dEval set songid=',num2str(ddata.id),...
              ' WHERE id=',num2str(rdata(ii).id)];
         mysql(sql);
      else
         fprintf('%d. %s/%s (id %d) rating already exists for dup %s/%s (id %d)\n',...
                 ii,rdata(ii).artist,rdata(ii).title,rdata(ii).songid,...
                 ddata.artist,ddata.title,ddata.id);
         
      end
   end
end


return

sql=['select dMusic.* from dMusic,dEval where (dup or filemissing or crap) and dEval.songid=dMusic.id'];
rdata=mysql(sql);

for ii=1:length(rdata),
   sql=['SELECT * FROM dMusic WHERE id<>',num2str(rdata(ii).id),...
        ' AND title="',rdata(ii).title,'"',...
        ' AND artist="',rdata(ii).artist,'"'];
   linkdata=mysql(sql);
   if ~isempty(linkdata),
      fprintf('%d: %d points to bad %d\n',ii,linkdata.id,rdata(ii).id);
   end
end




sql=['select dm2.id,dm2.dup,dMusic.bitrate,dm2.bitrate,' ...
     ' dMusic.artist,dMusic.title,dMusic.source from dMusic,dMusic' ...
      ' dm2 where (dMusic.filemissing or dMusic.crap or dMusic.dup) and dm2.dup=dMusic.id;'];
rdata=mysql(sql);



for ii=1:length(rdata),
   fprintf('%d: dMusic.id in (%d,%d)\n',ii,rdata(ii).id,rdata(ii).dup);
end

   sql=['UPDATE dMusic set dup=0 WHERE id=',num2str(rdata(ii).id)]
   sql=['UPDATE dMusic set dup=',num2str(rdata(ii).id),...
        ' WHERE id=',num2str(rdata(ii).dup),...
        ' OR dup=',num2str(rdata(ii).dup)]
   sql=['UPDATE dEval SET songid=',num2str(rdata(ii).id),...
        ' WHERE songid=',num2str(rdata(ii).dup)]
