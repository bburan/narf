dbopen('polka','david','nine1997','music');

adata=mysql(['select artist,album,min(id) as albumid,count(id) as' ...
             ' trackcount FROM dMusic where not(crap) and not(dup)',...
             ' and album<>"" and not(album like "unknown")',....
             ' and not(album like "Complete Beetho%") and isnull(albumid)',...
             ' AND source="kathleen"',...
             ' group by artist,album']);

for ii=1:length(adata),
   yn='';
   try,
      tdata=mysql(['SELECT * FROM dMusic where albumid>0 and not(dup) and not(crap)',...
                ' and artist="',adata(ii).artist,'" AND album="',adata(ii).album,'"']);
   catch,
      tdata=[];
   end
   try,
      trackcount=adata(ii).trackcount;
   catch,
      trackcount=0;
   end
   if trackcount==0,
      disp('skipping, trackcount error');
   elseif length(tdata)>0,
      sql=['UPDATE dMusic SET albumid=',num2str(tdata(1).albumid),...
           ' WHERE artist="',adata(ii).artist,'" AND album="',adata(ii).album,'"'];
      disp(sql);
      try,
         mysql(sql);
      catch
         disp('mysql error');
         keyboard
      end


   elseif trackcount>1,
      fprintf('artist/album: %s/%s, %d tracks.\n',adata(ii).artist,adata(ii).album,trackcount);
      sql=['UPDATE dMusic SET albumid=',num2str(adata(ii).albumid),...
           ' WHERE artist="',adata(ii).artist,'" AND album="',adata(ii).album,'"'];
      disp(sql);
      try,
         mysql(sql);
      catch
         disp('mysql error');
      end
   elseif 0,
      % check that track hasn't already been taken care of
      sql=['SELECT * FROM dMusic WHERE id=', num2str(adata(ii).albumid),' AND isnull(albumid)'];
      checkdata=mysql(sql);
      % find candidate matches:
      sql=['SELECT * FROM dMusic WHERE not(crap) AND not(dup) AND album="',adata(ii).album,'" order by track'];
      try,
         tdata=mysql(sql);
      catch
         disp('mysql error');
         tdata=[];
      end
      if length(checkdata)>0 && length(tdata)>1,
         for jj=1:length(tdata)
            if isempty(tdata(jj).track),
               tdata(jj).track=0;
            end
            fprintf('%d %s/%s\n',tdata(jj).track,tdata(jj).artist,tdata(jj).title);
         end
         fprintf('mixed artist(?) album: %s, %d tracks?\n',adata(ii).album,length(tdata));
         yn=input('merge (y/[n])? ','s');
         if strcmp(yn,'y'),
            sql=['UPDATE dMusic SET albumid=',num2str(adata(ii).albumid),...
              ' WHERE not(crap) AND not(dup) AND album="',adata(ii).album,'"']
            mysql(sql);
         end
      end
   end
   if strcmp(yn,'q'),
      break
   end
end

