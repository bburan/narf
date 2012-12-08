% songdbdups.m  - find and mark duplicates
%
%function songdbdups(topdir);

dbopen;

if exist('afilt','var'),
   sql=['SELECT artist,title, count(id) as repcount FROM dMusic',...
        ' WHERE not(crap) and not(dup) AND not(filemissing)',...
        ' AND artist like "',afilt,'"',...
        ' GROUP BY artist,title HAVING repcount>1'];
else
   sql=['SELECT artist,title, count(id) as repcount FROM dMusic',...
        ' WHERE not(crap) and not(dup) AND not(filemissing)',...
        ' GROUP BY artist,title HAVING repcount>1'];
end
filedata=mysql(sql);

for ii=1:length(filedata);
   sql=['SELECT * FROM dMusic',...
        ' WHERE artist="',filedata(ii).artist,'"',...
        ' AND title="',filedata(ii).title,'"',...
        ' AND not(id3bad) and not(crap) and not(dup) AND not(filemissing)'];
   dupdata=mysql(sql);
   
   cmd='0';
   while ~isempty(cmd),
   
      fprintf('%s - %s\n',filedata(ii).artist,filedata(ii).title);
      for jj=1:length(dupdata),
         if dupdata(jj).dup>0,
            dupstr=sprintf(' ** dup-->%d **',dupdata(jj).dup);
         else
            dupstr='';
         end
         if isempty(dupdata(jj).track),
            trk=-1;
         else
            trk=dupdata(jj).track;
         end
         if isempty(dupdata(jj).source),
            src='';
         else
            src=dupdata(jj).source;
         end
         
         fprintf(' %d (%d) %s - Trk %d - %d sec (%s) %s\n',jj,dupdata(jj).id,...
                 dupdata(jj).album,trk,dupdata(jj).playsec,src,dupstr);
         fprintf('   %s\n',dupdata(jj).file);
      end
      
      cmd=input(sprintf('dup (1-%d,[q]):',length(dupdata)),'s');
      
      cmd=str2num(cmd);
      if cmd>0 & cmd<=length(dupdata),
         
         if dupdata(cmd).dup>0,
            dupdata(cmd).dup=0;
         else
            cmd2=input(sprintf('rep by (1-%d,[q]):',length(dupdata)),'s');
            cmd2=str2num(cmd2);
            if cmd2>0 & cmd2<=length(dupdata) & cmd2~=cmd,
               fprintf('%d is dup by %d\n',dupdata(cmd).id, ...
                       dupdata(cmd2).id);
               dupdata(cmd).dup=dupdata(cmd2).id;
            end
         end
      end
   end
   
   for jj=1:length(dupdata),
      if dupdata(jj).dup>0,
         fprintf('commiting changes for songid=%d\n',dupdata(jj).id);
         sql=['UPDATE dMusic set dup=',num2str(dupdata(jj).dup),...
              ' WHERE id=',num2str(dupdata(jj).id)];
         mysql(sql);
         
         sql=['SELECT * FROM dEval WHERE songid=',num2str(dupdata(jj).id)];
         evaldata=mysql(sql);
         
         for kk=1:length(evaldata),
            sql=['SELECT * FROM dEval',...
                 ' WHERE songid=',num2str(dupdata(jj).dup),...
                 ' AND user="',evaldata(kk).user,'"'];
            e2data=mysql(sql);
            if length(e2data)==0,
               sql=['UPDATE dEval SET',...
                    ' songid=',num2str(dupdata(jj).dup),...
                    ' WHERE id=',num2str(evaldata(kk).id)];
               mysql(sql);
            end
         end
         
         sql=['UPDATE dPlaylist SET',...
              ' songid=',num2str(dupdata(jj).dup),...
              ' WHERE songid=',num2str(dupdata(jj).id)];
         mysql(sql);
      end
   end
end

return


%% tidying utility.

cmd='a%';
while ~isempty(cmd),
   dbopen;
   
   sql=['SELECT count(id) as cnt,artist FROM dMusic ',...
        ' WHERE artist like "',cmd,'"',...
        ' GROUP BY artist ORDER BY artist'];
   adata=mysql(sql);
   
   for ii=1:length(adata),
      fprintf('%3d %s (%d)\n',ii,adata(ii).artist,adata(ii).cnt);
   end
   cin=input(sprintf('\nEdit whom [quit]: '),'s');
   
   nn=str2num(cin);
   if isempty(nn),
      cmd=cin;
   elseif nn>0 & nn<=length(adata),
      oldstr=adata(nn).artist;
      newstr=input(sprintf('Change %s to: ',oldstr),'s');
      if ~isempty(newstr),
         sql=['UPDATE dMusic SET artist="',newstr,'"',...
              ' WHERE artist="',oldstr,'"'];
         mysql(sql);
      end
   end
end
