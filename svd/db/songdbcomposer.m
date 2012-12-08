dbopen;

sql=['SELECT id,file,artist,composer FROM dMusic where not(crap)',...
     ' AND not(filemissing) AND not(dup) AND isnull(composer)'];
filedata=mysql(sql);

for ii=1:length(filedata);
   ts=filedata(ii).file;
   [s,w]=unix(['id3v2 --list "',ts,'" |grep TCM']);
   
   if ~s && length(w)>0,
      w=strsep(w,':',1);
      w=strtrim(w{2});
      w=strrep(w,'"','\"');
      if max(double(w))>65000,
         w='';
      end
   else
      w='';
   end
   if isempty(w),
      %sql=['UPDATE dMusic SET composer=NULL WHERE id=',...
      %     num2str(filedata(ii).id)];
      %mysql(sql);
   else
      fprintf('current artist: %s -- found composer %s\n',...
              filedata(ii).artist,w);
      sql=['UPDATE dMusic SET composer="',w,'" WHERE id=',...
           num2str(filedata(ii).id)];
      mysql(sql);
   end
end

