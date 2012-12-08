% songdblens.m  fix broken id3 lengths

dbopen;

sql=['SELECT * FROM dMusic where not(crap)',...
     ' AND not(isnull(artist)) AND artist<>""',...
     ' AND (playsec<2 OR isnull(playsec))',...
     ' AND not(dup) AND not(filemissing)'];
filedata=mysql(sql);

for ii=1:length(filedata);
   ts=filedata(ii).file;
   
   fprintf('playsec empty for:\n%s\n',ts);
   
   extidx=sort(cat(2,findstr(ts,'!'),findstr(ts,'$'),...
                   findstr(ts,' '),...
                   findstr(ts,''''),findstr(ts,'&'),...
                   findstr(ts,'('),findstr(ts,')'),...
                   findstr(ts,'{'),findstr(ts,'}'),...
                   findstr(ts,'['),findstr(ts,']')...
                   ));
   for jj=1:length(extidx),
      ts=[ts(1:(extidx(jj)-1)),'\',ts(extidx(jj):end)];
      extidx((jj+1):end)=extidx((jj+1):end)+1;
   end
   
   [s,w]=unix(['mp3info -p "%S" ',ts]);
   
   if s==0 & w(1)~='/',
      playsec=str2num(w);
   else
      keyboard
      [s,w]=unix(['du -sk ',ts]);
      nn=min([find(w==' ' | w=='/') length(w)+1]);
      playsec=round(str2num(w(1:(nn-1)))/1024*60);
   end
   if isempty(playsec),
      playsec=0;
   end
   
   fprintf('setting playsec=%d\n',playsec);
   %keyboard
   sql=['UPDATE dMusic SET ',...
        ' playsec=',num2str(playsec),...
        ' WHERE id=',num2str(filedata(ii).id)];
   mysql(sql);
end

   

