% songdblens.m  fix broken id3 lengths

dbopen('polka','david','nine1997','music');

CHECKBR=1;

if CHECKBR,
   sql=['SELECT id,file,bitrate,playsec,year FROM dMusic where not(crap)',...
        ' AND not(filemissing) AND not(dup) AND cat1=0 AND id>30000',...
        ' ORDER BY id'];
  
else
   sql=['SELECT id,file,bitrate,playsec,year FROM dMusic where not(crap)',...
        ' AND (playsec<1 OR isnull(playsec) OR isnull(bitrate) OR bitrate=0 OR bitrate=-2 OR year=0 OR isnull(year))',...
        ' AND not(filemissing) AND not(dup) AND id>100000'];
end
filedata=mysql(sql);

fprintf('%d songs to check\n',length(filedata));

for ii=1:length(filedata);
   ts=filedata(ii).file;
   extidx=sort(cat(2,findstr(ts,'!'),findstr(ts,'$'),...
                   findstr(ts,' '),...
                   findstr(ts,';'),...
                   findstr(ts,''''),findstr(ts,'&'),...
                   findstr(ts,'('),findstr(ts,')'),...
                   findstr(ts,'{'),findstr(ts,'}'),...
                   findstr(ts,'['),findstr(ts,']')...
                   ));
   for jj=1:length(extidx),
      ts=[ts(1:(extidx(jj)-1)),'\',ts(extidx(jj):end)];
      extidx((jj+1):end)=extidx((jj+1):end)+1;
   end
   
   if isempty(filedata(ii).playsec) || filedata(ii).playsec<1,
      fprintf('playsec empty for:\n%s\n',filedata(ii).file);
      [s,w]=unix(['mp3info -p "%S" ',ts]);
      
      if s==0 & w(1)~='/',
         playsec=str2num(w);
      else
         %keyboard
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
   if isempty(filedata(ii).bitrate) || ...
         ismember(filedata(ii).bitrate,[-2 0 1])
      fprintf('bitrate empty for:\n%s\n',filedata(ii).file);
      
      [s,w]=unix(['mp3info -p "%r" ',ts]);
      
      if s==0 && w(1)~='/',
         bitrate=str2num(w);
         if isempty(bitrate) && ~isempty(findstr(lower(w),'variable')),
            bitrate=-1;
         elseif isempty(bitrate),
            bitrate=0;
         end
      else
         bitrate=-2;
      end
      
      fprintf('setting bitrate=%d\n',bitrate);
      filedata(ii).bitrate=bitrate;
      %keyboard
      sql=['UPDATE dMusic SET ',...
           ' bitrate=',num2str(bitrate),...
           ' WHERE id=',num2str(filedata(ii).id)];
      mysql(sql);
   end
   if ~CHECKBR && (isempty(filedata(ii).year) || filedata(ii).year<1),
      fprintf('year empty for:\n%s\n',filedata(ii).file);
      [s,w]=unix(['mp3info -p "%y" ',ts]);
      
      if s==0 && ~isempty(w) && w(1)~='/',
         year=str2num(w);
      else
         year=0;
      end
      
      if ~isempty(year) && year>0 && ~isinf(year),
         fprintf('setting year=%d\n',year);
         sql=['UPDATE dMusic SET ',...
              ' year=',num2str(year),...
              ' WHERE id=',num2str(filedata(ii).id)];
         %keyboard
         mysql(sql);
      end
   end
   if CHECKBR && filedata(ii).bitrate>0,
      m4file=strrep(filedata(ii).file,'.mp3','.m4a');
      m4file=strrep(m4file,'.MP3','.m4a');
      
      dm4a=dir(m4file);
      if ~isempty(dm4a),
         bm4a=dm4a.bytes;
         dmp3=dir(filedata(ii).file);
         bmp3=dmp3.bytes;
         guessorigbr=round(filedata(ii).bitrate./bmp3.*bm4a);
         if guessorigbr<80,
            betterguess=56;
         elseif guessorigbr<110,
            betterguess=96;
         elseif guessorigbr<150,
            betterguess=128;
         elseif guessorigbr<190,
            betterguess=160;
         elseif guessorigbr<220,
            betterguess=196;
         elseif guessorigbr<245
            betterguess=224;
         elseif guessorigbr>280
            betterguess=320;
         else
            betterguess=256;
         end
         if betterguess>filedata(ii).bitrate,
            fprintf('%d: m4a for %s has higher br?\n',...
                    filedata(ii).id,basename(filedata(ii).file));
            sql=['UPDATE dMusic set cat1=',num2str(betterguess),...
                 ' WHERE id=',num2str(filedata(ii).id)];
            %keyboard
            mysql(sql);
         elseif betterguess<filedata(ii).bitrate,
            fprintf('%d: m4a for %s has lower br: %d from %d\n',...
                    filedata(ii).id,...
                    basename(filedata(ii).file),betterguess,...
                    filedata(ii).bitrate);
            sql=['UPDATE dMusic set bitrate=',num2str(betterguess),...
                 ', cat1=',num2str(betterguess),...
                 ' WHERE id=',num2str(filedata(ii).id)];
            %keyboard
            mysql(sql);
            
         else
            fprintf('%d: m4a for %s has matching br\n',...
                    filedata(ii).id,basename(filedata(ii).file));
         end
      end
   end
end

   

