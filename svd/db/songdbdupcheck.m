
dbopen;
sql=['SELECT * FROM dMusic WHERE dup>1'];
dupdata=mysql(sql);

for ii=1:length(dupdata),
   sql=['SELECT * FROM dMusic WHERE id=',num2str(dupdata(ii).dup),...
        ' OR dup=',num2str(dupdata(ii).dup),' ORDER BY bitrate'];
   matchdata=mysql(sql);
   if mod(ii,100)==0,
      fprintf('.');
   end
   br=zeros(size(matchdata));
   dup=zeros(size(matchdata));
   playsec=zeros(size(matchdata));
   brguess=zeros(size(matchdata));
   wasm4a=zeros(size(matchdata));
   for jj=1:length(br),
      if isempty(matchdata(jj).bitrate),
         br(jj)=-2;
      else
         br(jj)=matchdata(jj).bitrate;
      end
      dup(jj)=matchdata(jj).dup;
      playsec(jj)=matchdata(jj).playsec;
      
     
      mp3checkfile=matchdata(jj).file;
      m4acheckfile=strrep(mp3checkfile,'.mp3','.m4a');
      %%d3=dir(mp3checkfile);
      %d4=dir(m4acheckfile);
      if exist(m4acheckfile,'file'),
         [w,s]=unix(['du -sk "',m4acheckfile,'"']);
         wasm4a(jj)=1;
      else
         [w,s]=unix(['du -sk "',mp3checkfile,'"']);
      end
      if ~w,
         s=strsep(s,9);
         s=s{1};
         brguess(jj)=s./playsec(jj).*8 .*1.1;
      else
         brguess(jj)=-2;
      end
   end
   
   if find(brguess<br),
      mm=find(brguess<br);
      if wasm4a(mm),
         disp('found mis-matched bitrate from m4a conversion!');
         likelybr=[56 96 128 160 192 256 320];
         gg=max(find(likelybr<brguess(mm)));
         fprintf('likely br=%d\n',likelybr(gg));
         sql=['UPDATE dMusic set bitrate=',num2str(likelybr(gg)),...
              ' WHERE id=',num2str(matchdata(mm).id)]
         disp('--pause--');
         pause;
         mysql(sql);
         br(mm)=likelybr(gg);
      else
         disp('found mis-matched bitrate');
         keyboard
      end
   end
   masterid=find(~dup);
   if br(masterid)<max(br),
      disp('found lower than max br');
      keyboard
   end
   
      
   
end

return   
   mb=min(find(br==max(br)));
   
   
   if br(~dup)>=br(mb),
      % check if really higher BR or artefact of m4a conversion
      
   else
      % if dup'ed actually is higher, then swap dup markings
   end  
   