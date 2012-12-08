

dbopen;
sql=['SELECT DISTINCT cellid from sRunData',...
     ' WHERE batch in (86)',...
     ' ORDER BY cellid'];
rundata=mysql(sql);


for ii=length(rundata):-1:2, % :length(rundata)-1,
   v1cellcomp(rundata(ii).cellid,[86 83 83 86 87],0,0,[4 1 2 2 4]);
   drawnow
   %print -f1 -dpsc2 -Plj2200
   %print -f2 -dps2 -Plj2200
   %print -f2 -dps2 -Pgcolor
   print -f3 -dps2 -Plj2200
   %print -f1 -dpsc2 -Pgcolor
end


dbopen;
sql=['SELECT DISTINCT cellid from sRunData',...
     ' WHERE batch in (24)',...
     ' ORDER BY cellid'];
rundata=mysql(sql);


for ii=11:-1:1, % :length(rundata)-1,
   v1cellcomp(rundata(ii).cellid,[24,27,30]);
   drawnow
   %print -f1 -dpsc2 -Plj2200
   print -f2 -dpsc2 -Plj2200
end
