function israddruns(runids)

dbopen;
for ii=1:length(runids),
   sql=['SELECT * FROM sRunData WHERE id=',num2str(runids(ii))];
   rundata=mysql(sql);
   if length(rundata)>0,
      isrmakequeue(rundata.cellid,rundata.batch);
   end
end
