% function r=dbrunrevert(runidx,changeby)
function r=dbrunrevert(runidx,changeby)

dbopen;
if ~exist('changeby'),
   changeby=-1;
end

for rr=1:length(runidx),
   
   sql=['SELECT * FROM tRunData WHERE id=',num2str(runidx(rr))];
   rundata=mysql(sql);
   
   if length(rundata)==0,
      fprintf('ERROR: runidx=%d not listed in tRunData.\n',runidx(rr));
   else
      resfileold=rundata.resfile;
      
      i1=max(findstr(resfileold,'.res.'));
      i2=max(findstr(resfileold,'.mat'));
      if ~isempty(i1),
         i1=i1+5;
         
         oldcount=str2num(resfileold(i1:i2-1));
         
         newcount=oldcount+changeby;
         resfilenew=[resfileold(1:i1-1) num2str(newcount) resfileold(i2:end)];
         
         sql=['UPDATE tRunData set resfile="',resfilenew,'"'...
              ' WHERE id=',num2str(runidx(rr))];
         mysql(sql);
         fprintf('%d: Changed %s to %s.\n',runidx(rr),resfileold,resfilenew);
      end
   end
end

