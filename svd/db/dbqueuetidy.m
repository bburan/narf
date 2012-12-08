function [totaljobs,shiftedjobs]=dbqueuetidy;

dbopen;

maxdays=28;

% delete dead and completed jobs that haven't been changed for 28 days
sql=['SELECT *,(TO_DAYS(NOW())-TO_DAYS(lastdate)) AS daysago',...
     ' FROM tQueue',...
     ' WHERE (TO_DAYS(NOW())-TO_DAYS(lastdate))>',num2str(maxdays),...
     ' AND allowqueuemaster',...
     ' AND complete in (1,2)',...
     ' ORDER BY id'];
oqueuedata=mysql(sql);
deletedjobs=length(oqueuedata);

fprintf('dbqueuetidy.m: deleting %d jobs (%d days idle)\n',...
        deletedjobs,maxdays);

for ii=1:deletedjobs,
   sql=['DELETE FROM tQueue WHERE id=',num2str(oqueuedata(ii).id)];
   mysql(sql);
end


% find lowest id of an active or not-started job
sql=['SELECT count(id) as totaljobs,(max(id)+1) as maxid FROM tQueue'];
alldata=mysql(sql);
totaljobs=alldata.totaljobs;
fprintf('%d total job count\n',totaljobs);

disp('future possibility: shift queue ids?');
return
disp('need to figure out how to deal with logs!');

% find lowest id of an active or not-started job
sql=['SELECT min(id) as minid FROM tQueue WHERE complete in (-1,0)'];
lodata=mysql(sql);
lid=lodata.minid;

% none found, just find the max id in tQueue
if length(lid)==0,
   lid=alldata.maxid;
end

% find all entries with id < lid
sql=['SELECT * FROM tQueue WHERE allowqueuemaster',...
     ' AND complete in (1,2) AND id<',num2str(lid),...
     ' ORDER BY id'];
lodata=mysql(sql);

% shift numbers down
shiftedjobs=0;
for cminid=1:length(lodata),
   if lodata(cminid).id>cminid,
      sql=['UPDATE tQueue SET id=',num2str(cminid),...
           ' WHERE id=',num2str(lodata(cminid).id)];
      mysql(sql);
      shiftedjobs=shiftedjobs+1;
   end
end
fprintf('%d jobs shifted to lower ids\n',totaljobs);

