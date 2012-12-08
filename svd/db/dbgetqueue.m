% function queuedata=dbgetqueue(queueid)
%
% created SVD 4/1/03
%
function queuedata=dbgetqueue(queueid)

dbopen;
sql=['SELECT *,((TO_DAYS(NOW())-TO_DAYS(lastdate))*86400+',...
     ' TIME_TO_SEC(NOW())-TIME_TO_SEC(lastdate)) AS secago',...
     ' FROM tQueue WHERE id=',num2str(queueid)];
queuedata=mysql(sql);

