function dbprunebatch(batchid)

sql=['select * from sBatch where id=',num2str(batchid)];
batchdata=mysql(sql);

fprintf('You have selected to delete batch %d (%s)!\n',...
        batchid,batchdata.name);

yn=input('\nContinue (y/[n])? ','s');
if strcmp(yn,'y'),  % yes, go ahead and delete run
   disp('ok, deleting!');
   sql=['DELETE FROM sResults WHERE batch=',num2str(batchid)];
   mysql(sql);
   sql=['DELETE FROM sRunData WHERE batch=',num2str(batchid)];
   mysql(sql);
   sql=['DELETE FROM sBatch WHERE id=',num2str(batchid)];
   mysql(sql);
end
