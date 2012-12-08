%function dbchangerootpath(oldroot,newroot,table,field);
%
function dbchangerootpath(oldroot,newroot,table,field);

dbopen;

sql=['SELECT id,',field,' as oldpath, concat("',newroot,'",',...
     'substring(',field,',',num2str(length(oldroot)+1),')) as newpath',...
     ' FROM ',table,' WHERE ',field,' like "',oldroot,'%"'];
pathdata=mysql(sql);

if length(pathdata)==0,
   disp('No matches.');
   return;
end

for ii=1:length(pathdata),
   fprintf('%d: %s --> %s\n',pathdata(ii).id,pathdata(ii).oldpath,...
           pathdata(ii).newpath);
end

yn=input('Continue with update (y/[n])??? ','s');
if length(yn)==0 | ~strcmp(yn(1),'y'),
   return;
end

sql=['UPDATE ',table,' SET ',field,'=concat("',newroot,'",',...
     'substring(',field,',',num2str(length(oldroot)+1),'))',...
     ' WHERE ',field,' like "',oldroot,'%"'];
mysql(sql);

