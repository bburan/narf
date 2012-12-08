% function v=dbdelrec(stable,value,sfield[='id'])
function v=dbdelrec(stable,value,sfield,force)

v=0;

if ~exist('value','var'),
   disp('SYNTAX: dbdelrec(stable,value,sfield[=''id''])');
   return
end
if ~exist('force','var'),
   force=0;
end

if ~exist('sfield','var'),
   sfield='id';
end

if isnumeric(value),
   sval=num2str(value);
else
   sval=['"',value,'"'];
end

dbopen;
sql=['SELECT * FROM ',stable,' WHERE ',sfield,'=',sval];
data=mysql(sql);

if length(data)==0,
   fprintf('%s=%s not found in table %s.\n',sfield,sval,stable);
   return
end

if force
   yn='y';
else
   data
   yn=input('Continue with delete (y/[n])? ','s');
end

if length(yn)>0 & strcmp(yn(1),'y'),
   sql=['DELETE FROM ',stable,' WHERE ',sfield,'=',sval];
   [res,v]=mysql(sql);
end


