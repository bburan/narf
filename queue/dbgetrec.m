% function r=dbgetrec(stable,value,sfield[='id'])
function r=dbgetrec(stable,value,sfield)

v=0;

if ~exist('value','var'),
   disp('SYNTAX: dbdelrec(stable,value,sfield[=''id''])');
   return
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
r=mysql(sql);



