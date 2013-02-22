% function r=dbset(table,id,field,value)
% 
% set table.field=value where table.id=id
%
% returns number of affected rows.
% 
% created SVD 2002
%
function r=dbset(table,id,field,value)

r=0;
if not(exist('value','var')),
   disp('SYNTAX: dbset(table,id,field,value)');
   return
end

dbopen;

if isnumeric(value),
   sql=['UPDATE ',table,' SET ',field,'=',num2str(value),' WHERE id=',num2str(id),';'];
else
   sql=['UPDATE ',table,' SET ',field,'=''',value,''' WHERE id=',num2str(id),';'];
end

[res,r]=mysql(sql);

