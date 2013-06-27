% function r=dbget(table,id,field)
%
% select table.field where table.id=id
% if field not specified, return all fields for table.id=id
% if id not specified, return entire table
%
% created SVD 2002
%
function r=dbget(table,id,field)

v=-1;  % value to return on failure
if not(exist('table','var')),
   disp('SYNTAX: dbget(table,id,field)');
   return
end

dbopen;

if exist('field','var') & ~isempty(field),
   sql=['SELECT ',field,' FROM ',table,' WHERE id=',num2str(id),';'];
   res=mysql(sql);
   r=getfield(res,field);

elseif exist('id','var'),
   sql=['SELECT * FROM ',table,' WHERE id=',num2str(id),';'];
   r=mysql(sql);
   
else
   sql=['SELECT * FROM ',table,';'];
   r=mysql(sql);
end




