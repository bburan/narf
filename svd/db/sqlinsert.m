function [affected,insertid]=sqlinsert (table,varargin);
%SQLINSERT - insert records into an SQL database
%
%    [AFFECTED,INSERTID] = SQLINSERT (TABLE,'valuename',value,...);
%
%    See also: MYSQL

sqlc=['INSERT INTO ' table '('];
for ii=1:2:length(varargin)-2,
   sqlc=[sqlc varargin{ii} ','];
end;
ii=length(varargin)-1;
sqlc=[sqlc varargin{ii} ') VALUES ('];
for ii=2:2:length(varargin)-2,
   if ischar(varargin{ii}),
      sqlc=[sqlc '"' varargin{ii} '"' ','];
   else
      sqlc=[sqlc num2str(varargin{ii}) ','];
   end;
end;
ii=length(varargin);
if ischar(varargin{ii}),
   sqlc=[sqlc '"' varargin{ii} '"' ')'];
else
   sqlc=[sqlc num2str(varargin{ii},20) ')'];
end;
[res,affected,insertid]=mysql(sqlc);
