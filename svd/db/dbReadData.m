% function [parm,perf]=dbReadData(rawid);
%
% created SVD 2006-02-23
%
function [parm,perf]=dbReadData(rawid);

if ~exist('rawid','var'),
   error('parameter rawid required');
end

dbopen;
global DB_USER

sql=['SELECT * FROM gData WHERE rawid=',num2str(rawid),...
     ' ORDER BY id'];
data=mysql(sql);

parm=[];
perf=[];

for ii=1:length(data),
   switch data(ii).datatype,
    case 0,
     val=data(ii).value;
    case 1,
     val=eval(data(ii).svalue);
    case 2,
     val=data(ii).svalue;
   end
   
   if data(ii).parmtype==0,
      parm=setfield(parm,data(ii).name,val);
   else
      perf=setfield(perf,data(ii).name,val);
   end
end
