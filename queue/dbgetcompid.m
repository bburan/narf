% function computerid=dbgetcompid(machinename)
%
% machine name is 'name.ext....' (eg, lsd.jlg.berkeley.edu or lsd.jlg);
% if machinename not specified, use unix('hostname'). ie, the name
% of the machine where dbgetcompid is called.
%
% created SVD 7/26/03
%
function computerid=dbgetcompid(machinename);

if ~exist('machinename','var'),
   machinename=getenv('MYHOST');
   if isempty(machinename),
      machinename=fileread('/etc/hostname');
      machinename=strsep(machinename,char(10));
      machinename=machinename{1};
   end
   
   machinename=strtrim(deblank(machinename));
end

t=strsep(machinename,'.');
host=t{1};
if length(t)>1,
   ext=t{2};
else
   ext='%';
end

dbopen;
sql=['SELECT * FROM tComputer where name="',host,'"',...
     ' AND ext like "',ext,'"'];
compdata=mysql(sql);
if isempty(compdata),
   warning('computer not found in celldb!');
   computerid=-1;
   return
end

computerid=compdata(1).id;














