% function [computerid,compdata]=dbgetcomp(machinename)
%
% created SVD 7/26/03
%
% machine name is 'name.ext....' (eg, lsd.jlg.berkeley.edu or lsd.jlg);
%
function [computerid,compdata]=dbgetcomp(machinename);

dbopen;

t=strsep(machinename,'.',1);
host=strtrim(t{1});

if length(t)>1,
   ext=t{2};
   sql=['SELECT * FROM tComputer where name="',host,'"',...
        ' AND ext="',ext,'"'];
else
   ext='%';
   sql=['SELECT * FROM tComputer where name="',host,'"'];
end

compdata=mysql(sql);


if length(compdata)==1,
   computerid=compdata(1).id;
elseif length(compdata)>1,
   disp('>1 matching computers, compidx assigned randomly');
else
   % computer doesn't exist in db, add it
   if strcmp(ext,'%'),
      ext='jlg';
   end
   computerid=-1;
   return
   
   dbcomputeredit(host,ext);
   
   % now get the info about the newly added computer
   sql=['SELECT * FROM tComputer where name="',host,'"',...
        ' AND ext like "',ext,'"'];
   compdata=mysql(sql);
   computerid=compdata(1).id;
   
end

