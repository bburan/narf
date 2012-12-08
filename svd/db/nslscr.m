
dbopen;

rawdata=mysql('SELECT * FROM gDataRaw');

for ii=1:length(rawdata),
   if isempty(rawdata(ii).parmfile),
      parmfile=[rawdata(ii).respfile ,' (guessing respfile)'];
   else
      parmfile=rawdata(ii).parmfile;
   end
   
   if findstr(rawdata(ii).resppath,'/afs/glue.umd.edu/'),
      fprintf('AFS: %s%s\n',rawdata(ii).resppath, ...
              parmfile);
   elseif length(rawdata(ii).resppath)>=7 &...
          strcmpi(rawdata(ii).resppath(1:7),'m:/daq/'),
      fprintf('metal: %s%s\n',rawdata(ii).resppath, ...
              parmfile);
   else
       fprintf('other: %s%s\n',rawdata(ii).resppath, ...
              parmfile);
   end
end