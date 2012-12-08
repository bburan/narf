dbopen;

sql='select distinct respfile,parmfile from gDataRaw where not(training) and respfile like "/auto/%"';

pathdata=mysql(sql);


for ii=2001:length(pathdata),
   testpath=[pathdata(ii).respfile];
   [pp,parmfile,ee]=fileparts(pathdata(ii).parmfile);
   testfile=[testpath parmfile '001.map'];
   if exist(testfile,'file'),
      %ls([testpath parmfile '*'])
      disp(['cd ',testpath]);
      cd(testpath)
      disp(['tar czvf ' parmfile '.tgz ' parmfile '*map --remove-files'])
      unix(['tar czvf ' parmfile '.tgz ' parmfile '*map --remove-files'])
   end
end
