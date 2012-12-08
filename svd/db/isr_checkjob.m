function [status,extended]=isr_checkjob(queueidx);

queuedata=dbgetqueue(queueidx);

extended=[];
if ~isempty(queuedata),

  
   cmd=['ssh svd@seil.umd.edu checkjob -A ',num2str(queuedata.pid)];
   [w,s]=unix(cmd);
   status=1-w;
else
   fprintf('no queue entry %d\n',queueidx);
   status=0;
end

if status>0,
   s=strsep(strtrim(s),';');
   for ii=1:length(s),
      ts=strsep(s{ii},'=');
      if length(ts)>1,
         extended=setfield(extended,ts{1},ts{2});
      elseif length(ts)>0,
         extended=setfield(extended,ts{1},[]);
      end
   end
end