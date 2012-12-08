function pid=getpid()

pid=getenv('MYPID');
if strcmp(pid,''),
   pid=getenv('PBS_JOBID');
   if ~isempty(pid),
      pid=strsep(pid,'.',1);
      pid=pid{1};
   end
end

if isempty(pid),
   pid=0;
else
   pid=str2num(pid);
end

