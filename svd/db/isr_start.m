% function isr_start(queueidx)

function isr_start(queueidx)

isrhome='/homes/svd/';

queuedata=dbgetqueue(queueidx);

if ~isempty(queuedata),
   fprintf('sending job %d to cluster\n',queueidx);
   
   localqueuefile=sprintf('/tmp/mat.%d.sh',queueidx);
   fid=fopen(localqueuefile,'w');
   fprintf(fid,'#PBS -lwalltime=32:00:00\n');
   fprintf(fid,'#PBS -lnodes=1:ppn=%d\n\n',queuedata.priority);
   fprintf(fid,'matlab -r "queuerun(%d)"\n',queueidx);
   fclose(fid);
   cmd=['scp ',localqueuefile,' svd@seil.umd.edu:',isrhome,'queue/'];
   unix(cmd);
   
   queuefile=sprintf('%squeue/mat.%d.sh',isrhome,queueidx);
   
   %cmd=['ssh svd@seil.umd.edu qsub -q narrow-long ',queuefile];
   cmd=['ssh svd@seil.umd.edu qsub ',queuefile];
   unix(cmd);
   
   sql=['UPDATE tQueue SET startdate=now(),complete=-1,progress=0,pid=0,machinename=""',...
        ' WHERE id=',num2str(queueidx)];
   mysql(sql);
   
   delete(localqueuefile);
   
else
   fprintf('no queue entry %d\n',queueidx);
end
