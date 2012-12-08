function logfig(figid,logfile)

if ~exist('figid','var'),
   figid=gcf;
end
if ~exist('logfile','var'),
   logfile=['/auto/k5/',getenv('USER'),'/log/daily-log-',num2str(year(now)),...
            '-',sprintf('%.2d',month(now)),...
            '-',sprintf('%.2d',day(now)),'.ps'];
end

[b,p]=basename(logfile);
if ~exist(p,'dir'),
   % log directory
   
   disp(['creating log directory ',p]);
   unix(['mkdir -p ',p]);
end

if ~exist(logfile,'file'),
   % log file doesn't exist yet.  just print to that file
   
   fprintf('creating log file %s\n',logfile);
   fprintf('logging figure %d\n',figid);
   tp=[tempname,'.ps'];  %figure
   print(['-f',num2str(figid)],tp);
   [w,s]=unix(['\mv ',tp,' ',logfile]);
else
   % log file already exists. print figure to a temp file and then
   % concatenate onto the end of the existing log file using a2ps
   
   fprintf('logging figure %d\n',figid);
   tp=[tempname,'.ps'];  %figure
   tp2=[tempname,'.ps']; %temp concatenated file
   print(['-f',num2str(figid)],tp);
   %['a2ps --columns=1 -o ',tp2,' ',logfile,' ',tp]
   [w,s]=unix(['a2ps --columns=1 -o ',tp2,' ',logfile,' ',tp]);
   %['\mv ',tp2,' ',logfile]
   [w,s]=unix(['\mv ',tp2,' ',logfile]);
   delete(tp);
end
   