% function r=dblog(string,fprintfargs)
function r=dblog(varargin)

persistent LASTDATE
LOGCOUNT=3;
outfile='/auto/data/tmp/queue/queuemasterlog.txt';
%outfile='/data/tmp/queue/queuemasterlog.txt';

if ~isempty(LASTDATE) & ~strcmp(LASTDATE,date),
   disp('rotating log files');
   delete([outfile,'.',num2str(LOGCOUNT)]);
   for ii=LOGCOUNT:-1:2,
      unix(['mv ',outfile,'.',num2str(ii-1),' ',outfile,'.',num2str(ii)]);
   end
   unix(['mv ',outfile,' ',outfile,'.1']);
   
   dbqueuetidy;
end

LASTDATE=date;

s=varargin{1};
if length(varargin)>1,
   fprintfargs={varargin{2:end}};
else
   fprintfargs={};
end

if not(exist(outfile,'file')),
   fid = fopen(outfile,'w');
   fprintf(fid,'%s - DB LOG FILE STARTED\n',datestr(now));
   fclose(fid);
   unix(['chmod a+w ',outfile]);
end

fid = fopen(outfile,'a');
fprintf(fid,['%s - ',s,'\n'],datestr(now),fprintfargs{:});
fclose(fid);

% display to screen too for the time being
fprintf(['%s - ',s,'\n'],datestr(now),fprintfargs{:});

