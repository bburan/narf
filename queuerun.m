% queuerun.m:  wrapper for running a generic matlab program from
% the master queue
%
% basic flow: 
% 1. figure out current queue id from enviroment (set by dbqueuemaster)
% 2. load up parameters from tQueue
% 3. execute parameters in tQueue(QUEUEID).parmsstring
%
% created SVD 6/7/03 -- ripped off of cellxcqueue
%
function queuerun(queueid)

[s,host]=unix('hostname');

if exist('narf_set_path','file'),
    narf_set_path;
end
if exist('baphy_set_path','file'),
    baphy_set_path;
else
    addpath ~/code/baphy
    baphy_set_path
end

dbopen;

global BATQUEUEID
if exist('queueid','var'),
    BATQUEUEID=queueid;
else
    BATQUEUEID=str2num(getenv('QUEUEID'));
end
QUEUEOUTPATH=getenv('QUEUEOUTPATH');
if isempty(QUEUEOUTPATH),
    QUEUEOUTPATH='/tmp';
end

if isempty(BATQUEUEID),
   disp('syntax error: cellxcmaster(queueid) parameter required');
   disp('              or QUEUEID environment variable must be set');
   return
end

queuedata=dbgetqueue(BATQUEUEID);
if isempty(queuedata)
   fprintf('queue id=%d not found!\n',BATQUEUEID);
   return
end

if ~isempty(findstr(lower(host),'seil.umd.edu')),
   dbsetqueue(BATQUEUEID,1,-1);
end

% for debugging
fprintf('QUEUEID=%d\n',BATQUEUEID);
fprintf('MYHOST=%s\n',getenv('MYHOST'));

% run the matlab commands specified in parmstring. presumably this
% takes care of all the output.
disp(['RUNNING: ' char(queuedata.parmstring)]);
eval(char(queuedata.parmstring));

% record that we're done with this queue entry
% (may've been cleared, so requery environment)
BATQUEUEID=str2num(getenv('QUEUEID'));

if isempty(findstr(lower(host),'seil.umd.edu')),
    s=get(0,'Children');
    s=s(s<=3);
    figpath=sprintf('%s/%d',QUEUEOUTPATH,floor(BATQUEUEID/1000)*1000);
    if ~exist(figpath,'dir'),
        unix(['mkdir ',figpath]);
    end
    
    for ii=s(:)',
        figure(ii);
        fullpage portrait
        print(['-f' num2str(ii)],'-djpeg','-r75',sprintf('%s/%d.%d.jpg',figpath,BATQUEUEID,ii));
    end
end

global DB_SERVER BAPHY_LAB
if strcmpi(BAPHY_LAB,'lbhb'),
    DB_SERVER='hyrax.ohsu.edu';
    dbopen(1);
end

dbsetqueue(BATQUEUEID,1000,1);


