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

if ~exist('mysql','file'),
    addpath ~/code/baphy
    baphy_set_path
    
    %addpath ~svd/code/db/
    %addpath ~svd/code/mutils/
    %addpath ~svd/code/gen/
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
   figure(1);
   fullpage portrait
   print('-f1','-djpeg',sprintf('%s/%d.1.jpg',QUEUEOUTPATH,BATQUEUEID));
   figure(2);
   fullpage portrait
   print('-f2','-djpeg',sprintf('%s/%d.2.jpg',QUEUEOUTPATH,BATQUEUEID));
   figure(3);
   fullpage portrait
   print('-f3','-djpeg',sprintf('%s/%d.3.jpg',QUEUEOUTPATH,BATQUEUEID));
end

dbsetqueue(BATQUEUEID,1000,1);


