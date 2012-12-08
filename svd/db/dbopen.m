% function [r,errmsg]=dbopen(dbserver,dbuser,dbpassword,dbname,force)
%
% open connection to mysql database (requires mysql mex file)
%
% force - optional.  if force==1, reopen connection, even if global variable
%         DBISOPEN is already non-zero (allows connecting to a
%         different database)
%
% created SVD 6/2001
% modified - SVD - 2005-10-10 - persistent, optional db parameters
%
function [r,errmsg]=dbopen(dbserver,dbuser,dbpassword,dbname,force)

global DBISOPEN DBUSER
global DB_USER DB_SERVER DB_PASSWORD DB_NAME
global USEDB

if ~USEDB,
   r=1;
   return;
end

% force backwards compatibility
if exist('dbserver','var'),
   if isnumeric(dbserver),
      force=dbserver;
      clear dbserver
   end
end

% figure out platform-specific settings for mysql client program
dbcodepath=fileparts(which(mfilename));
if strcmp(computer,'PCWIN'),
    addpath([dbcodepath filesep 'db_win']);
    global MYSQL_BIN_PATH
    MYSQL_BIN_PATH=[dbcodepath filesep 'db_win' filesep];
    global LOCAL_DATA_ROOT
    LOCAL_DATA_ROOT = 'N:\';
else
    addpath([dbcodepath filesep 'db_linux']);
end

% set defaults on first call
if isempty(DB_SERVER), DB_SERVER='metal.isr.umd.edu'; end
if isempty(DB_USER), DB_USER='david'; end
if isempty(DB_PASSWORD), DB_PASSWORD='nine1997'; end
if isempty(DB_NAME), DB_NAME='cell'; end

DBUSER=getenv('USER');
if isempty(DBUSER), DBUSER=getenv('user'); end
if isempty(DBUSER), DBUSER=DB_USER; end

% use remembered values if connection parameters not passed
if ~exist('dbserver','var'), dbserver=DB_SERVER; end
if ~exist('dbuser','var'), dbuser=DB_USER; end
if ~exist('dbpassword','var'), dbpassword=DB_PASSWORD; end
if ~exist('dbname','var'), dbname=DB_NAME; end
if ~exist('force','var'),
   force=0;
end

% save as defaults if parameters not passed next time
DB_SERVER=dbserver;
DB_USER=dbuser;
DB_PASSWORD=dbpassword;
DB_NAME=dbname;

if ~force,
   try
      mysql(['use ',dbname]);
      DBISOPEN=1;
      errmsg='';
   catch
      force=1;
      DBISOPEN=0;
      errmsg=lasterr;
      r=0;
      if strcmp(computer,'PCWIN'),
         % this is the only way to open on the windows system, so
         % just return error.
         return;
      end
   end
end

% connection wasn't open.  try to open one.
if isempty(DBISOPEN) | not(DBISOPEN) | force==1,
    try,
        hostname=getenv('HOSTNAME');
        if strcmp(hostname,dbserver),
            mysql('open','localhost',dbuser,dbpassword);
        else
            mysql('open',dbserver,dbuser,dbpassword);
        end

        mysql(['use ',dbname]);  % ie, cell
        DBISOPEN=1;
    catch
        DBISOPEN=0;
        errmsg=lasterr;
    end
end

r=DBISOPEN;





