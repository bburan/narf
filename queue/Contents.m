% Stephen's matlab/mysql programs
% Last update: 2/11/03
%
% Basic database i/o:
%   dbopen - open a connection to the cell database
%   dbget - read a value for a row in a field of a table
%   dbset - set the value for a row in a field of a table
%
% Queue operations:
%   dbqueue - get information about entries currently in the queue
%   dbaddqueue - create a new queue entry
%   dbsetqueue - change the status of a queue entry
%   dbgetnextqueue - find next available queue entry and mark as
%                    actively being processed (complete=-1)
%   dbdeletequeue - remove a queue entry
%   dbkillqueue - mark an active queue entry for cancellation
%   dbcleanqueue - scan for dead entries and mark as unprocessed
%   dbloadmon - graphical report on queue load and completion status
%   dbgetload - non-graphical queue status
%   dbupdateload - update load information for an active computer
%                  processing queue data (called by dbgetload)
%

