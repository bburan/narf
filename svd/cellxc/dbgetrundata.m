% function [rundata,batchdata]=dbgetrundata(batchid)
%
% find all entries in sRunData with batch=batchid
% rundata -  structure array, one entry per data set
% batchdata - parameters used for this batch
%
% created SVD 10/31/04 - hacked from kvaparms.m
%
function [rundata,batchdata]=dbgetrundata(batchid);

dbopen;
sql=['SELECT * FROM sRunData',...
     ' WHERE batch=',num2str(batchid),...
     ' AND not(cellid like "m0000")',...
     ' AND not(cellid like "model%")',...
     ' ORDER BY cellid'];
rundata=mysql(sql);

if length(rundata)==0,
   fprintf('no runs found for batch %d.\n',batchid);
end

batchdata=dbget('sBatch',batchid);

