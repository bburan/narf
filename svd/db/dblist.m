% function dblist(cellid,batch)
%
% list files in scellfile for cell that match critera for batch
%
function dblist(cellid,batch)


[cellfiledata,times,batchdata]=cellfiletimes(cellid,batch);

fprintf('Cell %s Batch %d:\n',cellid,batch);

for fidx=1:length(cellfiledata),
   fprintf('stim: %s%s\nresp: %s%s\nfmt:  %s / %s\n',...
           cellfiledata(fidx).stimpath,cellfiledata(fidx).stimfile,...
           cellfiledata(fidx).path,cellfiledata(fidx).respfile,...
           cellfiledata(fidx).stimfilefmt,cellfiledata(fidx).respfilefmt);
end

