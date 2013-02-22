% function [masterid,rawid,penid,singleid,singlerawid]=dbgetcellid(cellid)
%
% masterid - index of site in gCellMaster
% rawid - index of all entries for this cell in gDataRaw
% penid - index of cell's gPenetration entry
% singleid - index of cell's gSingleCell entry
% singlerawid - index of cell's gSingleRaw entries... link to gDataRaw
%
function [masterid,rawid,penid,singleid,singlerawid]=dbgetcellid(cellid);

dbopen;

sql=['SELECT * FROM gSingleCell WHERE cellid="',cellid,'"'];
celldata=mysql(sql);
if length(celldata)==0,
   fprintf('error: cellid %s not found in gCellMaster\n',cellid);
   masterid=[];
   rawid=[];
   penid=[];
   singleid=[];
   return
end
masterid=celldata(1).masterid;
singleid=celldata(1).id;

sql=['SELECT id,rawid FROM gSingleRaw WHERE singleid=',num2str(celldata(1).id),];
rawdata=mysql(sql);
rawid=cat(1,rawdata.rawid);
singlerawid=cat(1,rawdata.id);

sql=['SELECT id FROM gPenetration WHERE id=',num2str(celldata(1).penid)];
pendata=mysql(sql);
penid=cat(1,pendata.id);

return

dbopen;

sql=['SELECT id,penid FROM gCellMaster WHERE cellid="',cellid,'"'];
celldata=mysql(sql);
if length(celldata)==0,
   fprintf('error: cellid %s not found in gCellMaster\n',cellid);
   masterid=[];
   rawid=[];
   penid=[];
   return
end
masterid=celldata(1).id;

sql=['SELECT id FROM gDataRaw WHERE cellid="',cellid,'"'];
rawdata=mysql(sql);
rawid=cat(1,rawdata.id);

sql=['SELECT id FROM gPenetration WHERE id=',num2str(celldata(1).penid)];
pendata=mysql(sql);
penid=cat(1,pendata.id);

sql=['SELECT id FROM gSingleCell WHERE cellid="',cellid,'"'];
singledata=mysql(sql);
singleid=cat(1,singledata.id);



