% function d=rfdiam(cellid,[pixref]);
%
% get cell rf diameter in degrees. if pixref is given, return the
% degrees corresponding to that pixel value for this cell
%
function d=rfdiam(cellid,pixref);

dbopen;

sql=['SELECT sCellFile.*,gPenetration.etudeg,gCellMaster.rfppd',...
     ' FROM sCellFile,gCellMaster,gPenetration',...
     ' WHERE gCellMaster.id=sCellFile.masterid',...
     ' AND gPenetration.id=gCellMaster.penid',...
     ' AND stimwindowsize>0',...
     ' AND stimfilecrf>0',...
     ' AND sCellFile.cellid="',cellid,'"'];
cellfiledata=mysql(sql);

if length(cellfiledata)>0,
   ii=1;
   
   if ~exist('pixref','var'),
      pixref=cellfiledata(ii).stimwindowsize./cellfiledata(ii).stimfilecrf;
   end
   if isempty(cellfiledata(ii).etudeg) | cellfiledata(ii).etudeg==0,
      d=pixref./cellfiledata(ii).rfppd;
   else
      d=pixref./cellfiledata(ii).etudeg;
   end

else
   error('no data for this cell in sCellFile with valid stimwindowsize ');
end
