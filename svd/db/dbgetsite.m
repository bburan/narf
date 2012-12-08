% function [rawdata,site,celldata,rcunit,rcchan]=dbgetsite(siteid)
%
% siteid is string name of site. 
% 
% returns: all entries from gDataRaw (rawdata), gCellMaster (site)
%           and gSingleCell (celldata) that match siteid=<siteid>
%          rcunit and rcchan- unit number (sorted on the electrode)
%           and channel number (electrode) for each raw cell and raw
%           file
%
% created SVD 2005-09-01 -- ripped off dbgetscellfile
%
function [rawdata,site,celldata,rcunit,rcchan]=dbgetsite(siteid)

dbopen;

sql=['SELECT * FROM gCellMaster WHERE siteid="',siteid,'"'];
site=mysql(sql);

if length(site)>0,
   sql=['SELECT * FROM gDataRaw WHERE masterid=',num2str(site.id),...
        ' AND not(bad) ORDER BY parmfile,respfile,id'];
   rawdata=mysql(sql);
   rawids=cat(1,rawdata.id);
   
   sql=['SELECT * FROM gSingleCell WHERE masterid=',num2str(site.id),...
        ' ORDER BY cellid'];
   celldata=mysql(sql);
   cellids=cat(1,celldata.id);
   
else
   fprintf('%s: site %s doesn''t exist in celldb\n',...
           mfilename,siteid);
   rawdata=[];
   celldata=[];
end

rawcount=length(rawdata);
cellcount=length(celldata);

rcunit=zeros(rawcount,cellcount).*nan;
rcchan=zeros(rawcount,cellcount).*nan;

for ii=1:rawcount,
   sql=['SELECT * FROM gSingleRaw WHERE rawid=', ...
        num2str(rawids(ii))];
   singrawdata=mysql(sql);
   for jj=1:length(singrawdata),
      jji=find(singrawdata(jj).singleid==cellids);
      rcunit(ii,jji)=singrawdata(jj).unit;
      rcchan(ii,jji)=singrawdata(jj).channum;
   end
end
