% function dbsetsite(siteid,rcunit,rcchan,[singleids],[rawids])
%
% save sorted identities of single units across data files for a
% recording site to celldb
%
% rcunit, 
% rcchan - R x C matrices, where R is the number of raw data files
%          and C is the number of cells.  values of rcunit and
%          rcchan indicate the unit and channel (electrode),
%          respectively for each (rawfile,cell) pair
%          rcunit() values of -1 or nan mean the cell did not exist
%          in that rawfile (marked as "crap" in db)
%
% rawids and singleids are vectors that point to entries in
% gDataRaw and gSingleCell, respectively. if they aren't passed,
% it's assumed that they are the same order as would be returned by
% dbgetsite.
%
% singleids==-1 means replace all existing cells with new data
%
% CREATED SVD 2005-09-02
%
function dbsetsite(siteid,rcunit,rcchan,singleids,rawids)

dbopen;
global DBUSER

sql=['SELECT * FROM gCellMaster WHERE siteid="',siteid,'"'];
site=mysql(sql);

if length(site)==0,
   fprintf('%s: site %s doesn''t exist in celldb\n',...
           mfilename,siteid);
   return
end

if ~exist('rawids','var'),
   sql=['SELECT * FROM gDataRaw WHERE masterid=',num2str(site.id),...
        ' ORDER BY parmfile,respfile,id'];
   rawdata=mysql(sql);
   rawids=cat(1,rawdata.id);
end

if length(rawids)~=size(rcunit,1),
   fprintf('%s: mismatch in rcunit size and rawid count\n',mfilename);
   return
end

if ~exist('singleids','var'),
   sql=['SELECT * FROM gSingleCell WHERE masterid=',num2str(site.id),...
        ' ORDER BY cellid'];
   celldata=mysql(sql);
   singleids=cat(1,celldata.id);
else
   if length(singleids)==1 & singleids==-1,
      % -1 means new single cell
      singleids=-ones(size(rcunit,2),1);
   end
end


if length(singleids)~=size(rcunit,2),
   fprintf('%s: mismatch in rcunit size and cell count\n',mfilename);
   return
end

% ok, a couple possibilities here:
% 1. all single cell infrastructure exists, just update r-c info
% 2. need to create new r-c entr(ies)

rawcount=length(rawids);
cellcount=length(singleids);

% confirm existence of each single cell. create ones for -1's if
% necessary
cellids={};
chanstr={'a','b','c','d','e'};
for ii=1:cellcount,
   if ~isnan(rcunit(1,ii))
      if max(rcchan)==1,
         cellids{ii}=sprintf('%s-%d',site.siteid,rcunit(1,ii));
      else
         cellids{ii}=sprintf('%s-%s%d',siteid,chanstr{rcchan(1,ii)},rcunit(1,ii));
      end
      
      if singleids(ii)==-1,
         % create new cell
         [aff,newidx]=sqlinsert('gSingleCell',...
                                'siteid',site.siteid,...
                                'cellid',cellids{ii},...
                                'masterid',site.id,...
                                'penid',site.penid,...
                                'unit',rcunit(1,ii),...
                                'channum',rcchan(1,ii),...
                                'addedby',DBUSER,...
                                'info','dbsetsite.m');
         if aff,
            fprintf('created cell %s (%d)\n',cellids{ii},newidx);
            singleids(ii)=newidx;
         end
      else,
         sql=['UPDATE gSingleCell set',...
              ' cellid="',cellids{ii},'",',...
              ' unit=',num2str(rcunit(1,ii)),',',...
              ' channum=',num2str(rcchan(1,ii)),...
              ' WHERE id=',num2str(singleids(ii))];
         [r,aff]=mysql(sql);
         if aff,
            fprintf('updated cell %s (%d)\n',cellids{ii},singleids(ii));
         end
      end
      
      for rr=1:rawcount,
         sql=['SELECT * FROM gSingleRaw',...
              ' WHERE singleid=',num2str(singleids(ii)),...
              ' AND rawid=',num2str(rawids(rr))];
         srdata=mysql(sql);
         
         % figure out if this unit exists for this file (mark crap if
         % doesn't exist)
         if isnan(rcunit(rr,ii)) | rcunit(rr,ii)==-1,
            crap=1;
            rcunit(rr,ii)=0;
            rcchan(rr,ii)=0;
         else
            crap=0;
         end
         
         if length(srdata)==0,
            % doesn't exist, add it
            [aff,newidx]=sqlinsert('gSingleRaw',...
                                   'cellid',cellids{ii},...
                                   'singleid',singleids(ii),...
                                   'masterid',site.id,...
                                   'penid',site.penid,...
                                   'rawid',rawids(rr),...
                                   'unit',rcunit(rr,ii),...
                                   'crap',crap,...
                                   'channum',rcchan(rr,ii),...
                                   'addedby',DBUSER,...
                                   'info','dbsetsite.m');
            if aff,
               fprintf('created cell %s link to rawid %d\n',...
                       cellids{ii},rawids(rr));
            end
         else
            sql=['UPDATE gSingleRaw set',...
                 ' cellid="',cellids{ii},'",',...
                 ' unit=',num2str(rcunit(rr,ii)),',',...
                 ' channum=',num2str(rcchan(rr,ii)),',',...
                 ' crap=',num2str(crap),...
                 ' WHERE id=',num2str(srdata.id)];
            [r,aff]=mysql(sql);
            if aff,
               fprintf('updated cell %s link to rawid %d\n',...
                       cellids{ii},rawids(rr));
            end
         end
      end
   end
end

% string with ids of singleids to preserve
singleidstr='(-1,';
for ii=singleids(singleids(:)>-1)',
   singleidstr=[singleidstr,num2str(ii),','];
end
singleidstr(end)=')';

sql=['SELECT * FROM gSingleCell',...
     ' WHERE masterid=',num2str(site.id),...
     ' AND not(id in ',singleidstr,')'];
exdata=mysql(sql);
if length(exdata)>0,
   fprintf('%s: WARNING! deleting %d cells that were not specified!\n',...
           mfilename,length(exdata));
   disp('press a key to go ahead');
   pause
   
   sql=['DELETE FROM gSingleRaw',...
     ' WHERE masterid=',num2str(site.id),...
     ' AND not(singleid in ',singleidstr,')'];
   mysql(sql);
   sql=['DELETE FROM gSingleCell',...
     ' WHERE masterid=',num2str(site.id),...
     ' AND not(id in ',singleidstr,')'];
   mysql(sql);
end


   
