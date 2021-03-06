% function [cellfileids]=dms2scellfile(params)
%
% given a bunch of information about imsm/resp.dat files produced from
% dms pre-processing, insert a new entry into sCellFile. assumption
% is that a for any gDataRaw entry, there should only be a single
% entry in sCellFile per respfmtcode.
%
% params.respfmtcode [default 0] psth 
%                             1  pfth
% params.nosync [default 0] if 1, no sync correction made for spike times
%
% CREATED SVD 8/1/03  (ripped off dbraw2scellfile.m)
% modified SVD 6/29/04 (made more generic. got rid of a_state)
% modified SVD 12/20/04 (added support for mult simul neurons)
%
function cellfileid=dms2scellfile(params)

dbopen;
rcsetstrings;

% this assumes non-fixation triggered data.
% if it's PFTH, respfmtcode should be set to 1
params.respfmtcode=getparm(params,'respfmtcode',0);
params.nosync=getparm(params,'nosync',0);

sql=['SELECT * FROM sCellFile',...
     ' WHERE singlerawid=',num2str(params.singlerawid),...
     ' AND respfmtcode=',num2str(params.respfmtcode)];
cellfiledata=mysql(sql);

sql=['SELECT gSingleRaw.*,gDataRaw.runclassid,gDataRaw.stimspeedid',...
     ' FROM gSingleRaw',...
     ' INNER JOIN gDataRaw ON gSingleRaw.rawid=gDataRaw.id',...
     ' WHERE gSingleRaw.id=',num2str(params.singlerawid)];
rawdata=mysql(sql);

sql=['SELECT * FROM gSingleCell',...
     ' WHERE id=',num2str(params.singleid)];
celldata=mysql(sql);

if celldata.rfsize==0,
   disp('error!  rfsize cannot be zero!  dms2scellfile aborted.');
   return
end

if length(cellfiledata)>0,
   fprintf('sCellFile %d entr(ies) exist for singlerawid=%d. o-writing.\n',...
           length(cellfiledata),params.singlerawid);
   for ii=1:length(cellfiledata),
      sql=['DELETE FROM sCellFile WHERE id=',num2str(cellfiledata(ii).id)];
      mysql(sql);
   end
end

respfilefmt=respfilefmtstr{params.respfmtcode+1};
stimfmtcode=2;
stimfilefmt=stimfilefmtstr{stimfmtcode+1};
[framecount,iconside,stimfile]=imfileinfo([params.stimpath,params.stimfile],1);
iconside=sprintf('%d,%d',iconside(1),iconside(2));

stimfilecrf=params.stimwindowsize./celldata.rfsize;

[aff,cellfileid]=...
    sqlinsert('sCellFile',...
              'cellid',celldata(1).cellid,...
              'masterid',rawdata.masterid,...
              'rawid',params.rawid,...
              'singleid',params.singleid,...
              'singlerawid',params.singlerawid,...
              'runclassid',rawdata.runclassid,...
              'path',params.path,...
              'respfile',params.respfile,...
              'respvarname',params.respvarname,...
              'resplen',params.resplen,...
              'repcount',params.repcount,...
              'spikes',params.spikes,...
              'respfilefmt',respfilefmt,...
              'respfmtcode',params.respfmtcode,...
              'nosync',params.nosync,...
              'stimfile',params.stimfile,...
              'stimpath',params.stimpath,...
              'stimwindowsize',params.stimwindowsize,...
              'stimfilefmt',stimfilefmt,...
              'stimfmtcode',stimfmtcode,...
              'stimfilecrf',params.stimfilecrf,...
              'stimspeedid',rawdata.stimspeedid,...
              'stimiconside',iconside,...
              'addedby',params.addedby,...
              'info',params.info);
fprintf('added scellfile entry %d\n',cellfileid);
