% function [cellfileids,cellfiledata]=dbraw2scellfile(rawids,fmt)
%
% (if not already done) convert to pype files to matlab-readable
% format and translate entries from gDataRaw to sCellFile
%
% if the conversion has already been done, simply return the
% cellfileids corresponding to the rawids and fmt specified.
%
% for WG: use fmt='parm'
% for NR: use fmt='sfft' for phase separated fourier domain
%         use fmt='wav' for wavelet domain
%
function [cellfileids,cellfiledata]=dbraw2scellfile(rawids,fmt)

dbopen;

% stimulus imsm archive path may vary with the well.
% after the end of a session, path should be updated to show archived locations of data.
stimpathdata.reese{13}='/auto/data/archive/stimarchive/REESE13/imsm/';
stimpathdata.reese{14}='/auto/data/archive/stimarchive/REESE14/imsm/';
stimpathdata.reese{15}='/auto/data/archive/stimarchive/REESE15/imsm/';
stimpathdata.reese{16}='/auto/data/archive/stimarchive/REESE16/imsm/';
stimpathdata.reese{17}='/auto/data/stimarchive/REESE17/imsm/';
stimpathdata.mac{11}='/auto/data/archive/stimarchive/MAC11/imsm/';
stimpathdata.van{1}='/auto/c1/archive/stimarchive/VAN01/imsm/';

if exist('fmt','var'),
   fmtstr=[' AND stimfilefmt="',fmt,'"'];
else
   fmt='';
   fmtstr='';
end

if ~exist('rawids'),
   RUNCLASSID=2;
   sql=['SELECT gDataRaw.* FROM gDataRaw',...
        ' WHERE isnull(sCellFile.rawid)',...
        ' AND gDataRaw.runclassid=',num2str(RUNCLASSID),...
        ' AND bad=0',...
        ' ORDER BY gDataRaw.cellid,gDataRaw.id'];
   rawdata=mysql(sql);
   
   rawids=cat(1,rawdata.id);
   GENERALUPDATE=1;
   keyboard
else
   srawid='(';
   for ii=1:length(rawids),
      srawid=[srawid,num2str(rawids(ii)),','];
   end
   srawid(end)=')';
   sql=['SELECT gDataRaw.*',...
        ' FROM gDataRaw',...
        ' WHERE gDataRaw.id in ',srawid,...
        ' AND bad=0',...
        ' ORDER BY gDataRaw.cellid,gDataRaw.id'];
   rawdata=mysql(sql);
   GENERALUPDATE=0;
end
cellfileids=[];

for rawidx=1:length(rawdata),
   
   sql=['SELECT * FROM gCellMaster where id=',...
        num2str(rawdata(rawidx).masterid)];
   celldata=mysql(sql);
   
   sql=['SELECT * FROM sCellFile',...
        ' WHERE rawid=',num2str(rawdata(rawidx).id),...
        fmtstr];
   cellfiledata=mysql(sql);
   
   if ~GENERALUPDATE & length(cellfiledata)>0,
      cellfileids=[cellfileids;cat(1,cellfiledata.id)];
   else
      if strcmp(rawdata(rawidx).cellid,'model'),
         matfile=rawdata(rawidx).matlabfile;
      elseif strcmp(rawdata(rawidx).cellid(1),'R') & ...
            ~isempty(rawdata(rawidx).matlabfile),
         if ~isempty(findstr(rawdata(rawidx).matlabfile,'0.5')) | ...
               ~isempty(findstr(rawdata(rawidx).matlabfile,'0.3')) | ...
               ~isempty(findstr(rawdata(rawidx).matlabfile,'0.2')),
            matfile=[];
            disp('skipping');
         else
            matfile=[rawdata(rawidx).resppath,rawdata(rawidx).matlabfile];
         end
      else
         keyboard
         matfile=processpypedata([rawdata(rawidx).resppath,...
                    rawdata(rawidx).respfile]);
      end
      
      if ~findstr(rawdata(rawidx).info,'CELL'),
         
         %ie, this is an old (SGI-collected) reese cell
         
         z=load(matfile);
         resplen=sum((z(:,1)>-1));
         if strcmp(rawdata(rawidx).cellid,'model')
            respvarname='psth';
         else
            respvarname=rawdata(rawidx).cellid;
         end
         respfmtcode=0;
         respfilefmt='PSTH'; % rather than PFTH
         respfile=rawdata(rawidx).matlabfile;
         resppath=rawdata(rawidx).resppath;
         
         stimpath=resppath;
         stimfile={rawdata(rawidx).stimfile};
         stimfmtcode={0};
         stimfilefmt={'pixel'};
         [rlen,stimwindowsize]=imfileinfo([rawdata(rawidx).resppath,...
                    rawdata(rawidx).stimfile]);
         stimwindowsize=stimwindowsize(1);
         stimspeedid=0;
      else
         % for CELLDB cells
         
         z=load(matfile);
         resplen=sum(~isnan(z.psth(:,1)));
         respvarname='psth';
         respfmtcode=0;
         respfilefmt='PSTH'; % rather than PFTH
         respfile=[rawdata(rawidx).respfile,'.mat'];
         resppath=rawdata(rawidx).resppath;
         stimspeedid=rawdata(rawidx).stimspeedid;
         
         try,
           stimpath=z.h.play_dir0;
           indexfile=z.h.index0;
         catch,
           stimpath=z.h.moviedir;
           indexfile=z.h.indexfile;
         end
         
         % kludgey hack thing SVD 9/30/02
         % for situations where header of index file doesn't
         % have these parameters
         if ~isfield(z.h,'stim'),
            % guess a value for the stim type from runclassid
            if rawdata(rawidx).runclassid==1,
               z.h.stim='gr';
            else
               z.h.stim='nr';
            end
         end
         if ~isfield(z.h,'pix'),
            z.h.pix=0;
         end
         
         if strcmp(z.h.stim(1:2),'wg'),
            % wavy gravy
            stimfile={respfile};
            stimfmtcode={6};
            stimfilefmt={'parm'}; % paramtric RC
            stimpath=resppath;
            stimwindowsize=z.h.pix;
         elseif strcmp(z.h.stim(1:2),'gr') | strcmp(z.h.stim(1:2),'nc'),
            % gratrev
            stimfile={respfile};
            stimfmtcode={6};
            stimfilefmt={'parm'}; % paramtric RC
            stimpath=resppath;
            stimwindowsize=z.h.pix;
         elseif strcmp(z.h.stim(1:2),'be') | strcmp(z.h.stim(1:2),'sy'),
            
            if stimbase(end)=='/',
               stimbase=stimbase(1:(end-1));
            end
            ii=max(findstr(stimbase,'/'));
            if ii-length(stimbase)<5,
               stimbase(ii)='-';
               ii=max(findstr(stimbase,'/'));
            end
            stimfile={[stimbase((ii+1):end),'.',indexfile,'.pix.imsm']};
            stimfmtcode={0};
            stimfilefmt={'pixel'};
            sql=['SELECT * FROM gCellMaster where id=',...
                 num2str(rawdata(rawidx).masterid)];
            celldata=mysql(sql);
            stimpathidx=celldata.well;
            stimpathlist=getfield(stimpathdata,celldata.animal);
            
            if stimpathidx<=length(stimpathlist),
               stimpath=stimpathlist{stimpathidx};
            else
               fprintf('stimpath not found for %s!\n',stimfile{1});
               input('enter path: ',stimfile);
            end
            stimwindowsize=z.h.pix;
         else
            % natrev/review data.
            
            % figure out construction of pre-computed imsm file:
            try,
              stimbase=z.h.play_dir0;
            catch,
              stimbase=z.h.moviedir;
            end
            if stimbase(end)=='/',
               stimbase=stimbase(1:(end-1));
            end
            ii=max(findstr(stimbase,'/'));
            
            % for new natrev naming scheme where dir is only size
            if ~isnan(str2double(stimbase((ii+1):end)))
               stimbase(ii)='-';
               ii=max(findstr(stimbase,'/'));
            end
            
            stimbase=[stimbase((ii+1):end),'.',indexfile,'.'];
            
            % 3 formats for each NR file: downsamp, fft, wav
            stimfmtcode={2,5,7,8};
            stimfilefmt={'downsamp','wav','sfft','pfft'};
            stimfile={[stimbase,'pix.imsm'],[stimbase,'wav.imsm'],...
                      [stimbase,'fft.imsm'],[stimbase,'pfft.imsm']};
            
            sql=['SELECT * FROM gCellMaster where id=',...
                 num2str(rawdata(rawidx).masterid)];
            celldata=mysql(sql);
            stimpathidx=celldata.well;
            stimpathlist=getfield(stimpathdata,celldata.animal);
            
            if stimpathidx<=length(stimpathlist),
               stimpath=stimpathlist{stimpathidx};
            else
               fprintf('stimpath not found for %s!\n',stimfile{1});
               stimfile{1}=input('enter path: ');
            end
            
            % this following thing SHOULD work, but above hack is passable
            stimpix=z.h.pix;
            stimwindowsize=round(stimpix/2);
            
         end
      end
      
      if resppath(end)~='/',
         resppath=[resppath,'/'];
      end
      
      for stimfmtidx=1:length(stimfmtcode),
         [aff,rawdata(rawidx).cellfileid]=...
             sqlinsert('sCellFile',...
                       'cellid',rawdata(rawidx).cellid,...
                       'masterid',rawdata(rawidx).masterid,...
                       'rawid',rawdata(rawidx).id,...
                       'runclassid',rawdata(rawidx).runclassid,...
                       'path',resppath,...
                       'resplen',resplen,...
                       'respfile',respfile,...
                       'respvarname',respvarname,...
                       'respfiletype',1,...
                       'respfilefmt',respfilefmt,...
                       'respfmtcode',respfmtcode,...
                       'stimfile',stimfile{stimfmtidx},...
                       'stimpath',stimpath,...
                       'stimwindowsize',stimwindowsize,...
                       'stimfilefmt',stimfilefmt{stimfmtidx},...
                       'stimfmtcode',stimfmtcode{stimfmtidx},...
                       'stimspeedid',stimspeedid,...
                       'addedby','david',...
                       'info','dbraw2scellfile.m');
         if GENERALUPDATE | strcmp(stimfilefmt{stimfmtidx},fmt),
            cellfileids=[cellfileids;rawdata(rawidx).cellfileid];
         end
      end
      
      % update rawdata stuff if not already there
      if isempty(rawdata(rawidx).stimfile),
         try,
           stimpath=z.h.play_dir0;
           stimfile=z.h.index0;
         catch,
           stimpath=z.h.moviedir;
           stimfile=z.h.indexfile;
         end
         if ~strcmp(stimpath(end),'/'),
            stimpath=[stimpath,'/'];
         end
         sql=['UPDATE gDataRaw SET',...
              ' stimpath="',stimpath,'",',...
              ' stimfile="',stimfile,'"',...
              ' WHERE id=',num2str(rawdata(rawidx).id)];
         mysql(sql);
      end
   end
end

[cellfiledata,cellids,cellfileids]=dbgetscellfile('rawid',rawids,'fmt',fmt);
