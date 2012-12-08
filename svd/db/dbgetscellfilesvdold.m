% function [cellfiledata,cellids,cellfileids]=dbgetscellfile([field,filter],...)
%
% return the cellfiledata and cellfileids masked by filter specified
% by parameters. parameter format is ('field',value).
%
% EXAMPLES:
% dbgetscellfile('cellid','R110A') - return all data for cell R110A.
% dbgetscellfile('cellid','R110A','runclassid',0) - return all
%                                    review data for cell R110A.
%
% PARAMETERS:
% 'rawid' - id in gDataRaw table (select a single data file)
% 'cellid' - name of cell, can be partial. eg, 'R' will return all
%            cells with the name "R*" (R110A, R210A, R225C, etc)
% 'speed' - speed of stimulus (in Hz)
% 'runclassid' - code of task
%             -1  all (default)
%              0 | NAT   | natural vision movie (Review)           |
%              1 | DGRAT | dynamic grating sequence                |
%              2 | DNAT  | dynamic natural image sequence (NatRev) |
%              3 | FS    | free view search                        |
%              4 | AFS   | attention free view search (obsolete)   |
%              5 | WG    | wavy gravy                              |
%              6 | SR    | dynamic sparse noise                    |
%              7 | MPEG  | mpeg freeview                           |
%              8 | EYEC  | eye calibration                         |
%              9 | IMV   | imview                                  |
%           | 10 | DMS   | delayed match-to-sample         |
%             11 | --    | non-task data file              |
%           | 12 | OPTI  | opti nat                        |
%           | 13 | MSEQ  | m-sequence                      |
% 'fmt' - file of preprocessed stimulus:
%            'pixel' - raw pixel domain (all reese data <= R225C)
%            'parm' - parametric stim (gratrev and wavygrat)
%            'sfft' for phase separated fourier domain (new NR)
%            'wav' for wavelet domain
%            'pfft' for fourier power
%            'downsamp' for 16x16 pixel downsampled
% 'respfmtcode' - response format
%                 0: PSTH (binned in time 14 or 16 ms)
%                 1: PFTH (fixation aligned response)
% 'area' - match area field in gSingleCell (% is the wildcard) so
%          'area','V1%' will return data for all cells whose area
%          field begins with 'V1'
%
% RETURNS (example):
% cellfiledata( ).id: 795                         % id in sCellFile
%             cellid: 'R110A'
%           masterid: 82                          % id of cell in gCellMaster
%              rawid: 209                         % id in gDataRaw
%         celldataid: 183
%         runclassid: 1                              % ie, gratrev
%               path: '/auto/k5/david/data/R110A/'   % path to respfile
%            resplen: 21467
%           repcount: 1
%           respfile: 'R110A.gratsize.1.all.d-svsizemult1.0-.psth.neg_flag'
%        respvarname: 'R110A'
%       respfiletype: 1
%             nosync: 0
%        respfilefmt: '14ms PSTH'
%        respfmtcode: 0
%           stimfile: 'dat.R110A.gratrev.imsm'
%       stimfiletype: 1
%       stimiconside: '32,32'
%        stimfilecrf: 2
%     stimwindowsize: 32
%        stimfilefmt: 'pixel'                    % 'fmt' name
%        stimfmtcode: 0
%            addedby: 'david'
%               info: 'dbcrap.m - from tCellFile'
%            lastmod: '20030220152821'
%           stimpath: '/auto/k5/david/data/R110A/'   % path to stimfile
%        stimspeedid: 0
%             spikes: 13270
%            multrep: 0
% cellids: cell array of unique cellids matched in search
%
% CREATED SVD 2/28/03
%
function [cellfiledata,cellids,cellfileids]=dbgetscellfile(varargin)

dbopen;

narg = 1;
wherestr='WHERE 1';
while narg <= length(varargin)
   switch varargin{narg}
    case 'rawid'
     rawids = varargin{narg + 1};
     srawid='(';
     for ii=1:length(rawids),
        srawid=[srawid,num2str(rawids(ii)),','];
     end
     srawid(end)=')';
     wherestr=[wherestr,' AND sCellFile.rawid in ',srawid];
     
    case 'cellid'
     wherestr=[wherestr,' AND sCellFile.cellid like "',varargin{narg + 1},'%"'];

    case 'fmt'
     wherestr=[wherestr,' AND sCellFile.stimfilefmt="',varargin{narg + 1},'"'];

    case 'runclassid'
     wherestr=[wherestr,' AND sCellFile.runclassid=',num2str(varargin{narg + 1})];
    
    case 'speed'
     wherestr=[wherestr,' AND sCellFile.stimspeedid=',num2str(varargin{narg + 1})];
    
    case 'respfmtcode'
     wherestr=[wherestr,' AND sCellFile.respfmtcode=',num2str(varargin{narg + 1})];
    
    case 'area'
     wherestr=[wherestr,' AND gSingleCell.area like "',(varargin{narg + 1}),'"'];
    
    otherwise
     error(sprintf('unknown option: %s', varargin{narg}));
     narg = narg - 1;
   end
   narg = narg + 2;
end

sql=['SELECT sCellFile.*,area,gDataRaw.bad',...
     ' FROM (sCellFile INNER JOIN gDataRaw ON sCellFile.rawid=gDataRaw.id)',...
     ' INNER JOIN gSingleCell' ...
     ' ON sCellFile.singleid=gSingleCell.id ',...
     wherestr,...
     ' AND not(gDataRaw.bad) AND not(gSingleCell.crap)',...
     ' ORDER BY cellid'];

cellfiledata=mysql(sql);
cellfileids=cat(1,cellfiledata.id);

sql=['SELECT DISTINCT sCellFile.cellid',...
     ' FROM sCellFile INNER JOIN gSingleCell',...
     ' ON sCellFile.singleid=gSingleCell.id ',...
     wherestr,' ORDER BY cellid'];
celliddata=mysql(sql);
cellids={celliddata.cellid};

