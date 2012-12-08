% script to insert a new batch in sBatch and associated sRunData entries.

clear batch

MINLEN=0; % minimum total length of respfiles allowed into
             % batchdata. 0 lets everyone in.

rcsetstrings;

% RUNCLASSID: stimulus class
%|  0 | NAT   | natural vision movie                  |
%|  1 | DGRAT | dynamic grating sequence              |
%|  2 | DNAT  | dynamic natural image sequence        |
%|  3 | FS    | free view search                      |
%|  4 | AFS   | attention free view search (obsolete) |
%|  5 | WG    | wavy gravy                            |
%|  6 | SR    | dynamic sparse noise                  |
%|  7 | MPEG  | mpeg freeview                           |
%|  8 | EYEC  | eye calibration                         |
%|  9 | IMV   | imview                                  |

ADDRUNDATA=1;
ONLINE=0;

batchdata.runclassid=7;
batchdata.stimspeedid=0;   % 0 = default/fast, 1=slow (8 Hz NR and WG)

%| stimfmtcode | stimfilefmt |
%|           0 | pixel       |
%|           1 | logfft2     |
%|           2 | downsamp    |
%|           4 | mfilt       |
%|           5 | wav         |
%|           6 | parm        |
%|           7 | sfft        |
batchdata.stimfmtcode=2;
batchdata.respfmtcode=0;  %  0=PSTH, 1=PFTH
batchdata.attstate=0;     % 0 = none, 1=feature

batchdata.stimwindowcrf=2;

%batchdata.kernfmt='space';
batchdata.kernfmt='pfft';
%batchdata.kernfmt='fft';
%batchdata.kernfmt='wav';
%batchdata.kernfmt='lin2';
%batchdata.kernfmt='mfilt';
%batchdata.kernfmt='wg';

% stim filter - eg, convert to phase-sep fourier domain.
batchdata.stimfiltercmd='';
batchdata.stimfilterparms='';
if ismember(batchdata.stimfmtcode,[0 2]) & strcmp(batchdata.kernfmt,'pfft'),
   batchdata.stimfiltercmd='movpower';
   batchdata.stimfilterparms='0,0,0,1,0'; %[startidx,stopidx,bnorm,bpower,bzdc]
elseif ismember(batchdata.stimfmtcode,[0 2]) & strcmp(batchdata.kernfmt,'fft'),
   batchdata.stimfiltercmd='movphasesep';
   batchdata.stimfilterparms='0,0,0,1,0'; % [startidx,stopidx,bnorm,bpower,bzdc]
elseif ismember(batchdata.stimfmtcode,[0 2]) & strcmp(batchdata.kernfmt,'wav'),
   batchdata.stimfiltercmd='mov2wav';
   batchdata.stimfilterparms='4,4,3';
elseif ismember(batchdata.stimfmtcode,[0 2]) & strcmp(batchdata.kernfmt,'lin2'),
   batchdata.stimfiltercmd='movlinpower';
   batchdata.stimfilterparms='0,0,0,1,0';
end

if ONLINE,
   batchdata.nloutparm=1;
   batchdata.decorrspace=3;
   batchdata.sfscount=30;
   batchdata.fitfrac=0.15;
   batchdata.predfrac=0.0;
else
   batchdata.nloutparm=3;
   batchdata.decorrspace=2;
   batchdata.sfscount=50;
   batchdata.fitfrac=0.0;
   batchdata.predfrac=0.1;
end

batchdata.sfsstep=1;
batchdata.decorrtime=1;
batchdata.resampfmt=1;   % ie jackknifing
batchdata.resampcount=20;
batchdata.sffiltsigma=7;

if batchdata.stimspeedid==0,
   sspeed='';
else
   sspeed=['.',num2str(batchdata.stimspeedid)];
end
if batchdata.attstate==0,
   satt='';
else
   satt=['/',attstr{batchdata.attstate+1}];
end

batchdata.name=[runclassstr{batchdata.runclassid+1},sspeed,'/',...
            respfilefmtstr{batchdata.respfmtcode+1},'/',...
            batchdata.kernfmt,satt];

batchdata.details=input(['Details for ',batchdata.name,': '],'s');


% resp load cmd. format: resploadcmd(respfile,p1,p2,...)
if batchdata.attstate==0 & batchdata.runclassid==7,
   batchdata.resploadcmd='loadimfile';
   batchdata.resploadparms='';
   
   % resp filter - eg, rebin 
   batchdata.respfiltercmd='';
   batchdata.respfilterparms='';
elseif batchdata.attstate==0 & batchdata.stimfmtcode==0,
   batchdata.resploadcmd='respload';
   batchdata.resploadparms=[',1,1,1'];
   
   % resp filter - eg, rebin 
   batchdata.respfiltercmd='';
   batchdata.respfilterparms='';
elseif batchdata.attstate==0,
   batchdata.resploadcmd='loadpypedata';
   batchdata.resploadparms=[',1,0'];  % path,
   
   % resp filter - eg, rebin 
   batchdata.respfiltercmd='';
   batchdata.respfilterparms='';
else
   batchdata.resploadcmd='resploadatt';
   if batchdata.respfmtcode==0,
      batchdata.resploadparms=['psth,',num2str(batchdata.attstate),',0,1']; ...
   else
      batchdata.resploadparms=['pfth,',num2str(batchdata.attstate),',0,1'];
   end
   batchdata.respfiltercmd='respvarbins';
   batchdata.respfilterparms=['[1 51 101 101 101 101 151 151 151 200 200 250],',...
                          '[100 100 150 200 250 300 200 250 300 250 300 300]'];
end

if batchdata.stimfmtcode==6,
   % stim load cmd. should be of format
   batchdata.stimloadcmd='loadwgstim';
   batchdata.stimloadparms='4,4,5';  % obincount,sbincount,spacebincount
else
   batchdata.stimloadcmd='loadimfile';
   
   if ismember(batchdata.runclassid,[7]),
      batchdata.stimloadparms='16,8,16';   % mpeg, already downsampled
   elseif ismember(batchdata.stimfmtcode,[2]),
      batchdata.stimloadparms='16,17,16';   % already downsampled
   elseif ismember(batchdata.stimfmtcode,[0]),
      % v1 analysis, figure out scale parameter during analysis
      batchdata.stimloadparms='0,8,16';
   else
      batchdata.stimloadparms='';
   end
end

if batchdata.respfmtcode==0,
   batchdata.minlag=-6;    % these lags seem to work for v1
   batchdata.maxlag=13;
else
   batchdata.minlag=0;
   batchdata.maxlag=0;
end
batchdata.hfiltsigma=2;
batchdata.srfiltsigma=0;  % obsolete?
batchdata.sffiltthresh=0; % obsolete?
batchdata.sffiltsmooth=0; % obsolete?
batchdata.predsmoothsigma=1;
batchdata.predtype=0;
%batchdata.sfscount=50;
%batchdata.sfsstep=1;

nstr='(';
vstr=' VALUES (';
flist=fields(batch);
for ii=1:length(flist),
   nstr=[nstr,flist{ii},','];
   t=getfield(batch,flist{ii});
   if isnumeric(t),
      vstr=[vstr,num2str(t),','];
   else
      vstr=[vstr,'"',t,'",'];
   end
end
nstr(end)=')';
vstr(end)=')';

rcsetstrings;
resbase=['.',runclassstr{batchdata.runclassid+1},sspeed,...
         '.',batchdata.kernfmt, ...
         '.',respfilefmtstr{batchdata.respfmtcode+1},'.res.mat'];
kernbase=['.',runclassstr{batchdata.runclassid+1},sspeed,...
          '.',batchdata.kernfmt, ...
          '.',respfilefmtstr{batchdata.respfmtcode+1},'.kern.mat'];
sql=['SELECT DISTINCT masterid,cellid,path,sum(resplen) as totlen',...
     ' FROM sCellFile',...
     ' WHERE runclassid=',num2str(batchdata.runclassid),...
     ' AND stimspeedid=',num2str(batchdata.stimspeedid),...
     ' AND respfmtcode=',num2str(batchdata.respfmtcode),...
     ' AND stimfmtcode=',num2str(batchdata.stimfmtcode),...
     ' AND cellid<>"R240A"',...
     ' AND cellid<>"model"',...
     ' GROUP BY cellid',...
     ' ORDER BY cellid,path'];
filedata=mysql(sql);

disp('paused before db action.');
keyboard

%batchid=1;
sql=['INSERT INTO sBatch ',nstr,vstr];
[result,affected,batchid]=mysql(sql);

if ADDRUNDATA,
   
   for ii=1:length(filedata),
      if filedata(ii).totlen>MINLEN,
         sql=['INSERT INTO sRunData ',...
              ' (masterid,cellid,batch,respath,resfile,kernfile)',...
              ' VALUES (',num2str(filedata(ii).masterid),',',...
              '"',filedata(ii).cellid,'",',...
              num2str(batchid),',',...
              '"',filedata(ii).path,'",',...
              '"',[filedata(ii).cellid,resbase],'",',...
              '"',[filedata(ii).cellid,kernbase],'")'];
         mysql(sql);
      end
   end
   
   fprintf('Created batch id %d, runcount=%d, resbase=%s\n',...
           batchid,length(filedata),resbase);

else
   fprintf('Created batch id %d\n',batchid);
end