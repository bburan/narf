% script to insert a new batch in sBatch and associated sRunData entries.

clear batch

MINLEN=5000; % minimum total length of respfiles allowed into
             % batch. 0 lets everyone in.

rcsetstrings;

batch.runclassid=3;
%|  0 | NAT   | natural vision movie                  |
%|  1 | DGRAT | dynamic grating sequence              |
%|  2 | DNAT  | dynamic natural image sequence        |
%|  3 | FS    | free view search                      |
%|  4 | AFS   | attention free view search (obsolete) |
%|  5 | WG    | wavy gravy                            |
%|  6 | SR    | dynamic sparse noise                  |

batch.stimspeedid=0;   % 0 = default/fast, 1=slow (8 Hz NR and WG)

batch.stimfmtcode=2;
%| stimfmtcode | stimfilefmt |
%|           0 | pixel       |
%|           1 | logfft2     |
%|           2 | downsamp    |
%|           4 | mfilt       |
%|           5 | wav         |
%|           6 | parm        |
%|           7 | sfft        |
batch.respfmtcode=0;  %  0=PSTH, 1=PFTH
batch.attstate=1;     % 0 = none, 1=feature

batch.kernfmt='pfft';
%batch.kernfmt='mfilt';
%batch.kernfmt='wav';
%batch.kernfmt='gr';
%batch.kernfmt='fftgr';
%batch.kernfmt='wg';

% stim filter - eg, convert to phase-sep fourier domain.
%batch.stimfiltercmd='';
%batch.stimfilterparms='';
%batch.stimfiltercmd='mov2gr';
%batch.stimfilterparms='8,8,0,0,0,1,1';
%batch.stimfiltercmd='movcrop';
%batch.stimfilterparms='[1 1 4]';
%batch.stimfiltercmd='mov2wav';
%batch.stimfilterparms='4,4,3';
batch.stimfiltercmd='movpower';
batch.stimfilterparms='0,0,0,1,1'; % [startidx,stopidx,bnorm,bpower,bzdc]
%batch.stimfiltercmd='movphasesep';
%batch.stimfilterparms='0,0,0,1,1'; % [startidx,stopidx,bnorm,bpower,bzdc]

batch.decorrspace=2;
batch.decorrtime=1;
batch.sfscount=40;
batch.sfsstep=4;


if batch.stimspeedid==0,
   sspeed='';
else
   sspeed=speedstr{batch.stimspeedid+1};
end
if batch.attstate==0,
   satt='';
else
   satt=['/',attstr{batch.attstate+1}];
end

batch.name=[runclassstr{batch.runclassid+1},sspeed,'/',...
            respfilefmtstr{batch.respfmtcode+1},'/',...
            batch.kernfmt,satt];

batch.details=input(['Details for ',batch.name,': '],'s');


% resp load cmd. format: resploadcmd(respfile,p1,p2,...)
if batch.attstate==0,
   batch.resploadcmd='loadpypedata';
   batch.resploadparms=[',1,0'];  % path,
   
   % resp filter - eg, rebin 
   batch.respfiltercmd='';
   batch.respfilterparms='';
else
   batch.resploadcmd='resploadatt';
   if batch.respfmtcode==0,
      batch.resploadparms=['psth,',num2str(batch.attstate),',0,1']; ...
   else
      batch.resploadparms=['pfth,',num2str(batch.attstate),',0,1'];
   end
   batch.respfiltercmd='respvarbins';
   batch.respfilterparms=['[1 51 101 101 101 101 151 151 151 200 200 250],',...
                          '[100 100 150 200 250 300 200 250 300 250 300 300]'];
end

if batch.stimfmtcode==6,
   % stim load cmd. should be of format
   batch.stimloadcmd='loadwgstim';
   batch.stimloadparms='4,4,5';  % obincount,sbincount,spacebincount
else
   batch.stimloadcmd='loadimfile';
   
   if strcmp(batch.kernfmt,'gr'),
      batch.stimloadparms=''; 
   elseif ismember(batch.stimfmtcode,[2]),
      batch.stimloadparms='10,17,10'; % downsample for better snr
   else
      batch.stimloadparms='';
   end
end

if batch.respfmtcode==0,
   batch.minlag=-8;
   batch.maxlag=15;
else
   batch.minlag=0;
   batch.maxlag=0;
end
batch.resampcount=20;
batch.resampfmt=1;   % ie bootstrapping
batch.fitfrac=0.1;
batch.predfrac=0.1;
batch.hfiltsigma=2;
batch.sffiltsigma=7;
batch.srfiltsigma=0;  % obsolete?
batch.sffiltthresh=0; % obsolete?
batch.sffiltsmooth=0; % obsolete?
batch.predsmoothsigma=1;
batch.predtype=0;
%batch.sfscount=50;
%batch.sfsstep=1;

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

%batchid=1;
sql=['INSERT INTO sBatch ',nstr,vstr];
[result,affected,batchid]=mysql(sql);

rcsetstrings;
resbase=['.',runclassstr{batch.runclassid+1},sspeed,...
         '.',batch.kernfmt, ...
         '.',respfilefmtstr{batch.respfmtcode+1},'.res.mat'];
kernbase=['.',runclassstr{batch.runclassid+1},sspeed,...
          '.',batch.kernfmt, ...
          '.',respfilefmtstr{batch.respfmtcode+1},'.kern.mat'];
sql=['SELECT DISTINCT masterid,cellid,path,sum(resplen) as totlen',...
     ' FROM sCellFile',...
     ' WHERE runclassid=',num2str(batch.runclassid),...
     ' AND stimspeedid=',num2str(batch.stimspeedid),...
     ' AND respfmtcode=',num2str(batch.respfmtcode),...
     ' AND stimfmtcode=',num2str(batch.stimfmtcode),...
     ' AND cellid>"m006"',...
     ' AND cellid<>"R240A"',...
     ' GROUP BY cellid',...
     ' ORDER BY cellid,path'];
filedata=mysql(sql);

fprintf('Created batch id %d, runcount=%d, resbase=%s\n',...
        batchid,length(filedata),resbase);

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


