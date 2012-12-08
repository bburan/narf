% function [filedata,params]=lfpfiletimes(siteid,batchid/params);
%
% first parameter can be a cellid, in which case all sCellFile
% entries matching cellid and specs in params are loaded. if
% it's an array of numbers, all sCellFile entries with matching
% rawids (and in format specified by params) are loaded.
%
% second parameter is either batchid or data structure from a
% particular batch.  if batchid is supplied, params is loaded
% from sBatch with matching id.
%
% if you specify params explicitly, you need to provide [defaults]:
% params.resploadcmd    [='respload']
%          .resploadparms  [={'',1,1,1}]
%          .respfiltercmd  [='']
%          .respfilterparms [={}];
%          .predfrac [=0]  % frac of data for validation
%          .fitfrac [=0]   % fraction for intermediate fit
%          .runclassid [=2 (natrev)] can be vector of mult. values
%          .stimspeedid [=60] (hz of update rate) can be vector
%          .stimfmtcode [=0 (pixel domain)] 2=downsamp
%          .respfmtcode [=0 (psth)] 1=pfth
%
% created SVD 6/02
% modified SVD 1/03 - added batchid option
% 
function [cellfiledata,params]=lfpfiletimes(cellid,params,bquiet);

MAXPREDLEN=7000; % plenty data for fitting regularization/threshold parms?

dbopen;

if ~exist('bquiet','var'),
   bquiet=0;
end

if isnumeric(params),
   batchid=params;
   sql=['SELECT * FROM sBatch WHERE id=',num2str(batchid)];
   params=mysql(sql);
   if length(params.parmstring)>0,
      eval(char(params.parmstring));
   end
else
   %batchid=params.id;
end

params.stimspeedid=getparm(params,'stimspeedid',60);
params.runclassid=getparm(params,'runclassid',1);

if length(params.runclassid)>1,
   batchstr='runclassid in (';
   for ii=1:length(params.runclassid),
      batchstr=[batchstr,num2str(params.runclassid(ii)),','];
   end
   batchstr(end)=')';
else
   batchstr=['runclassid=',num2str(params.runclassid)];
end

if isnumeric(cellid),
   swhere=' WHERE rawid in (';
   for ii=1:length(cellid),
      swhere=[swhere,num2str(cellid(ii)),','];
   end
   swhere(end)=')';
else
   if params.stimspeedid==0,
      speedtest='';
   elseif length(params.stimspeedid)>1,
      speedtest=' AND stimspeedid in (';
      for ii=1:length(params.stimspeedid),
         speedtest=[speedtest,num2str(params.stimspeedid(ii)),','];
      end
      speedtest(end)=')';
      
   elseif params.stimspeedid>=3,
      speedtest=[' AND stimspeedid>=',num2str(params.stimspeedid)];
   else
      speedtest=[' AND stimspeedid=',num2str(params.stimspeedid)];
   end
   
   if isempty(params.stimsnr),
      snrtest='';
   elseif params.stimsnr>=100,
      snrtest=[' AND stimsnr>=',num2str(params.stimsnr)];
   else
      snrtest=[' AND stimsnr=',num2str(params.stimsnr)];
   end
   
   swhere=[' WHERE cellid="',cellid,'"',...
           ' AND ',batchstr, speedtest, snrtest];
   sql=['SELECT *,(repcount>2) as multrep FROM sCellFile',...
        swhere,...
        ' AND stimfmtcode=',num2str(stimfmtcode),...
        ' AND respfmtcode=',num2str(respfmtcode),...
        ' ORDER BY stimspeedid,repcount,resplen DESC,respfile'];
   cellfiledata=mysql(sql);
   if length(cellfiledata)==0,
      speedtest='';
      swhere=[' WHERE cellid="',cellid,'"',...
              ' AND ',batchstr, snrtest];
      sql=['SELECT *,(repcount>2) as multrep FROM sCellFile',...
           swhere,...
           ' AND stimfmtcode=',num2str(stimfmtcode),...
           ' AND respfmtcode=',num2str(respfmtcode),...
           ' ORDER BY stimspeedid,repcount,resplen DESC,respfile'];
   end
end

%
% new file selection scheme
%
if 0 & params.runclassid>0
   sql=['SELECT *,(repcount>1) as multrep FROM sCellFile',...
        swhere,...
        ' AND stimfmtcode=',num2str(stimfmtcode),...
        ' AND respfmtcode=',num2str(respfmtcode),...
        ' ORDER BY multrep,resplen,respfile'];
else
   sql=['SELECT *,(repcount>2) as multrep FROM sCellFile',...
        swhere,...
        ' AND stimfmtcode=',num2str(stimfmtcode),...
        ' AND respfmtcode=',num2str(respfmtcode),...
        ' ORDER BY stimspeedid,repcount,resplen DESC,respfile'];
end
cellfiledata=mysql(sql);

if length(cellfiledata)==0 | nargout<=1,
   times=[];
   return
end

params.stimloadcmd=getparm(params,'stimloadcmd','loadimfile');
params.stimloadparms=getparm(params,'stimloadparms','');
%params.stimloadparms=strsep(params.stimloadparms,',');
params.stimfiltercmd=getparm(params,'stimfiltercmd','');
params.stimfilterparms=getparm(params,'stimfilterparms','');
%params.stimfilterparms=strsep(params.stimfilterparms,',');
params.resploadcmd=getparm(params,'resploadcmd','respload');
params.resploadparms=getparm(params,'resploadparms',',1,1,1');

params.resploadparms=strsep(params.resploadparms,',');
params.respfiltercmd=getparm(params,'respfiltercmd','');
params.respfilterparms=getparm(params,'respfilterparms','');
params.respfilterparms=strsep(params.respfilterparms,',');

% find out info about response files:
params.cellid=cellid;
params.resplen=zeros(1,length(cellfiledata));
params.respfiles={};
params.stimcrfs=zeros(1,length(cellfiledata));

params.repcount=cat(2,cellfiledata.repcount);
params.resplen=cat(2,cellfiledata.resplen);

for ii=1:length(cellfiledata),
   params.respfiles{ii}=[cellfiledata(ii).path,cellfiledata(ii).respfile];
   params.stimfiles{ii}=[cellfiledata(ii).stimpath,...
                    cellfiledata(ii).stimfile];
   params.stimcrfs(ii)=ifstr2num(cellfiledata(ii).stimfilecrf);
   if cellfiledata(ii).stimfilecrf==0,
      tstimpix=strsep(cellfiledata(ii).stimiconside,',');
      if length(tstimpix)>0,
         tstimpix=tstimpix{1};
      end
      sql=['SELECT * FROM gSingleCell WHERE id=',...
           num2str(cellfiledata(ii).singleid)];
      celldata=mysql(sql);
      if celldata.rfsize>0,
         params.stimcrfs(ii)=tstimpix./celldata.rfsize;
      end
   end
end

[s,host]=unix('hostname');
if onseil,
   disp('mapping stim/resp files to seil');
   for ii=1:length(cellfiledata),
      trin=params.respfiles{ii};
      tsin=params.stimfiles{ii};
      trout=strrep(trin,'/auto/data/','/homes/svd/data/');
      tsout=strrep(tsin,'/auto/data/','/homes/svd/data/');
      [bb,pp]=basename(trout);
      if strcmp(bb(1:5),'model'),
         disp('model cell, forcing recopy of response');
         ['\rm ' trout]
         unix(['\rm ' trout]);
      end
      if ~exist(trout,'file'),
         unix(['mkdir -p ',pp]);
         ['scp svd@bhangra.isr.umd.edu:',trin,' ',pp]
         unix(['scp svd@bhangra.isr.umd.edu:',trin,' ',pp]);
      end
      if ~exist(tsout,'file') & ~exist([tsout '.m'],'file')
         [bb,pp]=basename(tsout);
         unix(['mkdir -p ',pp]);
         ['scp svd@bhangra.isr.umd.edu:',tsin,' ',pp]
         unix(['scp svd@bhangra.isr.umd.edu:',tsin,' ',pp]);
         ['scp svd@bhangra.isr.umd.edu:',tsin,'.m ',pp]
         unix(['scp svd@bhangra.isr.umd.edu:',tsin,'.m ',pp]);
      end
      if ~isempty(findstr(tsin,'.par')),
         tsin2=strrep(tsin,'.par','.inf');
         tsout2=strrep(tsout,'.par','.inf');
         if ~exist(tsout2,'file'),
            [bb,pp]=basename(tsout2);
            ['scp svd@bhangra.isr.umd.edu:',tsin2,' ',pp]
            unix(['scp svd@bhangra.isr.umd.edu:',tsin2,' ',pp]);
         end
      end
      params.respfiles{ii}=trout;
      params.stimfiles{ii}=tsout;
   end
end

[times,params]=xcfilefracs(params);

