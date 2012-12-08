% function [hf,stephf,movformat,kern]=rcgronline(rawids/filename,fftflag);
%
% rawids - vector of id tags from gDataRaw in cellDB. raw id is
%          listed next to the file name for each data set
% fftflag - if 1 (default) convert gr parms into pfft domain,
%                (scaled to 8 cyc/2 crf standard)
%           if 0 do rc in gr parameter space (or,sf)
%
% returns:
% hf - STRF estimate in space X time dimensions
% stephf - step response of STRF (ie, cumsum(hf,2))
% movformat - parameters used for RC that can be used to guide
%             optinat stimulus generation
% kern - spatial only kernel. slice of stephf at maximum latency
%        (ie, greatest power before inhibition/adaptation kicks in)
%
% created SVD 2/5/03 -- ripped off of rconline.m
%
function [hf,stephf,movformat,kern]=rcgronline(rawids,fftflag);

% rawfileids gives entries into gDataRaw for data sets to analyze
if ~exist('rawids','var'),
   rawids=[627];
end
if ~exist('fftflag','var'),
   fftflag=1;
end

disp('rcgronline.m:');

% cellfileids give entries into sCellFile for data sets to analyze.
% if this is a new cell, dbraw2scellfile will load relevant info
% from gDataRaw and create appropriate entries in sCellFile in
% order to return valid indices
dbopen;
BATQUEUEID=[];

fmt='parm';

if isnumeric(rawids),
   pypefiles={};
   for ii=1:length(rawids),
      sql=['SELECT * FROM gDataRaw WHERE id=',num2str(rawids(ii))];
      rawdata=mysql(sql);
      if length(rawdata)>0,
         pypefiles{ii}=[rawdata(1).resppath,rawdata(1).respfile];
      end
   end
elseif iscell(rawids),
   pypefiles=rawids;
   cellid=strsep(basename(pypefiles{1}));
   cellid=cellid{1};
else
   pypefiles{1}=rawids;
   cellid=strsep(basename(pypefiles{1}));
   cellid=cellid{1};
end

filecount=length(pypefiles);
matfiles={};
for ii=1:filecount,
   matfiles{ii}=processpypedata(pypefiles{ii});
end   

fitfrac=0;
predfrac=0;

disp('Loading data...');
ii=1
r=load(matfiles{ii});

% stim is in columns: [sf or ph];
resp=r.psth(:,1);
resp=resp-nanmean(resp);
resp(find(isnan(resp)))=0;

GRPHASE=0;

if strcmp(r.h.stim(1:2),'gr') & ~fftflag,
   disp('Gratrev: binning in parameter space...');
   
   a=r.s(r.sid+1,3);
   f=r.s(r.sid+1,4);
   ph=r.s(r.sid+1,5);
   
   % assume f goes 0-ceil(max(f)), a goes 0 to pi
   fmax=ceil(max(f));
   fbins=8;
   obins=8;
   fdiv=fmax/fbins;
   odiv=pi/obins;
   f=ceil(f./fdiv);
   a=ceil(a./odiv);
   
   if GRPHASE,
      stim=[f a ph];
      stimdim=3;
      disp('stim is freq X orientation X phase');
   else
      stim=[f a];
      stimdim=2;
      disp('stim is freq X orientation');
   end
elseif strcmp(r.h.stim(1:2),'gr'),
   disp('Gratrev: Binning in 2d Fourier space...');
   
   wx=r.s(r.sid+1,4)./2 .* sin(r.s(r.sid+1,3));
   wy=r.s(r.sid+1,4)./2 .* cos(r.s(r.sid+1,3));
   
   ph=r.s(r.sid+1,5);
   
   % stupid round about way of shuffling bins to the correct
   % region of fourier space
   maxwx=find(wx==max(wx));
   maxwy=find(wy==max(wy));
   wx=round(wx);
   wy=round(wy);
   
   wx(find(wy<0))=-wx(find(wy<0));
   wy(find(wy<0))=-wy(find(wy<0));
   flipidx=find(wx>0 & wy==0);
   wx(flipidx)=-wx(flipidx);
   wx=wx+9;
   wy=wy+1;
   wx(find(maxwy))=16;
   wy(find(maxwy))=9;
   wx(find(wx>16))=0;
   wy(find(wy>9))=0;
   
   
   if GRPHASE,
      stim=[wx wy ph];
      stimdim=3;
      disp('stim is wx X wy X phase');
   else
      stim=[wx wy];
      stimdim=2;
      disp('stim is wx X wy');
   end
elseif strcmp(r.h.stim(1:2),'nc'),
   disp('NcartRev: binning in parameter space...');
   
   disp('need to figure out binning for ncartrev');
   keyboard
   
else
   fprintf('stimtype %s unknown!\n',r.h.stim);
   keyboard
end

stimwindowsize=r.h.pix;


disp('Computing STA...');
tbins=16;
L=length(resp);
f=stim(:,1);
fbins=max(f);
a=stim(:,2);
obins=max(a);
fprintf('(fbins,obins)=(%d,%d)\n',fbins,obins);

hf=zeros(fbins,obins,tbins);
scount=zeros(fbins,obins);
for ii=1:L-tbins+1,
   if f(ii)>0 & a(ii)>0,
      %size(squeeze(hf(f,a,:)))
      %size(resp(ii:ii+tbins-1))
      hf(f(ii),a(ii),:)=squeeze(hf(f(ii),a(ii),:))+resp(ii:ii+tbins-1);
      scount(f(ii),a(ii))=scount(f(ii),a(ii))+1;
   end
end

% first order normalization... ie, don't deal with correlations
% between spatial channels
scount(find(scount<20))=1000;
hf=hf./repmat(scount,[1 1 tbins]);

if fftflag,
   % flip kernel symmetrically to make it fit standard 2d FFT bins
   hf=cat(2,zeros(fbins,obins-2,tbins),hf);
   hf(1:end-1,1:obins-1,:)=hf(1:end-1,1:obins-1,:) + ...
       flipdim(flipdim(hf(2:end,obins-1:end-1,:),1),2);
   iconside=[16 16];
   hf=reshape(hf,prod(iconside),tbins);
   if 1,
      % convert to pfft domain for compatibility with standard kernels
      cfilt=gencfilt(iconside(1),iconside(2));
      hf=hf(cfilt,:);
      kernfmt='pfft';
   else
      kernfmt='pix';
   end
else
   iconside=[fbins obins];
   hf=reshape(hf,prod(iconside),tbins);
   kernfmt='pix';
end

stephf=cumsum(hf,2);
steptime=sum(stephf,1);
maxidx=min(find(steptime==max(steptime)));
kern=stephf(:,maxidx);

% save important parms to movformat
movformat.iconside=iconside;
movformat.kernfmt=kernfmt;
movformat.daterun=date;
movformat.spacecount=prod(iconside);
movformat.maxlag=[0 tbins-1];
movformat.resampcount=1;
movformat.stimfmtcode=6;
movformat.fitfrac=fitfrac;
movformat.predfrac=predfrac;
movformat.cellid=cellid;
movformat.stimwindowsize=stimwindowsize;

disp('Done!');

save /tmp/rconline.mat hf stephf movformat kern
fprintf('hf stephf movformat kern saved to /tmp/rconline.mat\n');

figure;
titles={sprintf('%s (rawid=%d): gratrev parametric kernel',...
                r.fext,rawids),...
        sprintf('%s (rawid=%d): gratrev parametric step response',...
                r.fext,rawids)};
showkern(cat(3,hf,stephf),kernfmt,iconside,titles,1,16);


