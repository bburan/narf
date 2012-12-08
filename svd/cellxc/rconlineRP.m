% function [hf,stephf,movformat,kern,modspace]=rconlineRP(respfiles,stimfiles,origpix,fmt,params,varargin);
%
% respfiles - single string or cell array of strings containing
%             file name(s) of pype data files with review/natrev
%             style data
% stimfiles - downsampled imsm files with the stimuli shown in the
%             corresponding respfiles
% origpix - number of pixels in original space spanned by stimfiles
% fmt - input transformation of stimulus to use for RC
%       for WG: use fmt='parm'
%       for NR: use fmt='sfft' for phase separated fourier domain
%               use fmt='wav' for wavelet domain
% params - additional cellxcnodb parameters to override defaults
%
% returns:
% hf - STRF estimate in space X time dimensions
% stephf - step response of STRF (ie, cumsum(hf,2))
% movformat - parameters used for RC that can be used to guide
%             optinat stimulus generation
% kern - spatial only kernel. slice of stephf at maximum latency
%        (ie, greatest power before inhibition/adaptation kicks in)
%
% created SVD 5/13/02
% modified SVD 3/4/04 for data acquisition independent of celldb
%                     also generates opti-synth stimuli as part of
%                     the package
%
% BW test data:
% respfiles='/auto/k5/david/tmp/online/natsync.000';
% stimfiles='/auto/k5/david/tmp/online/natrev2004.96.index60.1.pix16';
% origpix=160;
%
function [hf,stephf,movformat,kern,modspace]=rconlineRP(respfiles,stimfiles,origpix,fmt,params,varargin);

VER='synth2.0';
BGGRAY=73;
STDGRAY=51;
FRAMEDUR=8;
OUTFRAMES=30;

if ~iscell(respfiles),
   respfiles={respfiles};
end
if ~iscell(stimfiles),
   stimfiles={stimfiles};
end
if ~exist('origpix','var'),
   disp('must specify origpix!');
   disp('syntax: rconline(respfiles,stimfiles,origpix,fmt)');
   return
end
if ~exist('fmt','var'),
   fmt='sfft';
end

if ~isempty(varargin)
  pypestimpath=varargin{1};
else
  pypestimpath='/tmp';
end

if pypestimpath(end)~='/',
   pypestimpath=[pypestimpath,'/'];
end

unix(['mkdir -p ',pypestimpath]);

tr=basename(respfiles{1});
tr=strsep(tr);
cellid=tr{1};

fprintf('%d files to process...\n',length(respfiles));
matfiles={};
for ii=1:length(respfiles),
   matfiles{ii}=processpypedata(respfiles{ii},1,0);
end

dbopen;
global BATQUEUEID
BATQUEUEID=[];

params.cellid=cellid;
params.stimfiles=stimfiles;
params.respfiles=matfiles;

resplens=zeros(length(matfiles),1);

for ii=1:length(matfiles),
   params.resplens=imfileinfo(stimfiles{ii});
   params.stimcrfs(ii)=4;
end

params.fitfrac=getparm(params,'fitfrac',0.1);
params.predfrac=getparm(params,'predfrac',0.0);
%params.=getparm(params,'',);

params.respfmtcode=getparm(params,'respfmtcode',0);
params.stimfmtcode=getparm(params,'stimfmtcode',2);

% resp load cmd. format: resploadcmd(respfile,p1,p2,...)
params.resploadcmd=getparm(params,'resploadcmd','loadpypedata');
params.resploadparms=getparm(params,'resploadparms',{'',1,0});
params.respfiltercmd=getparm(params,'respfiltercmd','');
params.respfilterparms=getparm(params,'respfilterparms',{});

% stim load cmd. should be of format
%   stimloadcmd(stimfile,startframe,endframe,p1,p2,...)
[resplen,stimpix]=imfileinfo(stimfiles{1});
stimpix=stimpix(1);
params.stimloadcmd=getparm(params,'stimloadcmd','loadimfile');
params.stimloadparms=getparm(params,'stimloadparms',{stimpix,0,stimpix,0,1});
params.stimiconside=getparm(params,'stimiconside',[stimpix stimpix]);
params.stimwindowsize=getparm(params,'stimwindowsize',stimpix);

% stim filter - eg, convert to phase-sep fourier domain.

if strcmp(fmt,'sfft'),
   params.stimfiltercmd=getparm(params,'stimfiltercmd','movphasesep');
   params.stimfilterparms=getparm(params,'stimfilterparms',{0,0,1,1,0});
   params.kernfmt=getparm(params,'kernfmt','pfft+4');
   
elseif strcmp(fmt,'space'),
   params.stimfiltercmd=getparm(params,'stimfiltercmd','');
   params.stimfilterparms=getparm(params,'stimfilterparms',{});
   params.kernfmt=getparm(params,'kernfmt','space');
   
elseif strcmp(fmt,'pfft'),
   params.stimfiltercmd=getparm(params,'stimfiltercmd','movpower');
   params.stimfilterparms=getparm(params,'stimfilterparms',{0,0,1,1,0});
   params.kernfmt=getparm(params,'kernfmt','pfft');
   
else
   disp('fmt must be: "space", "sfft" or "pfft" for NR data');
   return
end

params.minlag=getparm(params,'minlag',-8);
params.maxlag=getparm(params,'maxlag',15);
params.resampfmt=getparm(params,'resampfmt',1);
params.resampcount=getparm(params,'resampcount',20);

params.sffiltsigma=getparm(params,'sffiltsigma',5);
params.sfscount=getparm(params,'sfscount',30);
params.sfsstep=getparm(params,'sfsstep',6);

params.outfile=getparm(params,'outfile','/tmp/rconline.mat');
params.zipoutfile=getparm(params,'zipoutfile',0);

cellxcnodb(params);

load(params.outfile);
movformat=params;

nlidx=1;
hf=strf(nlidx).h;
xc=xc(:,:,nlidx,1);
predxc=predxc(:,nlidx,1);
sigfit=strf(nlidx).parms.sigfit;
sfsfit=strf(nlidx).parms.sfsfit;

stephf=cumsum(hf,2);
steptime=sum(stephf,1);
maxidx=min(find(steptime==max(steptime)));
kern=stephf(:,maxidx);

fprintf('peak of step response in bin %d. (nlidx=%d) kern=stephf(:,%d)\n',...
        maxidx,nlidx,maxidx);

if strcmp(fmt,'sfft'),
   modspace=movphasesepinv(kern);
   filt=(hanning(stimpix,'periodic') * hanning(stimpix,'periodic')');
   kern1=kern;
elseif strcmp(fmt,'space'),
   modspace=kern;
   filt=ones(stimpix,stimpix);
   kern1=kern(:);
elseif strcmp(fmt,'pfft'),
   modspace=movphasesepinv(repmat(kern,[4 1]));
   filt=(hanning(stimpix,'periodic') * hanning(stimpix,'periodic')');
   kern1=repmat(kern,[4,1]);
   movformat.stimfiltercmd='movphasesep';
end

fprintf('Generating %d good/%d bad images. origpix=%d BG=%d STD=%d\n',...
        OUTFRAMES,OUTFRAMES,origpix,BGGRAY,STDGRAY);
fprintf('loading %s\n',movformat.stimfiles{1});
mov=loadimfile(movformat.stimfiles{1})-BGGRAY;
if ~isempty(movformat.stimfiltercmd),
   pmov=feval(movformat.stimfiltercmd,mov,movformat.stimfilterparms{:});
else
   pmov=reshape(mov,stimpix*stimpix,size(mov,3));
end

predresp=(kern1'*pmov);
[tt mmax]=sort(-predresp);

peakframes=pmov(:,mmax(1:OUTFRAMES));

if strcmp(fmt,'sfft') | strcmp(fmt,'pfft'),
   phasecount=4;
   cfilt=gencfilt(stimpix,stimpix);
   cc=round(stimpix/2+1);
   c0=find(cfilt==sub2ind([stimpix,stimpix],cc,cc));
   for phaseidx=0:phasecount-1,
      kern1(c0+(stimpix^2/2*phaseidx))=0;
   end
end

kernval=abs(kern1.*(kern1>0))./std(kern1);
kernval=1-0.33*kernval.^(2);
kernval(kernval<0)=0;

peakproj=peakframes.*(kernval*ones(1,OUTFRAMES));

if strcmp(fmt,'sfft') | strcmp(fmt,'pfft'),
   [imphase,peakmov]=movphasesepinv(peakframes);
   [imphase,peakout]=movphasesepinv(peakproj);
else
   peakmov=reshape(peakframes,stimpix,stimpix,OUTFRAMES);
   peakout=reshape(peakproj,stimpix,stimpix,OUTFRAMES);
end

predbest=(kern1'*peakframes)';
predout=(kern1'*peakproj)';

biggood=zeros(origpix,origpix,OUTFRAMES);
bigbad=zeros(origpix,origpix,OUTFRAMES);

figure
for ii=1:OUTFRAMES,
   % make sure edges are still hanning'ed
   peakout(:,:,ii)=peakout(:,:,ii).*filt;
   
   subplot(2,3,1);
   imagesc(peakmov(:,:,ii)+BGGRAY,[0 255]);
   axis image; axis off;
   
   subplot(2,3,2);
   imagesc(peakmov(:,:,ii)-peakout(:,:,ii)+BGGRAY,[0 255]);
   axis image; axis off;
   
   subplot(2,3,3);
   imagesc(peakout(:,:,ii)+BGGRAY,[0 255]);
   axis image; axis off;
   
   biggood(:,:,ii)=imresize(peakmov(:,:,ii)-peakout(:,:,ii),...
                            [origpix origpix],'bilinear');
   bigbad(:,:,ii)=imresize(peakout(:,:,ii),...
                           [origpix origpix],'bilinear');
   tt=biggood((1:round(origpix/2))+round(origpix/2),...
              (1:round(origpix/2))+round(origpix/2),ii);
   tstd=sqrt(mean(tt(:).^2));
   biggood(:,:,ii)=biggood(:,:,ii)./tstd.*STDGRAY+BGGRAY;
   tt=bigbad((1:round(origpix/2))+round(origpix/2),...
              (1:round(origpix/2))+round(origpix/2),ii);
   tstd=sqrt(mean(tt(:).^2));
   bigbad(:,:,ii)=bigbad(:,:,ii)./tstd.*STDGRAY+BGGRAY;
   
   subplot(2,3,4);
   imagesc(biggood(:,:,ii)+bigbad(:,:,ii)-BGGRAY,[0 255]);
   axis image; axis off;
   
   subplot(2,3,5);
   imagesc(biggood(:,:,ii),[0 255]);
   axis image; axis off;
   
   subplot(2,3,6);
   imagesc(bigbad(:,:,ii),[0 255]);
   axis image; axis off;
   
   colormap(gray);
   drawnow
end

fprintf('Saving opti stim to: %s\n',pypestimpath);

indexfile='index';
fid=fopen([pypestimpath,indexfile],'w');
fprintf(fid,'# stim %s\n# pix %d\n# framecount %d\n',...
        VER,origpix,FRAMEDUR*OUTFRAMES*20);
fprintf(fid,'# respfile %s\n# stimfile %s\n',...
        movformat.respfiles{1},movformat.stimfiles{1});

for ii=1:OUTFRAMES,
   fn=sprintf('image_%.4d.pgm',ii*2-1);

   % bw mar 18 2004. changed to use pgmWritebw so that pgms don't contain any zeros
   tmp = max(min(biggood(:,:,ii),255),1);
   if ~isempty(find(tmp(:)==0))
     'FUCK'
   end
   
   pgmWritebw(tmp,[pypestimpath,fn],[0 255]);   
   %pgmWrite(biggood(:,:,ii),[pypestimpath,fn],[0 255]);

   fprintf(fid,'%s %d\n',fn,FRAMEDUR);
   
   fn=sprintf('image_%.4d.pgm',ii*2);

   % bw mar 18 2004. changed to use pgmWritebw so that pgms don't contain any zeros
   tmp = max(min(bigbad(:,:,ii),255),1);
   pgmWritebw(tmp,[pypestimpath,fn],[0 255]);
   if ~isempty(find(tmp(:)==0))
     'FUCK'
   end
   %pgmWrite(bigbad(:,:,ii),[pypestimpath,fn],[0 255]);

   fprintf(fid,'%s %d\n',fn,FRAMEDUR);
end

rand('state',0);

for jj=1:9
  r=randperm(OUTFRAMES*2);
  for ii=1:OUTFRAMES*2
    fn=sprintf('image_%.4d.pgm',r(ii));
    fprintf(fid,'%s %d\n',fn,FRAMEDUR);
  end
end

fclose(fid);

save /tmp/rconline.mat hf stephf movformat kern modspace
fprintf('hf stephf movformat kern modspace saved to /tmp/rconline.mat\n');



%respfile='/auto/k5/david/tmp/online/r0284.rev.002'
%stimfile='/auto/data/archive/stimarchive/REESE16/imsm/natrev_malik-96.index60.1.pix.imsm'


