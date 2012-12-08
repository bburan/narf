%
% LOAD ALL TARGETS AND OTHER PATCHES. CONVERT TO FORMAT USED FOR RC
%
% set parameters for conversion -- hacked out of free2mov03.m
movformat.cellid=params.cellid;
sql=['SELECT * FROM gCellMaster WHERE cellid="',movformat.cellid,'"'];
celldata=mysql(sql);

if isempty(celldata.rfsize) | celldata.rfsize<=0,
   movformat.stimwindowsize=32;
   movformat.rfsize=32;
elseif celldata.rfsize>200,
   movformat.stimwindowsize=200;
   movformat.rfsize=celldata.rfsize;
else
   movformat.stimwindowsize=celldata.rfsize;
   movformat.rfsize=celldata.rfsize;
end
if isempty(celldata.xoffset),
   movformat.offx=0;
   movformat.offy=0;
else
   movformat.offx=celldata.xoffset;
   movformat.offy=celldata.yoffset;
end

fprintf('loading patch data from %s...\n',params.respfiles{1});

% last file is longest, most likely to have all 4 targets noted
z=load(params.respfiles{end});
targlist=z.targlist;

fprintf('Targets id''d:');
for ii=1:length(targlist),
   fprintf(' %d',targlist(ii));
end
fprintf('\n');

disp('generating patches for matched filter analysis...');

if params.stimfmtcode==5 | params.stimfmtcode==4,
   iconside=[movformat.stimwindowsize movformat.stimwindowsize];
else
   iconside=strsep(cellfiledata(1).stimiconside,',');
   iconside=cat(2,iconside{:});
end

imagepix=size(z.patches,1);
patchcount=size(z.patches,3);
patchlist=1:patchcount;
if isfield(z,'alpha'),
   alphamask=z.alpha;
else
   disp('punting on alpha mask!');
   alphamask=ones(imagepix);
end
bgpix=z.movformat.blankpix;
targpatches=ones([iconside(:)',patchcount]) .* alphamask(1,1);

patchpix=movformat.stimwindowsize;
bigpatches=zeros(patchpix,patchpix,patchcount);

for ii=1:patchcount,
   
   % insert patch in window sized to rf.
   if imagepix>=patchpix,
      m=floor((imagepix+1)/2);
      timage=alphamask .* z.patches(:,:,ii) + ...
             (1-alphamask) .* bgpix;
      timage=timage(round(m-patchpix/2+1:m+patchpix/2),...
                    round(m-patchpix/2+1:m+patchpix/2));
   else
      m=floor((patchpix+1)/2);
      timage=ones(patchpix) .* z.patches(1,1,ii);
      timage(m-imagepix/2+1:m+imagepix/2,...
             m-imagepix/2+1:m+imagepix/2)=...
         alphamask .* z.patches(:,:,ii) + ...
         (1-alphamask) .* bgpix;
   end
   
   fmask=ones(patchpix);
   croprange=[1 1 patchpix patchpix];
   targpatches(:,:,ii)=movresize(timage,iconside(1),fmask,croprange,0,1);
   bigpatches(:,:,ii)=movresize(timage,patchpix,fmask,croprange,0,1);
end


% option: use smaller window if alpha-ed patch region is less than
% rf diameter
% find non-bg region of patches
if params.stimloadparms{1}==0,
   tt=abs(targpatches(:,:,1)-targpatches(1))>0.1;
   maskwid=sum(sum(tt,1)>0);
   totalwid=size(targpatches,1);
   if maskwid<totalwid/2,
      maskwid=ceil(totalwid/2);
   end
   params.stimloadparms{1}=round(params.stimloadparms{3}.*totalwid/maskwid);
   fprintf('set scale_to_pix to %d\n',params.stimloadparms{1});
end


% simulate what happened to stimulus during RC:

% does it need to be resized again???
if params.stimfmtcode==5,
   wavside=strsep(cellfiledata(1).stimiconside,',');
   wavside=cat(2,wavside{:});
   spacebincount=wavside(1);
   obincount=wavside(3);
   sbincount=wavside(4);
   targpatches=mov2wav(targpatches,spacebincount,obincount,sbincount);
   targpatches=reshape(targpatches,prod(wavside),patchcount);
elseif params.stimfmtcode==4,
   targpatches=eye(patchcount);
elseif length(params.stimloadparms)>0,
   % hacked from loadimfile...
   [fmask,crop]=movfmask(...
      iconside(1),params.stimloadparms{2}./params.stimloadparms{1},...
      params.stimloadparms{3}*iconside(1)./params.stimloadparms{1});
   
   targpatches=movresize(targpatches,params.stimloadparms{3},fmask,crop,0,1);
end

% and do the appropriate stim filter:
if ~isempty(params.stimfiltercmd),
   fpatches=feval(params.stimfiltercmd,targpatches,params.stimfilterparms{:});
else
   fpatches=targpatches;
end
iconside=size(fpatches);
iconside=iconside(1:(end-1));
if length(iconside)>=2,
   % reshape all spatial dims into one
   fpatches=reshape(fpatches,prod(iconside),...
                    size(fpatches,length(iconside)+1));
end
mpatches=[mean(fpatches,2),fpatches(:,targlist)];

% clear big variables to save memory
clear data z

