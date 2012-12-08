% genbgpic.m
%
% matlab script to update the current photo selection. current
% photo is stored in photodir/daily.jpg
%
% currently photodir is /auto/k1/david/stuff/photos/
%
% in order to actually display the new picture, you have to run
% /usr/X11R6/bin/xv -root -geometry 1600x1200 /auto/k1/david/stuff/photos/bg/daily.jpg -quit
% in an environment with DISPLAY set to :0.0
%

% randomize
rand('state',sum(100*clock))

photodir='/home/svd/stuff/photos/bg/';

photofile=[photodir,'bginfo.mat'];
if ~exist(photofile,'file'),
   bglist={};
   bgcount=[];
else
   load(photofile);
end

newlist=jls([photodir,'*.jpg ',photodir,'*.JPG']);
newcheck=ones(length(newlist),1);
oldcheck=zeros(length(bglist),1);
for ii=1:length(newlist),
   tt=strmatch(newlist{ii},bglist);
   if tt,
      newcheck(ii)=0;
      oldcheck(tt)=1;
   end
end

bglist={bglist{find(oldcheck)},newlist{find(newcheck)}};
bgcount=[bgcount(find(oldcheck)); zeros(length(find(newcheck)),1)];

bgfrac=1./(bgcount+1);
bgfrac=cumsum(bgfrac./sum(bgfrac));
useidx=min(find(bgfrac>=rand));

%loidx=find(bgcount<=min(bgcount)+1);
%useidx=loidx(ceil(rand*length(loidx)));

bgcount(useidx)=bgcount(useidx)+1;

im=imread(bglist{useidx});
xx=max(size(im));
imsc=160/xx;
imth=imresize(im,imsc,'bilinear',0);

save(photofile,'bglist','bgcount');

unix(['\cp ',bglist{useidx},' ',photodir,'daily.jpg']);

imwrite(imth,[photodir,'thumb/daily.th.jpg']);

fprintf('Using %s (x %d).\n',bglist{useidx},bgcount(useidx));

fid = fopen([photodir,'thumb/daily.txt'],'w');
fname=strsep(bglist{useidx},'/');
fname=fname{end};
fprintf(fid,'%s',fname);
