% smap = rcspot(pypefile, ..opts..)
%
%  pypefile = name of pype data file ([] for last plotted file)
%  optional parameters:
%    'fit', 'notfit' -> set rf fitting.
%    'lat' -> set latency (in ms)
%    'binsize' -> ms per time bin (default 16)
%    'lagcount' -> number of time bins 0 thru lagcount-1 (default 10)
%    'tcorr' -> correction for sync pulse lag (in ms)
%
% created SVD 2/8/03 - ripped off jamie's spotmap.m
%
function smap = rcspot(pypefile, varargin)

l = jls(pypefile);
length(l)

if length(l) > 1
  smap = [];
  for i=1:length(l)
    clf;
    s = spotmap(char(l(i)), varargin{:});
    s.file = char(l(i));
    smap = [smap s];
    getframe;
    fullpage;
    print -dpsc
  end
  return
end

narg = 1;
latency = 0;
winsize=0;
while narg <= length(varargin)
  switch varargin{narg}
   case 'lat'
    latguess = varargin{narg + 1};
    narg = narg + 1;
   case 'tcorr'
    tcorr = varargin{narg + 1};
    narg = narg + 1;
   case 'lagcount'
    lagcount = varargin{narg + 1};
    narg = narg + 1;
   case 'binsize'
    binsize = varargin{narg + 1};
    narg = narg + 1;
   otherwise
    error(sprintf('unknown option: %s', varargin{narg}));
  end
  narg = narg + 1;
end

% dump contents of pype file to text fmt for reading into matlab
tp=tempname;
c = sprintf(['pypenv dumpspotmap %s > %s'],pypefile,tp);
unix(c);
f = basename(pypefile);

s=load(tp);
delete(tp);

% figure out some basic specs of the spotmap
trialcount=max(s(:,1));
spotidx=find(s(:,4)<2);
xrange=unique(s(spotidx,5));
yrange=unique(s(spotidx,6));

if ~exist('tcorr','var'), % ms to correct spot times
   tcorr=0;
end
if ~exist('binsize','var'),
   binsize=16;
end
if ~exist('lagcount','var'),
   lagcount=10;
end
fprintf('binsize=%d ms\nlagcount=%d\nsync correction=%d ms\n',...
        binsize,lagcount,tcorr);

if exist('latguess','var'),
   latguess=floor(latguess/binsize)+1;
   fprintf('latency fixed to %d ms\n',latguess*binsize);
else
   latguess=0;
end

totallen=lagcount*binsize;

xcount=length(xrange);
ycount=length(yrange);
hfhi=zeros(ycount,xcount,lagcount);
hflo=zeros(ycount,xcount,lagcount);
hfall=zeros(ycount,xcount,lagcount);
nhi=zeros(ycount,xcount);
nlo=zeros(ycount,xcount);
nall=zeros(ycount,xcount);
spikelist=[];
for tridx=1:trialcount,
   events=s(find(s(:,1)==tridx & s(:,4)<2),2:end);
   spiketimes=s(find(s(:,1)==tridx & s(:,4)==2),2);
   spikelist=[spikelist; diff(spiketimes)];
   
   events(:,1:2)=events(:,1:2)+tcorr;
   
   for ii=1:size(events,1),
      iy=find(yrange==events(ii,5));
      ix=find(xrange==events(ii,4));
      
      validspikes=find(spiketimes>=events(ii,1) & ...
                       spiketimes<events(ii,1)+totallen);
      if ~isempty(validspikes),
         nextspikes=spiketimes(validspikes)-events(ii,1);
         histspikes=histc(nextspikes,0:binsize:totallen+1);
         histspikes=reshape(histspikes(1:end-1),1,1,lagcount);
      else
         histspikes=zeros(1,1,lagcount);
      end
      if events(ii,3)==0,
         hflo(iy,ix,:)=hflo(iy,ix,:)+histspikes;
         nlo(iy,ix)=nlo(iy,ix)+1;
      else
         hfhi(iy,ix,:)=hfhi(iy,ix,:)+histspikes;
         nhi(iy,ix)=nhi(iy,ix)+1;
      end
      
      hfall(iy,ix,:)=hfall(iy,ix,:)+histspikes;
      nall(iy,ix)=nall(iy,ix)+1;
   end
end

meanrate=1000./mean(spikelist);

% normalize and convert to Hz
nlo(find(nlo==0))=1000;
nhi(find(nhi==0))=1000;
nall(find(nall==0))=1000;

hflo=hflo./repmat(nlo,[1 1 lagcount])./binsize*1000;
hfhi=hfhi./repmat(nhi,[1 1 lagcount])./binsize*1000;
hfall=hfall./repmat(nall,[1 1 lagcount])./binsize*1000;

pixcount=ycount*xcount;
hf=cat(4,hflo,hfhi,hfall);
hf=reshape(hf,pixcount,lagcount,3);
disppypefile=pypefile;
disppypefile(find(disppypefile=='_'))='-';
titles={sprintf('%s OFF',disppypefile),sprintf('%s ON',disppypefile),...
        sprintf('%s BOTH',disppypefile)};
names={'OFF','ON','BOTH'};
figure(1);
showkern(hf,'space',[ycount,xcount],titles,1,binsize);
set(gcf,'PaperOrientation','landscape','PaperPosition',[0.25 0.25 10.5 8]);
colormap(gray);
set(get(gcf,'Children'),'YDir','normal');

figure(2)
clf
colormap(redblue(0.6));

rowcount=3;
for ii=1:rowcount,
   
   % plot summary kernel stuff - take svd of kernel to get spatial
   % and temporal marginals
   tsf=hf(:,:,ii);
   
   %tsf=tsf-repmat(mean(tsf,2),[1 size(tsf,2)]);
   
   tsftime=sum(tsf,1)./(xcount*ycount);
   if latguess>0,
      lat=latguess;
   else
      lat=min(find((tsftime)==max((tsftime(1:min([lagcount 9]))))));
   end
   
   sp=tsf(:,lat);
   
   ir=tsftime;
   
   spacebincount=length(sp);
   tsf0=reshape(sp,ycount,xcount);
   
   % plot spatial signature
   h=subplot(rowcount,3,ii*3-2);
   pcsummax=max(abs(tsf0(:)));
   if 0 & ii==1, % ie OFF
      imagesc(xrange,yrange,-tsf0,[-pcsummax pcsummax]);
   else
      imagesc(xrange,yrange,tsf0,[-pcsummax pcsummax]);
      %imagesc(xrange,yrange,tsf0,[0 pcsummax]);
   end
   axis image
   axis xy 
   colorbar('horiz');
   
   q=zeros(3,1);
   q2=zeros(3,1);
   if 1, 
      xx=-2:2;
      gsf=exp(-(xx./0.5).^2/2);
      gsf=gsf./sum(gsf(:));
      gor=exp(-(xx./0.5).^2/2);
      gor=(gor./sum(gor(:)))';
      
      astd=std(tsf0(:));
      ta=conv2(tsf0,gsf,'same');
      ta=conv2(ta,gor,'same');
      tas=size(ta);
      
      hold on
      c=contour(xrange,yrange,ta,[ astd*1.0   astd*1.0],'k--');
      if length(c)>0,
	 [cx,cy]=contour2points(c);
	 q = fit_circ('fit',cx,cy);
      end
      
      c=contour(xrange,yrange,ta,[ astd*2.0   astd*2.0],'k-');
      if length(c)>0,
         [cx,cy]=contour2points(c);
         q2 = fit_circ('fit',cx,cy);
	 plot(q2(1),q2(2),'kx');
      end
      hold off
   end
   
   %convert fit rad to diam:
   q(3)=q(3)*2;
   q2(3)=q2(3)*2;
   
   if ii==1,
      title([disppypefile,' ',names{ii}]);
   else
      title(names{ii});
   end
   xlabel('x pos (pix)');
   ylabel('y pos (pix)');
   
   % plot temporal response
   subplot(rowcount,3,ii*3-1);
   tt=(1:length(ir))*binsize-binsize;
   h=plot(tt,ir,'k-','LineWidth',1);
   hold on
   plot(tt,meanrate.*ones(1,length(ir)),'k--');
   hold off
   axis([0 max(tt) min(ir)*0.8 max(abs(ir)*1.1)]);
   axis square
   title('temporal response');
   xlabel('time lag (ms)');
   ylabel('resp (Hz)');
   
   subplot(rowcount,3,ii*3);
   axis off
   text(0,1,{sprintf('{\\bf%s}',names{ii}),...
             sprintf('Big: (x,y)=(%.0f,%.0f) D=%.0f',q),...
             sprintf('Small: (x,y)=(%.0f,%.0f) D=%.0f',q2),...
             sprintf('Lat max: %d ms (bin %d)',(lat-1)*binsize,lat),...
            },'VerticalAlign','top');
   smap(ii).file=pypefile;
   smap(ii).spot=names{ii};
   smap(ii).fit1=q;
   smap(ii).fit2=q2;
   smap(ii).lat=lat;
end

set(gcf,'PaperOrientation','portrait','PaperPosition',[0.25 0.25 8 10.5]);

return

keyboard

off = s(find(s(:,1) == 0), 2:end);
on = s(find(s(:,1) == 1), 2:end);
all = s(find(s(:,1) == 2), 2:end);

smap.xmin = min(s(:,2));
smap.xmax = max(s(:,2));
smap.ymin = min(s(:,3));
smap.ymax = max(s(:,3));

smap.latency = latency;
smap.winsize = winsize;

x0=off(:,1); y0=off(:,2); z0=off(:,3);
[x, y, z1, cx, cy] = priv_xyz(x0, y0, z0, smooth, 1, contourval, 0);
x0=on(:,1); y0=on(:,2); z0=on(:,3);
[x, y, z2, cx, cy] = priv_xyz(x0, y0, z0, smooth, 1, contourval, 0);
x0=all(:,1); y0=all(:,2); z0=all(:,3);
[x, y, z3, cx, cy] = priv_xyz(x0, y0, z0, smooth, 1, contourval, 0);
z = [z1(:) z2(:) z3(:)];
zmax = max(z(:));



%%%%%%%% OFF %%%%%%%%%%
x0=off(:,1);
y0=off(:,2);
z0=off(:,3);

[x, y, z, cx, cy] = priv_xyz(x0, y0, z0, smooth, 1, contourval, zmax);
q = fit_circ('fit',cx,cy);
smap.c_off.x = q(1);
smap.c_off.y = q(2);
smap.c_off.r = q(3);
if ~noplot
  if isempty(plotpos)
    subplot(4,2,1);
  else
    mysubplot(plotpos(1:4));
  end
  priv_xyz(x0, y0, z0, smooth, 0, contourval, zmax);
  %colormap(blueyellow);
  colormap(hotcold(1));
  [a,b]=meshgrid(x0,y0);
  colorbar;
  %title({sprintf('{\\bf%s}', f), 'Off Response (sp/sec)'});
  title(sprintf('{\\bf%s}', f));
  ylabel('OFF');
  hold on;
  set(plot(smap.c_off.x, smap.c_off.y, 'g+'), ...
      'MarkerSize', 10);
  %circle([smap.c_off.x smap.c_off.y], smap.c_off.r, 'w-');
  plot(a, b, 'k.');
  hold off;
  %squareup;
  axis image
end


%%%%%%%% ON %%%%%%%%%%
x0=on(:,1);
y0=on(:,2);
z0=on(:,3);

[x, y, z, cx, cy] = priv_xyz(x0, y0, z0, smooth, 1, contourval, zmax);
q = fit_circ('fit',cx,cy);
smap.c_on.x = q(1);
smap.c_on.y = q(2);
smap.c_on.r = q(3);

if ~noplot
  if isempty(plotpos)
    subplot(4,2,3);
  else
    mysubplot(plotpos(5:8));
  end
  priv_xyz(x0, y0, z0, smooth, 0, contourval, zmax);
  %colormap(blueyellow);
  colormap(hotcold(1));
  [a,b]=meshgrid(x0,y0);
  colorbar;
  %title({sprintf('{\\bf%s}', f), 'On Response (sp/sec)'});
  ylabel('ON');
  hold on;
  set(plot(smap.c_on.x, smap.c_on.y, 'g+'), ...
      'MarkerSize', 10);
  %circle([smap.c_on.x smap.c_on.y], smap.c_on.r, 'w-');
  plot(a, b, 'k.');
  hold off;
  %squareup;
  axis image
end


%%%%%%%% ALL %%%%%%%%%%
x0=all(:,1);
y0=all(:,2);
z0=all(:,3);

[x, y, z, cx, cy] = priv_xyz(x0, y0, z0, smooth, 1, contourval, zmax);
q = fit_circ('fit',cx,cy);
smap.c_all.x = q(1);
smap.c_all.y = q(2);
smap.c_all.r = q(3);

if ~noplot
  if isempty(plotpos)
    subplot(4,2,5);
  else
    mysubplot(plotpos(9:12));
  end
  priv_xyz(x0, y0, z0, smooth, 0, contourval, zmax);
  %colormap(blueyellow);
  colormap(hotcold(1));
  [a,b]=meshgrid(x0,y0);
  colorbar;
  %title({sprintf('{\\bf%s}', f), 'All Response (sp/sec)'});
  ylabel('COMP');
  hold on;
  set(plot(smap.c_all.x, smap.c_all.y, 'g+'), ...
      'MarkerSize', 10);
  %circle([smap.c_all.x smap.c_all.y], smap.c_all.r, 'w-');
  plot(a, b, 'k.');
  hold off;
  %squareup;
  axis image
end

if ~noplot & isempty(plotpos)
  subplot(4,2,2);
  cla;
  axis off;
  text(0, 0.5, { ...
      sprintf('{\\bf OFF (pixels)}'), ...
      sprintf('x: %.0f', smap.c_off.x), ...
      sprintf('y: %.0f', smap.c_off.y), ...
      sprintf('d: %.0f', (smap.c_off.x^2+smap.c_off.y^2)^0.5), ...
      sprintf('s: %.0f', smap.c_off.r), ...
      sprintf('cont=%.0f%%max', contourval*100), ...
      sprintf('smooth=%.1f', smooth)});
  subplot(4,2,4); axis off;
  cla;
  text(0, 0.5, { ...
      sprintf('{\\bf ON (pixels)}'), ...
      sprintf('x: %.0f', smap.c_on.x), ...
      sprintf('y: %.0f', smap.c_on.y), ...
      sprintf('d: %.0f', (smap.c_on.x^2+smap.c_on.y^2)^0.5), ...
      sprintf('s: %.0f', smap.c_on.r), ...
      sprintf('cont=%.0f%%max', contourval*100), ...
      sprintf('smooth=%.1f', smooth)});
  subplot(4,2,6); axis off;
  cla;
  text(0, 0.5, { ...
      sprintf('{\\bf ALL (pixels)}'), ...
      sprintf('x: %.0f', smap.c_all.x), ...
      sprintf('y: %.0f', smap.c_all.y), ...
      sprintf('d: %.0f', (smap.c_all.x^2+smap.c_all.y^2)^0.5), ...
      sprintf('s: %.0f', smap.c_all.r), ...
      sprintf('cont=%.0f%%max', contourval*100), ...
      sprintf('smooth=%.1f', smooth)});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [newx, newy, zi, cx, cy] = priv_xyz(x, y, z, ...
					     smooth, noplot, p, zmax)

xv = unique(sort(x));
yv = unique(sort(y));

newx = (min(xv):max(diff(xv)):max(xv));
newy = (min(yv):max(diff(yv)):max(yv));

[xi,yi]=meshgrid(newx, newy);
zi = griddata(x,y,z,xi,yi);

if smooth > 0
  zi = smooth2d(zi, 0, 0, smooth);
end

if ~noplot
  contourf(xi, yi, zi, 10);
  caxis([0 zmax]);
end

z = zi - min(zi(:));
z = z ./ max(z(:));
c = contourc(newx, newy, z, [p p]);
[cx,cy]=contour2points(c);
if ~noplot
  hold on;
  l = plot(cx, cy, 'ro');
  set(l, 'MarkerSize', 5, 'MarkerFaceColor', 'y', 'MarkerEdgeCOlor', 'y');
  hold off;
end



