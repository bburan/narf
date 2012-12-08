function h=xticks(varargin)
%XTICKS  enables rich formatting of tickmarks
%  XTICKS(X,LABEL) creates tickmarks at positions X
%  X should be a vector, 
%  LABEL may be a cell array of the same size, or char matrix with lenght(X) rows
%  any TeX notation can be used.
% 
%  XTICKS(X) creates labels automaticly
%
%  XTICKS(...,'Property',Value,...) passes TEXT properties
%       
%  H=XTICKS(...) returns handles of TEXT objects
%        
%   Note: 
%  Tickmarks created this way are sensitive to YLIM: if it is changed, tickmarks 
%  may move. Calling XTICKS without any arguments after changing YLIM restores them.
%
%   Example:
%  XTICKS([0 pi 2*pi],{'Zero','\pi','2\pi'},'rotation',45,'fontsize',14) 

% 6/4/01 ashcherbina@ucsd.edu

fig = gcbo;     %if executed as a callback
if isempty(fig), fig=gcf;end

u = findobj(fig, 'Tag','SpecialXTick');
if ~isempty(u)
   axs = get(u,'parent');
   axs=cat(1,axs{:});
else
   axs=[];
end

% if moving, just update all ticks
if (isempty(varargin) & ~isempty(u))
   if length(u)==0 set(gcbo,'ResizeFcn',[]);return; end;
   for ax=unique(axs)'
      tm=u(axs==ax); %ticks belonging to this axis
      tm=sort(tm);
      pos=get_positions(ax);
      set(tm,{'position'},pos);
   end
   return;
end
% get ticks and labels
ax=gca;
xpos=get(ax,'Xtick');
str=cellstr(get(ax,'XtickLabel'));
if ~isempty(varargin) & ~isstr(varargin{1}), 
   xpos=varargin{1};
   set(ax,'Xtick',xpos);
   str=cellstr(get(ax,'XtickLabel'));
   if length(varargin)>1 & length(varargin{2})>1,
      str=varargin{2};
      if ~iscell(str),  str=cellstr(str);end;
      varargin={varargin{3:end}};
   else 
      varargin={varargin{2:end}};
   end
   xpos=xpos(:);
   str={str{:}}';
   visible=(xpos>=min(xlim) & xpos<=max(xlim));
   xpos=xpos(visible);
   str={str{visible}}';
   set(ax,'Xtick',xpos);

end

if ~isempty(axs)
   tm=u(axs==ax); %ticks belonging to the current axis
   if ~isempty(tm)
      % ticks exist - delete them
      delete(tm);  
   end
end

% fix limits & tickmarks
ylim(ylim);
xlim(xlim);

set(ax,'XtickMode','MANUAL','XtickLabelMode','MANUAL','Xticklabel',[]);

pos=get_positions(ax);
N=length(pos);

%create text objects
hh=text(zeros(N,1),zeros(N,1),str);
hh=sort(hh);
set(hh,{'position'},pos,...
   {'string'},str,...
        'horizontalalignment','center',...
   'clipping','off',...
   'Tag','SpecialXTick',...
   varargin{:});
set(gcf,'ResizeFcn','xticks');  % establish figure callback
if nargout>0, h=hh;end;

%%%%%%%%%%%%

function pos=get_positions(ax)
YSPACE=10;       %spacing between the axis and the labels, points
unit=get(ax,'units');
set(ax,'units','points');
rect=get(ax,'Position'); % actual_pos(ax);
set(ax,'units',unit);

xpos=get(ax,'Xtick');
yl=get(ax,'Ylim');
y=yl(1)-YSPACE/rect(4)*diff(yl);
ypos=repmat(y,1,length(xpos));
zpos=0*ypos;
pos=num2cell([xpos(:),ypos(:),zpos(:)],2);

