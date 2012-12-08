% function r=npclog(fid,string,fprintfargs)
function r=npclog(fid,varargin)

s=varargin{1};
if length(varargin)>1,
   fprintfargs={varargin{2:end}};
else
   fprintfargs={};
end

fprintf(fid,[s,'\n'],fprintfargs{:});

% display to screen too for the time being
fprintf([s,'\n'],fprintfargs{:});

