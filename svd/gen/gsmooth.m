% function gx=gsmooth(x,sigma,boundcond,offset)
%
% gaussian smoothing along first non-singleton dimension of x.
% filter has std dev equal to sigma (default 1)
% if sigma=[srow scol], smooth with 2D gaussian, srows across rows
%
% boundcond (default=0) 0-return gx same length as x
%                       1-return gx with tails
%                       2-return valid region of gx only
%
function gx=gsmooth(x,sigma,boundcond,offset)

if ~exist('sigma','var'),
   sigma=1;
end
if ~exist('boundcond','var'),
   boundcond=0;
end
if ~exist('offset','var'),
   offset=0;
end

%fprintf('gsmooth(...,%.2f,%d,%d):\n',sigma(1),boundcond,offset);

if abs(offset)>0,
   ll=size(x,2);
   x=[x fliplr(x(:,(ll-offset+1):ll))];
end
if length(sigma)==1,
   tt=(floor(-sigma*4):ceil(sigma*4));
   pfilt=exp(-tt.^2/(2*sigma.^2))./(sqrt(2*pi)*sigma);
   pfilt=pfilt./sum(pfilt);
   if size(x,1)>1,
      pfilt=pfilt';
   end
else
   [xx,yy]=meshgrid(floor(-sigma(2)*4):ceil(sigma(2)*4),...
                    floor(-sigma(1)*4):ceil(sigma(1)*4));
   pfilt=exp(-xx.^2/(2*sigma(2).^2)-yy.^2/(2*sigma(1).^2));
   pfilt=pfilt./sum(pfilt(:));
end

if boundcond==0,
   gx=rconv2(x,pfilt);
elseif boundcond==1,
   gx=conv2(x,pfilt,'full');
elseif boundcond==2,
   gx=conv2(x,pfilt,'valid');
end

if abs(offset)>0,
   gx=gx(:,(offset+1):end);
end
