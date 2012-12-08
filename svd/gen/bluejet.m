% function pal=bluejet(stepsize);
%
% colormap based on jet that emphasizes the blue
% 
function pal=bluejet(stepsize);

if ~exist('stepsize','var'),
   stepsize=3;
end

samples=fliplr(round(31:-stepsize:1));
boundary=32-length(samples);

pal=jet;
pal(boundary:31,:)=pal(samples,:);
pal(1:(boundary-1),:)=repmat([0 0 0.5],[boundary-1,1]);


