% function xcshowkern(strf)
%
% display xc-format strf
%
% created SVD 3/18/03
%
function xcshowkern(strf)

kerncount=length(strf(:));

titles={strf(:).name};
if isfield(strf(1),'tbinms'),
   tbinms=strf(1).tbinms;
else
   tbinms=16;
end

for kernidx=1:kerncount,
   tsf=ones([size(strf(kernidx).h) kerncount])*nan;
   tsf(:,:,kernidx)=strf(kernidx).h;
   
   showkern(tsf,strf(kernidx).parms.kernfmt,...
            strf(kernidx).parms.iconside,titles,1,tbinms);
   
end
