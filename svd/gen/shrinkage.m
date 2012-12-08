% function hf=shrinkage(mH,eH,sigrat,thresh[=0])
function hf=shrinkage(mH,eH,sigrat,thresh)

if ~exist('thresh','var'),
   thresh=0;
end

if ~sigrat,
   hf=mH;
   return
end

smd=abs(mH)./(eH+eps.*(eH==0)) ./ sigrat;

if thresh,
   hf=mH.*(smd>1);
else
   smd=(1-smd.^(-2));
   smd=smd.*(smd>0);
   smd(find(isnan(smd)))=0;
   hf=mH.*smd;
end
