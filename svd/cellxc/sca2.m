%sca.m 
%created 5/6/03 JP
%
% modified SVD 5/7/03

function smov = sca2(mov,componentCount,nlCode,doSparsen)

loadSCAmatrix;

if ~exist('doSparsen','var'),
   doSparsen=1;
end
% set in loadSCAmatrix
%if ~exist('nlCode','var'),
%   nlCode=0;
%end
%if ~exist('componentCount','var'),
%   componentCount=size(S,2);
%end

nDim = size(mov,1)*size(mov,2);
nTimes = size(mov,3);
mov = reshape(mov,nDim, nTimes);
if max(mov(1:1000))>1,
   mov=mov-128;
end

if doSparsen,
   [smov, noise] = estimateSignal(mov,S(:,1:componentCount),kappa);
   smov=full(smov);
else
   smov=S'*mov;
end

if nlCode==1,
   smov=abs(smov);
elseif nlCode==2,
   % do nl +- separation
   smov=cat(1,smov.*(smov>0), -smov.*(smov<0));
end
