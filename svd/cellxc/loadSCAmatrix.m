% loadSCAmatrix
%
% this is a script that can be called from other programs to get
% the current "best" S and kappa values for sca reverse correlation stuff.

kappa=2.5;

load /auto/k1/jess/ssa/data/ds4XvanhTrain400_10.mat


if ~exist('nlCode','var'),
   if exist('params','var') & length(params.stimfilterparms)>=2,
      nlCode=params.stimfilterparms{2};
   else
      nlCode=0;
   end
end
if ~exist('componentCount','var'),
   componentCount=size(S,2);
end

S=S(:,1:componentCount);
if nlCode==2,
   S=cat(2,S,-S);
end
