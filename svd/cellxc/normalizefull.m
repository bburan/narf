% function [H,lambda,pctvar]=normalizefull(SR,sSA,tSA,nfactors,lambdaomag,topSR,smoothtime)
% 
% created SVD 4/12/02 - hacked from normalize.m
% modified SVD 4/15/02 - scale lambdas according to eigenvalues
% 
% Inputs (basically the output from movxc.m or summed output from
%        multiple runs of it):
%      N - number of spatial dimensions
%      U - total number of time lags (eg, length([maxlag(1):maxlag(2)]))
%      M - number of response sets (or response dimensions)
% SR - raw spike-triggered average (correlation between stim and
%      resp) (N x U x M)
% sSAfull - spatial autocorrelation matrix (NU x NU X M and use
%       the same one for each response set--speeds things up)
% nfactors - number of different regularization parameters to test (ranging
%            (10^0 to 10^-lambdaomag)
% lambdaomag - orders of magnitude to probe lambda (the
%              regularization parameter). or (if more than one
%              value, specific regulartization parameters)
%
% Outputs:
% H - decorrelated kernel estimate(s), N x U x nfactors x M matrix
% lambda - regularization parameters corresponding to each estimate
%          of H
% 
function [H,lambda,pctvar]=normalizefull(SR,sSA,nfactors,lambdaomag,topSR,smoothtime);

disp('normalizefull.m: reshaping for full space-time decorrelation');

s=[size(SR) 1 1];
SR=reshape(SR,s(1)*s(2),1,s(3),s(4));
[H,lambda,pctvar]=normalizereg(SR,sSA,1,nfactors,lambdaomag);
H=reshape(H,[s(1:2) size(H,3) size(H,4) size(H,5)]);



