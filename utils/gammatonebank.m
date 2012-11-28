function [bms,envs,instfs,delays]=gammatonebank(x,lowcf,highcf,numchans,fs,align)

% GAMMATONEBANK Bank of fourth order gammatone filters
% 
% Implements a bank of 4th order gammatone filters.
% 
% SYNTAX
% 
%     [BMS,ENVS,INSTFS]=GAMMATONEBANK(X,LOWCF,HIGHCF,NUMCHANS,FS)
%     [BMS,ENVS,INSTFS]=GAMMATONEBANK(X,LOWCF,HIGHCF,NUMCHANS,FS,ALIGN)
%     [BMS,ENVS,INSTFS,DELAY]=GAMMATONEBANK(...)
% 
% 
% DESCRIPTION
% 
% [BMS,ENVS,INSTFS]=GAMMATONEBANK(X,LOWCF,HIGHCF,NUMCHANS,FS) returns
% an array of NUMCHANS gammatone filterbank outputs in response to
% vector X (sampled at frequency FS). Centre frequencies are equally
% spaced on the ERB-rate scale, between LOWCF and HIGHCF. Instantaneous
% envelopes ENVS and instantaneous frequencies INSTFS may also be
% returned.
% 
% [BMS,ENVS,INSTFS]=GAMMATONEBANK(X,LOWCF,HIGHCF,NUMCHANS,FS,ALIGN)
% allows phase alignment to be applied. With ALIGN=false, no alignment
% is applied (default). With ALIGN=true, fine structure and envelope
% alignment is applied so that the impulse response peak occurs at t=0.
% 
% [BMS,ENVS,INSTFS,DELAY]=GAMMATONEBANK(...) returns the delay (in
% samples) removed by the phase alignment of each gammatone filter,
% i.e. DELAY=0 if ALIGN=FALSE.
% 
% Adapted from code by Guy Brown, University of Sheffield and Martin
% Cooke.
% 
% See also GAMMATONE

if nargin<5
    error('Not enough input arguments')
end
if nargin<6
    align = false;
end
if numel(x)~=max(size(x))
    error('x must be a vector')
end
x = reshape(x,1,length(x));

bms = zeros(numchans,length(x));
envs = zeros(numchans,length(x));
instfs = zeros(numchans,length(x));
delays = zeros(numchans,1);

cfs=MakeErbCFs(lowcf,highcf,numchans);
for c=1:numchans
  cf=cfs(c);
  [bm,env,instf,delay]=gammatone(x,fs,cf,align);
  bms(c,:) = bm;
  envs(c,:) = env;
  instfs(c,:) = instf;
  delays(c) = delay;
end