function [bm,env,instf,delay]=gammatone(x,fs,cf,align)

% GAMMATONE Fourth order gammatone filter
% 
% Implements the 4th order gammatone filter using the approach described
% in Martin Cooke's PhD thesis (impulse invariant transform)
% 
% SYNTAX
% 
%     [BM,ENV,INSTF]=GAMMATONE(X,FS,CF)
%     [BM,ENV,INSTF]=GAMMATONE(X,FS,CF,ALIGN)
% 
% 
% DESCRIPTION
% 
% [BM,ENV,INSTF]=GAMMATONE(X,FS,CF) calculates the simulated basilar
% membrane displacement BM, as estimated by the gammatone filterbank,
% for input vector X, sampled at frequency FS and at centre frequency
% CF. Instantaneous envelope ENV and instantaneous frequency INSTF
% may also be returned.
% 
% [BM,ENV,INSTF]=GAMMATONE(X,FS,CF,ALIGN) allows phase alignment
% to be applied. With ALIGN=false, no alignment is applied (default).
% With ALIGN=true, fine structure and envelope alignment is applied
% so that the impulse response peak occurs at t=0.
% 
% [BMS,ENVS,INSTFS,DELAY]=GAMMATONEBANK(...) returns the delay (in
% samples) removed by the phase alignment of the filter, i.e. DELAY=0
% if ALIGN=FALSE.
% 
% Adapted from code by Guy Brown, University of Sheffield, and Martin
% Cooke.
% 
% See also GAMMATONEBANK

if nargin<3
    error('Not enough input arguments')
end
if nargin<4
    align = false;
end
if ~islogical(align)
    error('Phase correction parameter must be logical')
end
if numel(x)~=max(size(x))
    error('x must be a vector')
end
x = reshape(x,1,length(x));

if align
    B=1.019*2*pi*erb(cf);
    envelopecomptime = 3/B;
else
    envelopecomptime = 0;
end
phasealign=-2*pi*cf*envelopecomptime;
phasealign=mod(phasealign,2*pi);
phasealign=phasealign/(2*pi*cf);
shift=envelopecomptime;
intshift=round(shift*fs);

bw=1.019*erb(cf); % bandwidth
wcf=2*pi*cf; % radian frequency
tpt=(2*pi)/fs;
a=exp(-bw*tpt);
gain=((bw*tpt)^4)/6; % based on integral of impulse response

x = [x zeros(1,intshift)];
kT=(0:length(x)-1)/fs;

q=exp(1i.*(-wcf.*kT)).*x; % shift down to d.c.
p=filter([1 0],[1 -4*a 6*a^2 -4*a^3 a^4],q); % filter: part 1
u=filter([1 4*a 4*a^2 0],[1 0],p); % filter: part 2
bm=gain*real(exp(1i*wcf*(kT(intshift+1:end)+phasealign)).*u(intshift+1:end)); % shift up in frequency
env = gain*abs(u(intshift+1:end));
instf=real(cf+[diff(unwrap(angle(u(intshift+1:end)))) 0]./tpt);

delay = intshift;

function y=erb(x)
y=24.7*(4.37e-3*x+1);
