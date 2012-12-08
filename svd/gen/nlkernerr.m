function err=nlkernerr(beta,stim,resp);

n=size(stim,2);
r=stim*beta(1:n);

r=r-beta(n+1);
r(r<0)=0;
%r=r.*beta(n+2);

% test smoothness



err=mean((r-resp).^2);


