function r=nlkern(beta,stim);

n=size(stim,2);
r=stim*beta(1:n);

r=r-beta(n+1);
r(r<0)=0;
%r=r.*beta(n+2);
