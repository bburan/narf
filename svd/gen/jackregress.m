function [b,berr]=jackregress(y,x,N);
    
    if ~exist('N','var'),
        N=min(20,length(y));
    end
    
    bb=length(y)./N;
    
    bj=zeros(size(x,2),N);
    for ii=1:N,
        jidx=[1:round((ii-1)*bb) round(ii*bb+1):length(y)];
        bj(:,ii)=regress(y(jidx),x(jidx,:));
    end
    b=mean(bj,2);
    berr=std(bj,0,2) * sqrt(N-1);
    