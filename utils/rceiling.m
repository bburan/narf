% function rnorm=rceiling(pred,resp);
%
% calculate noise ceiling on correlation coefficient for
% finite-trial spike data. to correct measured r^2 for single trial
% noise, use:
%
%  r_norm(pred,resp) = <r(pred,rsingle)> ./ <sqrt(r(rsingle,rsingle))>
%
% inputs:
% pred- vector of predictions, same length as height of resp
% resp- matrix of spike rates, one column per rep
%
% output:
% rnorm - upper bound on correlation coefficient between single
%        trial data (provided) and theoretical underlying psth
%
% created SVD 2013-10-24
%
function rnorm=rceiling(pred,resp,N)

global NARF_DEBUG

repcount=size(resp,2);
if repcount<=1,
    rnorm=0;
    return
end
if NARF_DEBUG
    keyboard
end

paircount=nchoosek(repcount,2);
pairs=zeros(paircount,2);
pp=0;
for p1=1:repcount,
    for p2=(p1+1):repcount,
        pp=pp+1;
        pairs(pp,:)=[p1 p2];
    end
end
if ~exist('N','var'),
   N=500;
end
N=min(N,size(pairs,1));
[~,sidx]=sort(rand(size(pairs,1),1));
pairs=pairs(sidx(1:N),:);

if N==1,
    % only two repeats, break up data in time to get a better
    % estimate of single-trial correlations
    N=10;
    bstep=size(pred,1)./N;
    rac=zeros(N,1);
    for nn=1:N,
        tt=round((nn-1).*bstep+1):round(nn*bstep);
        if ~isempty(tt) && std(resp(tt,1))>0 && std(resp(tt,2))>0,
            rac(nn)=xcov(resp(tt,1),resp(tt,2),0,'coeff');
        end
    end
    
else
    rac=zeros(N,1);
    for nn=1:N,
        tt=find(~isnan(resp(:,pairs(nn,1))) & ~isnan(resp(:,pairs(nn,2))));
        if ~isempty(tt) && std(resp(tt,pairs(nn,1)))>0 && std(resp(tt,pairs(nn,2)))>0,
            rac(nn)=xcov(resp(tt,pairs(nn,1)),resp(tt,pairs(nn,2)),0,'coeff');
        end
    end
end

rsingle=zeros(repcount,1);
for nn=1:repcount,
   tt=find(~isnan(resp(:,nn)));
   if ~isempty(tt) && std(resp(tt,nn))>0
       rsingle(nn)=xcov(pred,resp(tt,nn),0,'coeff');
    end
end

%figure;subplot(2,1,1);hist(rac);
%subplot(2,1,2);hist(rsingle);

rac=mean(rac);
if rac<0.05,
    rac=0.05;
end
rsingle=mean(rsingle);

rnorm=rsingle./sqrt(rac);

return

if 0,
    %test code follows:
    
    T=10000;
    N=10;
    p=randn(T,1);
    r=repmat(p,1,N);
    n=randn(size(r)).*2;
    
    r1=r+n;
    
    % r11
    
    r12=xcov(r1(:,1),r1(:,2),0,'coeff')
    %mean(r1(:,1).*r1(:,2))./sqrt(var(r1(:,1)).*var(r1(:,2)))
    
    num=var(p)+mean(n(:,1).*p)+mean(n(:,2).*p)+mean(n(:,1).*n(:,2));
    den=sqrt((var(p)+2.*mean(p.*n(:,1))+var(n(:,1))).*...
             (var(p)+2.*mean(p.*n(:,2))+var(n(:,2))));
    num/den
    
    var(p)./(var(p)+var(n(:,1)))
    
    %pr
    rpr=xcov(p,r1(:,1),0,'coeff')
    
    num=mean(p.*p)+mean(p.*n(:,1));
    den=sqrt(var(p).*(var(p)+2.*mean(p.*n(:,1))+var(n(:,1))));
    num/den
    
    var(p)./sqrt(var(p).*(var(p)+var(n(:,1))))
    
    rpr./sqrt(r12)
    
    rceiling(p,r1)
end

