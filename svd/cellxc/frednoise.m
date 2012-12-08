% function [ccr,ccrM]=frednoise(r2);
%
% to get expected asymptotic cc^2 from single trial cc:
%
% ccexp=ccp./ccr;
%
% (ccp = single trial cc^2)
% (ccexp = asymp cc^2)
%
% ccr - ratio of single trial cc^2 to infinite trial cc^2
% ccrM - ratio of single trial cc^2 to M trial cc^2
%
function [ccr,ccrM,ccre]=frednoise(r2);

r2=r2(find(~isnan(r2(:,1))),:);

M=size(r2,2);

ccp0=zeros(M,1);
ccr0=zeros(M,1);
for mm=1:M
   
   ra=mean(r2(:,1:round(M/2)),2);
   rb=mean(r2(:,round(M/2)+1:end),2);
   r2=shift(r2,[0 1]);
   
   if sum(abs(ra))>0 & sum(abs(rb))>0,
      ccr0(mm)=xcov(ra,rb,0,'coeff').^2;
   else
      ccr0(mm)=0;
   end
   
end

ccrM=(1+sqrt(1./mean(ccr0))) ./ (-M + M * sqrt(1./mean(ccr0)) + 2);
ccr=2./(2-M+M*sqrt(1./ccr0));

ccre=std(ccr,1).*sqrt(M-1);
ccr=mean(ccr);

