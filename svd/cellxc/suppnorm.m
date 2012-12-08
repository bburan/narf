% function r=suppnorm(nlparms,rsep)
%
% r= | (c0 + e) / (c1 + i) - c2 |+
%
% where nlparms= [c0 c1] -- c0 and c1 are scalars
% and   rsep= [ e i ]    -- e and i are T x 1 vectors
%
% created SVD 11/19/04
%
function r=suppnorm(nlparms,rsep)

c0=nlparms(1);
c1=nlparms(2);
c2=nlparms(3);
ee=rsep(:,1);
ii=rsep(:,2);

% protect against ii getting too small

minhib=0.25;
ii(find(c1+ii<minhib))=minhib-c1;

r= (c0 .* (ee-ii)) ./ (c1 + ii) - c2;
r(r<0)=0;


