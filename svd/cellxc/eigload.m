% function [eigmatrix,eigval]=eigload(cellid)
%
% created SVD 10/08/03 -- ripped off of eigcomp and movgetcheck
%
function [eigmatrix,eigval]=eigload(cellid)

eigpath='/auto/k2/share/data/eigdata/';
eigfile=sprintf('%s%s.eig.mat',eigpath,cellid);

if ~exist(eigfile,'file'),
   [eigmatrix,eigval]=eigcomp(cellid);
else
   load(eigfile);
end
