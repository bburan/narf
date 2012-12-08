function [ndx,Pi,cost]= grPartition(C,k,nrep);
%
% function [ndx,Pi,cost]= grPartition(C,k,nrep);
%
% Partitions the n-node undirected graph G defined by the matrix C
% 
% Inputs:
% C - n by n edge-weights matrix. In particular, c(i,j)=c(j,i) is equal 
%     to the cost associated with cuting the edge between nodes i and j.
%     This matrix should be symmetric and doubly stochastic. If this
%     is not the case, this matrix will be normalized to
%     satisfy these properties (with a warning).
% k - desired number of partitions
% nrep - number of repetion for the clustering algorithm 
%       (optional input, defaults to 1)
% 
% Outputs:
% ndx  - n-vector with the cluster index for every node 
%       (indices from 1 to k)
% Pi   - Projection matrix [see Technical report
% cost - cost of the partition (sum of broken edges)
%
% By Joao Pedro Hespanha, Copyright 2004
%
% Example:
%
% X=rand(200,2);               % place random points in the plane
% C=pdist(X,'euclidean');      % compute distance between points
% C=exp(-.1*squareform(C));    % edge cost is a negative exponential of distance
%
% k=6;                         % # of partitions
% [ndx,Pi,cost]= grPartition(C,k,30);
%
% colors=HSV(k);               % plots points with appropriate colors
% colormap(colors)
% cla
% line(X(:,1),X(:,2),'MarkerEdgeColor',[0,0,0],'linestyle','none','marker','.');
% for i=1:k
%   line(X(find(ndx==i),1),X(find(ndx==i),2),...
%       'MarkerEdgeColor',colors(i,:),'linestyle','none','marker','.');
% end
% title(sprintf('Cost %g',cost))
% colorbar

if nargin<3
  nrep=1;
end

[n,m]=size(C);
if n~=m
  error('grPartition: Cost matrix is not square'); 
end  

if ~issparse(C)
  C=sparse(C);  
end

% Test for symmetry
if any(any(C~=C'))
  warning('grPartition: Cost matrix was not symmetric. Working with symmetric compenent')
  % Make C symmetric  
  C=(C+C')/2;
end

% Test for double stochasticity
if any(sum(C,1)~=1)
   
  warning('grPartition: Cost matrix was not doubly stochastic. Working with normalized costs')
  % Make C double stochastic
  C=C/(max(sum(C)));
  C=C+sparse(1:n,1:n,1-sum(C));
end

if any(any(C<0))
  error('grPartition: Edge costs cannot be negative')
end  

% Spectral partition
options.issym=1;               % matrix is symmetric
options.isreal=1;              % matrix is real
options.tol=1e-6;              % decrease tolerance 
options.maxit=500;             % increase maximum number of iterations
options.disp=0;
[U,D]=eigs(C,k,'la',options);  % only compute 'k' largest eigenvalues/vectors

% Simple case
%[Q,R,E]=qr(U',0);
%Z=Q*(R(:,1:k)')^(-1);
%[dummy,ndx]=max(U*Z,[],2);
%Pi=sparse(1:length(ndx),ndx,1);

% Clustering
[ndx,zs]=kmeans(U,k,'Distance','cosine','Start','sample','Replicates',nrep);

if nargout>1
  Pi=sparse(1:length(ndx),ndx,1);
end  

if nargout>2
  cost=full(sum(sum(C))-trace(Pi'*C*Pi));
end

return

% Example:
%
X=rand(200,2);               % place random points in the plane
C=pdist(X,'euclidean');      % compute distance between points
C=exp(-.1*squareform(C));    % edge cost is a negative exponential of distance

k=6;                         % # of partitions
[ndx,Pi,cost]= grPartition(C,k,30);

colors=HSV(k);               % plots points with appropriate colors
colormap(colors)
cla
line(X(:,1),X(:,2),'MarkerEdgeColor',[0,0,0],'linestyle','none','marker','.');
for i=1:k
  line(X(find(ndx==i),1),X(find(ndx==i),2),...
      'MarkerEdgeColor',colors(i,:),'linestyle','none','marker','.');
end
title(sprintf('Cost %g',cost))
colorbar
