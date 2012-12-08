%showSCAKernels.m
function [kernelArray] = showSCAKernels2(kernel, S, numComponents)

%input: structural array containing kernels for each component @ each latency
%      - S, signal matrix 
%      - # components to visualize (for now we're trying 20)
%output: a series of weighted kernels at all latencies (one latency string for each component)

%number of latencies over which kernel is computed.
numLat = size(kernel,2);
tileSize = floor(size(S,1)^.5);

p = sum(kernel.^2,2);

[Y,I] = sort(p);
I = flipdim(I,1);

goodidx = I(1:numComponents,1);
figure; plotnum=1; colormap(gray)
for ii = 1:length(goodidx)
    for tau = 1:numLat 
        hw(ii,tau,:) = (S(:,goodidx(ii)).*kernel(goodidx(ii),tau))';
        
    end
end
kernelArray = hw;

ub = max(max(max(hw)));
for ii = 1:length(goodidx)
    for tau = 1:numLat
        subplot(numComponents,numLat,plotnum), imagesc((reshape(hw(ii,tau,:),tileSize,tileSize)),[-ub,ub]);
        plotnum=plotnum+1;
    end
    
end
