function plot_wc_around_BFs(batch, cellids, modelnames)
% Plots the WC coefs centered around the BF (largest WC weight) 

    function centered_weights = extract_wc_around_BFs(~)
        global STACK;
        mods = find_modules(STACK, 'weight_channels');
        w = mods{1}{1}.weights;
        [N, M] = size(w);
        
       centered_weights = zeros(9,M);         
        
       % ----------------------------------------------------------------
       % OPTION ONE: Center the weights at the best frequency
      
        for jj = 1:M
            idx = find(abs(w(:,jj)) == max(abs(w(:,jj))), 1);
    
            % Extract +/-4 coefs
            for ii = -4:4
                if idx+ii >= 1 && idx+ii <= N
                    centered_weights(ii+5,jj) = w(idx+ii, jj);
                else
                    centered_weights(ii+5,jj) = nan;
                end
            end
        
        
            % Normalize to the BF
            %centered_weights = centered_weights ./ centered_weights(5);
        
            % Normalize by the variance
            centered_weights(:,jj) = centered_weights(:,jj) ./ nanstd(centered_weights(:,jj));
        end        
        
        % ----------------------------------------------------------------
        % OPTION TWO: Use all the weights, without recentering
        % (Showed we have more data at some freqs than others, but was hard
        % to derive many other conclusions from this)
        
        % centered_weights = w;
        % centered_weights = centered_weights ./ repmat(nanstd(centered_weights), size(w, 1), 1);
        
        % ----------------------------------------------------------------
    end

stack_extractor = @extract_wc_around_BFs;

[params, ~, ~] = load_model_batch(batch, cellids, modelnames, ...
    stack_extractor);

data = cell2mat(params(~cellfun('isempty', params)))';

pos = repmat(1:size(data,2), size(data,1), 1);
jitter = -0.15+0.3*rand(size(pos));

data2 = excise(data); % Alternative #1: Use only the "completely valid" rows
% data2 = ... % Alternative #2: Replace NANs with the mirrored values

data2(any(isinf(data2)'),:) = [];
pc = princomp(data2);

% Project data2 onto the principal components
[U, S, V] = svd(data2);
data3 = data2 * V;

figure('Name', 'WCs Around BF', 'NumberTitle', 'off', 'Position', [20 50 900 900]);

title('Centered weights around BF');

subplot(3,2,1);
plot(pos + jitter, data, 'k.');
xlabel('Relative Index Around BF');
ylabel('Weight Strength');

subplot(3,2,2);
bar(pos(1,:), mean(data2));
xlabel('Relative Index Around BF');
ylabel('Weight Strength');

subplot(3,2,3);
bar(pos(1,:), pc(:,1));
textLoc(sprintf('Eigenvalue: %f', S(1,1)), 'NorthWest');
xlabel('Relative Index Around BF');
ylabel('Principal Component #1');

subplot(3,2,4);
bar(pos(1,:), pc(:,2));
textLoc(sprintf('Eigenvalue: %f', S(2,2)), 'NorthWest');
xlabel('Relative Index Around BF');
ylabel('Principal Component #2');

subplot(3,2,5);
bar(pos(1,:), pc(:,3));
textLoc(sprintf('Eigenvalue: %f', S(3,3)), 'NorthWest');
xlabel('Relative Index Around BF');
ylabel('Principal Component #3');

subplot(3,2,6);
scatter(data3(:,1), data3(:,2),'k.');
xlabel('Principal Component #1');
ylabel('Principal Component #2');

%keyboard;

end
