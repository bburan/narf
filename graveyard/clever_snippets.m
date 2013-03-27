% % ----------------------
% For plotting a test/train scatter plot
% % % Also, show the test/train curves
% fh_tt = figure;
% plot(1:length(filenames), score(:,1), ...
%      1:length(filenames), score(:,3));
% xticks(1:length(filenames));
% 
% % label the x axis with filenames at a 90 degree angle
% lab = char(sorted_filenames);
% hx = get(gca,'XLabel');
% set(hx, 'Units', 'data');
% pos = get(hx, 'Position');
% 
% for ii = 1:size(lab,1)
%     t(ii) = text(ii, pos(2), lab(ii,:), 'Interpreter', 'none');
% end
% 
% P=get(gca,'Position');
% set(t, 'Rotation', 90, 'HorizontalAlignment', 'right') 
% set(gca, 'Position', [0.1300    0.4100    0.7750    0.5150])
% 
% legend('Test Score', 'Training Score', 'Location', 'NorthWest');
% 
% title(filenames{1}, 'Interpreter', 'none');

% fh1 = fig;
% fh2 = fh_tt;

        % ---------------------------
% FIR_FILTER.m
        % TODO: PLACEHOLDER FOR ARBITRARY DIMENSION FILTERING         
%         % Build up the input matrix
%         M = x.dat.(sf).(mdl.input);
%         M_selcols = x.dat.(sf).([mdl.input '_selcols']);
%         
%         chan_idx = 1;
%         
%         % The fir filter acts on every channel
%         for chan_idx = 1:mdl.num_dims
%             cols = M_selcols('chan', chan_idx);
%             X = M(:, cols);
%             [T, N] = size(X);
%             tmp = zeros(T, N); 
%             % Apply the filter to every other column, one at a time
%             for d = 1:N,
%                 tmp(:, d) = filter(mdl.coefs(:, chan_idx), [1], X(:, d));
%             end
%             % Store the results back in the data structure.
%             x.dat.(sf).(mdl.output)(:, cols) = tmp; 
%         end
%        % ------------------------
