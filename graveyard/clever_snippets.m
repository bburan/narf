
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
