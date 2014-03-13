function plot_heatmap_old()
%     uicontrol('Parent', bottom_panel, 'Style', 'pushbutton', 'Units', 'pixels',...
%           'HorizontalAlignment', 'left', 'String', 'Heat Map', ...
%           'Position', [300+ButtonWidth*5 bh-ts ButtonWidth ts-pad], ...
%           'Callback', @heatmap_plot_callback);
% 
%     function heatmap_plot_callback(~,~,~)
%         data = compute_data_matrix({sel_metric});
%         if isempty(data)
%             return;
%         end
%         figure('Name', 'Heat Map Comparison', 'NumberTitle', 'off', ...
%                'Position', [10 10 900 900]);
%         
%         data = log(1./data);
%         hist(data(:,1), 50);
%         
%     end
% 
%     function heatmap_plot_callback_old(~,~,~)
%         data = compute_data_matrix({sel_metric});
%         if isempty(data)
%             return;
%         end
%         figure('Name', 'Heat Map Comparison', 'NumberTitle', 'off', ...
%                'Position', [10 10 900 900]);      
%         ns = 1:size(data,2);
%         D = sortrows([nanmean(data)' ns' data'], [-1]);
%         idxs = D(:,2);       
%         heatmap(D(:, 3:end), sel_cellids, sel_models(idxs),  '', 'TickAngle', 90, 'ShowAllTicks', true);         
%     end


%         n_pieces = 10;
%         
%         % Sort the data
%         Dmeans = nanmean(data);
%         Dcount = sum(~isnan(data));
%         D = [nanmean(data)' data'];
%         [sD, idxs] = sortrows(D, -1);
%         data = sD(:, 2:end)';
%         data = abs(data); % Neg correlations are still FINE, I say!
%         figure('Name', 'Bar Plot', 'NumberTitle', 'off', ...
%                'Position', [10 10 900 300]);
%         hold on;
%         len = length(sel_models);
%         
%         if (mod(size(data,1), n_pieces) == 0)
%             np = max(1, floor(size(data, 1) / n_pieces));        
%         else
%             np = max(1, floor(size(data, 1) / (n_pieces-1)));
%         end
%         deciles = conv_fn(sort(data, 1, 'descend'), 1, @(x) max(x(:)), np, 0);
%         decadiffs = diff(flipud(deciles));
%         decas = [deciles(end,:); decadiffs];
%         bar(1:len, decas', 'stacked');          
%         bar(1:len, nanmean(data), 0.3, 'r'); 
%         %errorbar(1:len, nanmean(data), var(data)); 
%         plot(1:len, data, 'k.');        
%         %errorbar(1:len, nanmean(data), nanvar(data), max(data) - nanmean(data), 'xk');
%         hold off;
%         set(gca,'XTick', 1:len);
%         axis([0 len+1 -0.1 max(data(:))*1.1]);
%         thelabels = {};
%         for ii = 1:len
%             thelabels{ii} = [ '(' num2str(Dcount(ii)) ')[' sprintf('%5.3f', Dmeans(ii)) ']' sel_models{ii} ];  
%         end
%         set(gca,'XTickLabel', thelabels(idxs));
%         set(gca,'CameraUpVector',[-1,0,0]);
%         title(sprintf('Model Performance and Deciles (Batch %d, %s)',...
%                       sel_batch,datestr(now)));
%         keyboard;
%         xt = get(gca, 'XTick');
%         set(gca,'CameraUpVector',[-1,0,0]);
%         
%         hline2 = line(t+.06,sin(t),...
%              'LineWidth',4,'Color',[.8 .8 .8]);