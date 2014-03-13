function ax = plot_bar_pretty(data, modelnames)
        
        barwidth = 0.8;

        % Sort data by mean to make significance algorithm work
        data = abs(data); 
        Dunsortedmean = nanmean(data);
        D = [nanmean(data)' data'];
        [sD, idxs] = sortrows(D, -1);
        data = sD(:, 2:end)';
               
        Dmean = nanmean(data);
        Dcount = sum(~isnan(data));
        Dstddev = sqrt(nanvar(data));
        len = length(modelnames);
        
        fig = figure('Name', 'Bar Plot', 'NumberTitle', 'off', ...
                     'Position', [10 10 1200 50+len*30]);  
        hold on;       
        bar(1:len, Dmean, barwidth, 'r'); 
        errorbar(1:len, Dmean, Dstddev, 'k.', 'LineWidth',2);                            
        
        ax = gca;
        set(gca,'XTick', 1:len);
        axis([0 len+0.75 -0.1 max(data(:))*1.05]);
        thelabels = {};
        for ii = 1:len
            thelabels{ii} = [ '(' num2str(Dcount(ii)) ')[' sprintf('%5.3f', Dunsortedmean(ii)) ']' modelnames{ii} ];  
        end
        set(gca,'XTickLabel', thelabels(idxs));
        set(gca,'CameraUpVector',[-1,0,0]);

        % Draw the jittered data points
        for ii = 1:len
            a = data(:, ii);
            jitter = -0.15*barwidth+0.3*barwidth*rand(size(a));
            plot(ii*ones(size(a))+jitter, a, 'k.'); % Plot jittered values
        end
            
        % Draw the significance lines       
        sep = 0.1 / len;      
        for ii=1:len-1            
            a = data(:, ii);
            z = 0;
            for jj = ii+1:len                
                b = data(:, jj);
                
                [~, p] = ttest(a,b); % Paired T test
                
                if p < 0.001 && z < 3
                    ls = '-';
                    z = 3;
                    r = -5/6*sep - sep*ii;
                    l = -6/6*sep - sep*ii;
                elseif p < 0.01 && z < 2
                    ls = '--';
                    z = 2;
                    r = -3/6*sep - sep*ii;
                    l = -4/6*sep - sep*ii;
                elseif p < 0.05 && z < 1
                    ls = ':';
                    z = 1;
                    r = -1/6*sep - sep*ii;
                    l = -2/6*sep - sep*ii;
                else 
                    continue;
                end
                
                line([ii ii jj jj], [r l l r], 'LineWidth', 2, 'Color', [0 0 0], 'LineStyle', ls);                
            end
        end 
    hold off;
end