function ax = plot_bar_pretty(data, modelnames, show_datapoints)
        
        if ~exist('show_datapoints', 'var')
            %show_datapoints = true;
            show_datapoints = false
        end

        barwidth = 0.8;

        % Sort data by mean to make significance algorithm work
        data = abs(data); 
        Dunsortedmean = nanmean(data);        
        Dunsortedcount = sum(~isnan(data));
        D = [nanmean(data)' data'];
        [sD, idxs] = sortrows(D, -1);
        data = sD(:, 2:end)';
               
        Dmean = nanmean(data);
        Dcount = sum(~isnan(data));
        Dstddev = sqrt(nanvar(data));
        Dstderr = Dstddev ./ sqrt(Dcount);
        len = length(modelnames);
        
        fig = figure('Name', 'Bar Plot', 'NumberTitle', 'off', ...
                     'Position', [10 10 1200 50+len*30]);  
        hold on;       
        bar(1:len, Dmean, barwidth, 'r'); 
        errorbar(1:len, Dmean, Dstderr, 'k.', 'LineWidth',2);                            
        
        ax = gca;
        set(gca,'XTick', 1:len);
        if show_datapoints
            axis([0 len+0.75 -0.1 max(data(:))*1.05]);            
            % Draw the jittered data points
            for ii = 1:len
                a = data(:, ii);
                jitter = -0.15*barwidth+0.3*barwidth*rand(size(a));
                plot(ii*ones(size(a))+jitter, a, 'k.'); % Plot jittered values
            end
        else            
            axis([0 len+0.75 -0.1 max(Dmean+Dstderr)*1.05]);
        end
        thelabels = {};
        for ii = 1:len
            thelabels{ii} = [ '(' num2str(Dunsortedcount(ii)) ')[' sprintf('%5.3f', Dunsortedmean(ii)) ']' modelnames{ii} ];  
        end
        set(gca,'XTickLabel', thelabels(idxs));
        set(gca,'CameraUpVector',[-1,0,0]);
            
        % Draw the significance lines       
        sep = 0.1 / (len-1);      
        for ii=1:len-1            
            a = data(:, ii);
            z = 0;
            for jj = ii+1:len                
                b = data(:, jj);
                
                p = randttest(a-b, zeros(size(a)), 1000, 0); % Paired T test
                %[~,p2] = ttest(a,b); % 
                
                if p < 0.001 && z < 3
                    ls = '-';
                    z = 3;
                    r = -5/6*sep - sep*ii + sep;
                    l = -6/6*sep - sep*ii + sep;
                elseif p < 0.01 && z < 2
                    ls = '--';
                    z = 2;
                    r = -3/6*sep - sep*ii + sep;
                    l = -4/6*sep - sep*ii + sep;
                elseif p < 0.05 && z < 1
                    ls = ':';
                    z = 1;
                    r = -1/6*sep - sep*ii + sep;
                    l = -2/6*sep - sep*ii + sep;
                else 
                    continue;
                end
                
                line([ii ii jj jj], [r l l r], 'LineWidth', 2, 'Color', [0 0 0], 'LineStyle', ls);                
            end
        end 
    hold off;
end