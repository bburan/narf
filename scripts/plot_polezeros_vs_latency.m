function plot_polezeros_vs_latency(batch, cellids, modelnames)
global XXX STACK;

    function w = pack_polezeros_and_latency(unused_stack)        
        px = [];
        py = [];
        zx = [];
        zy = [];
        
        fprintf('%s\n', XXX{1}.cellid);
        
        for ii = 1:length(STACK)
            mm = STACK{ii};
            nsplits = length(mm);
            if nsplits > 1
                error('I cannot support split/unify things yet.');
            end
            m = mm{1};
            
            if strcmp(m.name, 'pole_zeros')
                px = cat(1, px, repmat(abs(m.delays), length(m.poles), 1));
                py = cat(1, py, abs(m.poles(:)));
                zx = cat(1, zx, repmat(abs(m.delays), length(m.zeros), 1));
                zy = cat(1, zy, abs(m.zeros(:)));
            end
        end
        
        w = [];
        w.zx = zx;
        w.zy = zy;
        w.px = px;
        w.py = py;
    end

stack_extractor = @pack_polezeros_and_latency;

[params, ~, ~] = load_model_batch(batch, cellids, modelnames, ...       
                                  stack_extractor);

fig = figure('Name', 'Poles,Zeros,Delays', 'NumberTitle', 'off', 'Position', [20 50 900 900]);

% Build a big vector of coords for plotting on xy plane
allzx=[]; % Zero x coord (latency)
allzy=[]; % Zero y coord
allpx=[]; % Pole x coord (latency)
allpy=[]; % Pole y coord

for ii=1:length(params)
    w = params{ii};
    allpx = cat(1, allpx, w.px);
    allpy = cat(1, allpy, w.py);
    allzx = cat(1, allzx, w.zx);
    allzy = cat(1, allzy, w.zy);    
end

hold on;
plot(allzx, allzy, 'r.');
plot(allpx, allpy, 'b.');
set(gca,'XScale','log');
set(gca,'YScale','log');
axis tight;
xlabel('Latency');
ylabel('Pole or Zero');
hold off;

figure('Name', 'Poles/Zeros Hist', 'NumberTitle', 'off', 'Position', [20 50 900 900]);
data = log10(allpy(:));
data(:,2) = NaN;
data(1:length(allzy),2) = log10(allzy);
[y,b] = hist(data, 100);
bar(b,y, 'grouped', 'BarWidth', 1);
title('Pole (blue) and Zero (red) Histogram');

xlabel('Log10(pole) or log10(zero)');
ylabel('# of poles/zeros');
end