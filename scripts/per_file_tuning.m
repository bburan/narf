function per_file_tuning(batch, cellids, modelnames)

    function ret = getit(x, idx)
        if isempty(x)
            ret = nan;
        else
            ret = x(idx);
        end
    end
    
    function ret = getfitval(x)
        if ~isfield(x, 'fit_time') || ~isfield(x, 'perf_val_corr')
            ret = [nan nan];
        else
            ret = [x.fit_time x.perf_val_corr];
        end
    end

stack_extractor = @(x) x;
meta_extractor = @getfitval;
% stacks is modelnames X cellid
[stacks, metas, x0s] = load_model_batch(batch, cellids, modelnames, ...
                                        stack_extractor, meta_extractor);
cellcount=length(cellids);
modelcount=length(modelnames);
for cc=1:cellcount,
    figure;
    for mm=1:modelcount,
        subplot(modelcount,1,mm);
        
        [nlmod] = find_modules(stacks{mm,cc}, 'nonlinearity', false);
        nlmod=nlmod{end};
        if mm==1,
            xx=linspace(-100,200,300)';
            tyy=nlmod{1}.nlfn(nlmod{1}.phi,xx);
            dtyy=abs(diff(tyy));
            m1=min(find(dtyy>=max(dtyy)./100));
            m2=max(find(dtyy>=max(dtyy)./100));
            xx=linspace(xx(m1),xx(m2),100);
        end
        yy=zeros(length(xx),length(nlmod));
        for nn=1:length(nlmod),
            yy(:,nn)=nlmod{nn}.nlfn(nlmod{nn}.phi,xx);
        end
        plot(xx,yy);
        if length(nlmod)>1,
            legend(x0s{mm,cc}.filecodes,'Location','northwest');
        end
    end
    subplot(modelcount,1,1);
    title(cellids{cc});
end

return




times = cellfun(@(x) getit(x,1), metas);
corrs = cellfun(@(x) getit(x,2), metas);

figure('Name', 'Fit Time vs Corr', 'NumberTitle', 'off', 'Position', [20 50 900 900]);
plot(times', corrs', '.');
h = legend(modelnames);
set(h, 'Interpreter', 'none');
xlabel('Fit Time');
ylabel('Val. Set Correlation');
title('Val. Set Correlation vs Fit Time');
n = ceil(sqrt(size(times,1)));
figure;
for ii = 1:size(times, 1)
    subplot(n,n,ii);
	hist(times(ii,:), ceil(size(times,2)/5));
    title([modelnames{ii} '(' num2str(nanmean(times(ii))) ')'], 'Interpreter', 'none');
end

end