function plot_performance_scatter (batch, cellid, holdtokens, freetokens, ordinate)

if ~ismember(ordinate,{'r_test','r_fit','sparsity'})
    warning('invalid ordinate specified, using r_test instead');
    ordinate='r_test';
end

    
n_pts = 200;
token_count = length(freetokens);
x = zeros(n_pts, token_count);
y = zeros(n_pts, token_count);
mdlcount = 0;

models=cell(length(freetokens),1);
modelnames=cell(length(freetokens),1);
cellids=cell(length(freetokens),1);
for token_idx = 1:length(freetokens)
    
    models{token_idx} = db_get_models(batch, cellid, cat(2, holdtokens,freetokens{token_idx}));
    for jj=1:length(models{token_idx}),
        modelnames{token_idx}{jj} = char(models{token_idx}(jj).modelname);
        cellids{token_idx}{jj} = char(models{token_idx}(jj).cellid);
   end    
end

pmtx=zeros(length(freetokens));
for p1=1:length(freetokens),
    for p2=(p1+1):length(freetokens),
        teststr=strrep(modelnames{p1}{1},freetokens{p1},freetokens{p2});
        pmtx(p1,p2)=sum(strcmpi(teststr,modelnames{p2}));
    end
end
[p1set,p2set]=find(pmtx);

if isempty(p1set),
    disp('no valid token pairs');
    return
end

colcount=ceil(sqrt(length(p1set)));
rowcount=ceil(length(p1set)./colcount);
figure;
for pidx=1:length(p1set),
    p1=p1set(pidx);
    p2=p2set(pidx);
    vals=[];
    for jj=1:length(modelnames{p1}),
        teststr=strrep(modelnames{p1}{jj},freetokens{p1},freetokens{p2});
        
        ff=find(strcmp(teststr,modelnames{p2}) & ...
                strcmp(cellids{p1}{jj},cellids{p2}));
        if ~isempty(ff),
            vals=cat(1,vals,[models{p1}(jj).(ordinate) ...
                             models{p2}(ff).(ordinate)]);
            %fprintf('%s %.3f - %.2f\n',cellids{p1}{jj},vals(end,:));
        end
    end
    if ~isempty(vals),
        subplot(rowcount,colcount,pidx);
        plot([-0.2 1.0],[-0.2 1.0],'k--');
        hold on
        plot(vals(:,1),vals(:,2),'k.');
        hold off
        xlabel(sprintf('%s (mean %.3f)',freetokens{p1},mean(vals(:,1))), 'Interpreter', 'none');
        ylabel(sprintf('%s (mean %.3f)',freetokens{p2},mean(vals(:,2))), 'Interpreter', 'none');
        tstr=[ordinate ' ('];
        for ii=1:length(holdtokens),
            tstr=[tstr holdtokens{ii} ,','];
        end
        tstr(end)=')';
        
        title(tstr, 'Interpreter', 'none');
        axis tight square
    end
end


