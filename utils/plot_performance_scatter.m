function plot_performance_scatter (batch, cellid, holdtokens, modelname, ordinate)

if ~ismember(ordinate,{'r_test','r_fit','sparsity'})
    disp('invalid ordinate specified, using r_test instead');
    ordinate='r_test';
end

models = db_get_models(batch, cellid, holdtokens);

modelnames=cell(size(models));
for ii=1:length(models),
    modelnames{ii}=char(models(ii).modelname);
end

umodelnames=unique(modelnames);
umodelnames=setdiff(umodelnames,{modelname});

modelcount=length(umodelnames);
preddata=cell(modelcount,1);
keepidx=zeros(modelcount,1);
for ii=1:modelcount,
    sql=['SELECT NarfResults.r_test as r1,n1.r_test as r2,',...
         ' n1.cellid,n1.modelname',...
         ' FROM NarfResults INNER JOIN NarfResults n1',...
         ' ON NarfResults.cellid=n1.cellid',...
         ' AND NarfResults.modelname="',modelname,'"',...
         ' AND n1.modelname="',umodelnames{ii},'"',...
         ' AND NarfResults.batch=',num2str(batch),...
         ' AND n1.batch=',num2str(batch)];
    resdata=mysql(sql);
    preddata{ii}=[cat(1,resdata.r1) cat(1,resdata.r2)];
    if length(resdata)>1,
        keepidx(ii)=1;
    end
end

keepidx=find(keepidx);
umodelnames=umodelnames(keepidx);
preddata=preddata(keepidx);
modelcount=length(keepidx);


colcount=ceil(sqrt(modelcount));
rowcount=ceil(modelcount./colcount);
figure;
name1=modelname;
ff=find(name1=='_');
ff=min(ff(ff>20));
if ~isempty(ff),
    name1(ff)=char(10);
end
minax=-0.3
maxax=1;
for pidx=1:modelcount,
    vals=preddata{pidx};
    
    if ~isempty(vals),
        subplot(rowcount,colcount,pidx);
        plot([minax maxax],[minax maxax],'k--');
        hold on
        plot(vals(:,1),vals(:,2),'k.');
        hold off
        name2=umodelnames{pidx};
        ff=find(name2=='_');
        ff=min(ff(ff>20));
        if ~isempty(ff),
            name2(ff(1))=char(10);
        end
        
        title(sprintf('%s',name2),'Interpreter', 'none');
        aa=axis;
        text(minax+0.02,maxax-0.01,sprintf('%.3f\n%.3f', mean(vals(:,2))),...
             'VerticalAlignment','top');
        text(maxax,minax+0.02,sprintf('%.3f\n%.3f', mean(vals(:,1))),...
             'HorizontalAlignment','right','VerticalAlignment','bottom');
        axis tight square
    end
end

set(gcf,'Name',[name1 ' versus']);

return


models=cell(length(freetokens),1);
modelnames=cell(length(freetokens),1);
cellids=cell(length(freetokens),1);

for jj=1:length(models{token_idx}),
    modelnames{token_idx}{jj} = char(models{token_idx}(jj).modelname);
    cellids{token_idx}{jj} = char(models{token_idx}(jj).cellid);
end

pmtx=zeros(length(freetokens));
for p1=1:length(freetokens),
    for p2=(p1+1):length(freetokens),
        teststr=strrep(modelnames{p1}{1},freetokens{p1},freetokens{p2});
        pmtx(p1,p2)=sum(strcmpi(teststr,modelnames{p2}));
    end
end
[p1set,p2set]=find(pmtx);



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


