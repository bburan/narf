function plot_performance_bar (batch, cellid, holdtokens, freetokens, ordinate)

if ~ismember(ordinate,{'r_test','r_fit','sparsity'})
    fprintf(['\nplot_performance_bar.m: invalid ordinate specified, ' ...
             'using r_test\n']);
    ordinate='r_test';
end


n_pts = 200;
token_count = length(freetokens);
x = zeros(n_pts, token_count);
y = zeros(n_pts, token_count);
mdlcount = 0;

models=db_get_models(batch,'modelnames',holdtokens);
modelcount=length(models);
modelnames=cell(modelcount,1);
celldata=mysql(['SELECT DISTINCT cellid FROM NarfResults',...
                ' WHERE batch=',num2str(batch)]);
cellids={celldata.cellid};
ordinates=ones(length(cellids),modelcount).*nan;
for jj=1:length(models),
    modelnames{jj} = ...
        char(models(jj).modelname);
    celldata=mysql(['SELECT * FROM NarfResults WHERE batch=',num2str(batch),...
                    ' AND modelname="',modelnames{jj},'"']);
    for ii=1:length(celldata),
        cc=find(strcmp(celldata(ii).cellid,cellids));
        ordinates(cc,jj)=celldata(ii).(ordinate);
    end
end
[~,modelorder]=sort_nat(modelnames);

modelnames=modelnames(flipud(modelorder));
ordinates=ordinates(:,flipud(modelorder));

figure;

subplot('position',[0.4 0.55 0.57 0.42]);
imagesc(ordinates');
axis xy
set(gca,'YTick',1:modelcount,'YTickLabel',modelnames);
xlabel('cell');
colorbar

subplot('position',[0.4 0.05 0.57 0.42]);
ff=find(sum(isnan(ordinates),2)==0);
bb=mean(ordinates(ff,:));
ee=std(ordinates(ff,:)./sqrt(length(ff)));
for jj=1:modelcount,
    plot([bb(jj)-ee(jj) bb(jj)+ee(jj)],[jj jj]);
    hold on
end
barh(bb);
hold off
aa=axis;
axis([aa(1:2) 0 modelcount+1]);

set(gca,'YTick',1:modelcount,'YTickLabel',modelnames);
xlabel(['Mean ',ordinate],'Interpreter','none');


