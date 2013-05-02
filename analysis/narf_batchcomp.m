
function narf_batchcomp(batchid,modellist);
    
% query celldb to find out quantitive prediction performance of the
% difference models

% get linear model performance (firn)
modelname1='env100_log2b_firn_npfnl_boost_mses5';

modelcount=length(modellist);
cellset=cell(modelcount,1);
rset=cell(modelcount,1);
cellids={};
for ii=1:modelcount,
    
    sql=['SELECT * FROM NarfResults WHERE modelname="',...
         modellist{ii},'" and batch=',num2str(batchid),...
         ' ORDER BY cellid'];
    narfdata=mysql(sql);
    cellset{ii}={narfdata.cellid};
    rset{ii}=cat(1,narfdata.r_test);
    cellids=union(cellids,cellset{ii});
end

cellcount=length(cellids);
rmtx=zeros(cellcount,modelcount);

for cc=1:cellcount,
    for ii=1:modelcount,
    
        % find match between linear model database entry and celldata entry
        cid=find(strcmp(cellids{cc},cellset{ii}));
    
        % record model performance
        if isempty(cid),
            fprintf('no %s match for cellid %s\n',modellist{ii},cellids{cc});
            rmtx(cc,ii)=nan;
        else
            rmtx(cc,ii)=rset{ii}(cid);
        end
    end
end

pairs=nchoosek(1:modelcount,2);
for pidx=1:size(pairs,1),
    figure
    ff=find(~isnan(rmtx(:,pairs(pidx,1))) & ...
            ~isnan(rmtx(:,pairs(pidx,2))));
    
    plot(rmtx(ff,pairs(pidx,1)),rmtx(ff,pairs(pidx,2)),'.');
    hold on;
    plot([-0.2 0.8],[-0.2 0.8],'k--');
    hold off
    axis tight square
    xlabel(sprintf('%s (%.3f)',modellist{pairs(pidx,1)},...
                   mean(rmtx(ff,pairs(pidx,1)))),'Interpreter','none');
    ylabel(sprintf('%s (%.3f)',modellist{pairs(pidx,2)},...
                   mean(rmtx(ff,pairs(pidx,2)))),'Interpreter','none');
end

return

keyboard

    % find match between depression model database entry and celldata entry
    cid=find(strcmp(celldata(ii).cellid,depncells));
    
    % record model performance
    if isempty(cid),
        fprintf('no depn match for cellid %s\n',celldata(ii).cellid);
        celldata(ii).r_dep=0;
    else
        celldata(ii).r_dep=depn_data(cid).r_test;
    end
    
    fprintf('%-13s %3d %4d %6.2f %12s %6.3f %6.3f \n',...
        celldata(ii).cellid,celldata(ii).centermatch,...
        celldata(ii).BF,...
        celldata(ii).SNR,...
        celldata(ii).CenterStim,...
        celldata(ii).r_lin,...
        celldata(ii).r_dep);

cmatch=cat(1,celldata.centermatch);
grtr=find(cmatch>=60);
lesr=find(cmatch<60);
r_lin=cat(1,celldata.r_lin);
r_dep=cat(1,celldata.r_dep);

figure
plot(r_lin(grtr),r_dep(grtr),'r.');
hold on
plot(r_lin(lesr),r_dep(lesr),'g.')
minval=-0.2;
maxval=0.9;
plot([minval maxval],[minval maxval],'k--');
hold off
axis tight
xlabel(modelname1,'Interpreter','none');
ylabel(modelname2,'Interpreter','none');

r_lin_gmean=mean(r_lin(grtr));
r_dep_gmean=mean(r_dep(grtr));
r_lin_lmean=mean(r_lin(lesr));
r_dep_lmean=mean(r_dep(lesr));

figure
bar(r_lin_gmean,r_dep_gmean);
hold on
hold off
xlabel(modelname1,'Interpreter','none');
ylabel(modelname2,'Interpreter','none');

return

if 0,
    % make a template for celldata script
    for ii=1:length(firncells),
        fprintf('celldata(%d).cellid=''%s'';\n',ii,firncells{ii});
        fprintf('celldata(%d).centermatch=\n',ii);
        fprintf('celldata(%d).BF=\n',ii);
        fprintf('celldata(%d).STRF=\n',ii);
        fprintf('celldata(%d).CenterStim=\n',ii);
    end
end

 