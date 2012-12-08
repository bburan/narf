
close all

batch=99;
outpath=sprintf('/auto/k5/david/data/batch%d/output/',batch);

rundata=dbgetrundata(batch);

for ii=1:length(rundata);
   
   % load and process cell-specific results
   kvatune(rundata(ii).cellid,batch,'curves');
   
   drawnow;
   print('-f1','-depsc',[outpath,rundata(ii).cellid,'.beststim.eps']);
   print('-f2','-depsc',[outpath,rundata(ii).cellid,'.curves.eps']);
   print('-f3','-depsc',[outpath,rundata(ii).cellid,'.altbest.eps']);
   print('-f4','-depsc',[outpath,rundata(ii).cellid,'.altbestpfft.eps']);
   
end


batch=82;
outpath=sprintf('/auto/k5/david/data/batch%d/output/',batch);

rundata=dbgetrundata(batch);

for ii=1:length(rundata);
   
   % load and process cell-specific results
   kvatune(rundata(ii).cellid,batch,'curves');
   
   drawnow;
   print('-f1','-depsc',[outpath,rundata(ii).cellid,'.beststim.eps']);
   print('-f2','-depsc',[outpath,rundata(ii).cellid,'.curves.eps']);
   print('-f3','-depsc',[outpath,rundata(ii).cellid,'.altbest.eps']);
   print('-f4','-depsc',[outpath,rundata(ii).cellid,'.altbestpfft.eps']);
   
end
