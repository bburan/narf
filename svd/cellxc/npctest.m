cheatfile='/home2/svd/npc/02-fieldL-cheatpreds.mat';
cheatfile2='/home2/svd/npc/02-fieldL-cheat50preds.mat';
cheatfile3='/home2/svd/npc/02-fieldL-cheat25preds.mat';


cd /home2/data/npc/data/02_fieldL_song/val_data

cellinfo;

prediction=[];
for ii=1:length(celldata),
   pp=load(celldata(ii).datafile);
   prediction(ii).cellid=celldata(ii).cellid;
   prediction(ii).response=pp.resp;
end

save(cheatfile,'prediction');

for ii=1:length(celldata),
   % make each "prediction" 50% noise
   prediction(ii).response=(prediction(ii).response+...
       shuffle(prediction(ii).response))./2;
end

save(cheatfile2,'prediction');

for ii=1:length(celldata),
   % make each "prediction" 75% noise
   prediction(ii).response=(prediction(ii).response+...
       shuffle(prediction(ii).response))./2;
end

save(cheatfile3,'prediction');
