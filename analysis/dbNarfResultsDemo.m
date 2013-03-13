

% clear any existing test cell entries (be careful)
cellid='test001a-a1';
sql=['DELETE FROM NarfResults WHERE cellid="',cellid,'"'];
mysql(sql);

% insert sample data into NarfResults:
batch=240;
modelname='testmodel_A_B_C';
r_test=0.5;
r_fit=1;
modelpath='/save/path';
modelfile='testmodel_A_B_C.mat';
sqlinsert('NarfResults',...
          'cellid',cellid,...
          'batch',batch,...
          'r_test',r_test,...
          'r_fit',r_fit,...
          'modelname',modelname,...
          'modelpath',modelpath,...
          'modelfile',modelfile);

modelname='testmodel_D_E_F';
modelfile='testmodel_D_E_F.mat';
r_test=0.9;
sqlinsert('NarfResults',...
          'cellid',cellid,...
          'batch',batch,...
          'r_test',r_test,...
          'r_fit',r_fit,...
          'modelname',modelname,...
          'modelpath',modelpath,...
          'modelfile',modelfile);

% find all matching entries
sql=['SELECT * FROM NarfResults WHERE modelname like "%_A_%"'];
resultsdata=mysql(sql);

% "char" command required to translate "text" fields back to
% strings.  Not required for varchar fields.
modelname=char(resultsdata(1).modelname)


% change existing entry
sql=['UPDATE NarfResults SET r_test=0.7 WHERE id=',...
     num2str(resultsdata(1).id)];
mysql(sql);


% sort all data by r_test:
% find all matching entries
sql=['SELECT * FROM NarfResults WHERE cellid="',cellid,'"',...
     ' ORDER BY r_test'];
resultsdata=mysql(sql);

for ii=1:length(resultsdata),
    fprintf('model: %s  r_test: %.1f\n',...
            resultsdata(ii).modelname,resultsdata(ii).r_test);
end
