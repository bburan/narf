function db_insert_model()
% db_insert_model()
%
% Forcibly inserts the loaded model to the database table 'NarfResults',
% deleting any previous model that existed there. 
%
% Also creates a model plot using plot_model_summary().
%
% No arguments or return values.

global XXX META;

if ~isfield(META,'batch')
    META.batch = 0;
end

if ~isfield(META,'git_commit')
    META.git_commit = 'unknown';
end

% Don't add if matching batch, cellid, and modelname already exist
sql = ['SELECT * FROM NarfResults WHERE modelname="' META.modelname '"' ...
       ' AND batch=' num2str(META.batch) ...
       ' AND cellid="' XXX{1}.cellid '"'];
r=mysql(sql);

if length(r) == 1
    fprintf('NOTE: Deleting old NarfResults entry for %s/%d/%s\n',...
            XXX{1}.cellid,META.batch,META.modelname);
    sql = ['DELETE FROM NarfResults WHERE id=', num2str(r(1).id)];
    mysql(sql);
elseif length(r) > 1
    error('Duplicate values in DB found!');
end

r_test = XXX{end}.score_test_corr;

if isnan(r_test)
    r_test = 0;
end

% Otherwise, generate a model plot and insert the results into the DB
plotpath = plot_model_summary();
[affected, ~] = sqlinsert('NarfResults', ...
          'cellid',    XXX{1}.cellid,...
          'batch',     META.batch,...
          'r_fit',     XXX{end}.score_train_corr,...
          'r_test',    r_test,...
          'score',     XXX{end}.score, ...
          'sparsity',  XXX{end}.sparsity, ...
          'modelname', META.modelname, ...
          'modelpath', META.modelpath, ...
          'modelfile', META.modelfile, ... 
          'githash',   META.git_commit, ...
          'figurefile', plotpath);
      
if affected ~= 1
    error('The number of affected sql entries was not exactly 1!');
end
      
% fprintf('affected %d', affected);

end
