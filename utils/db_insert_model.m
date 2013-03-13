function db_insert_model()
% Inserts the current model to the database table 'NarfResults'

global STACK XXX META;

if ~isfield(META,'batch')
    META.batch = 242;
end

if ~isfield(META,'git_commit')
    META.git_commit = 'unknown';
end

% Don't add if matching batch, cellid, and modelname already exist
sql=['SELECT * FROM NarfResults WHERE modelname="' META.modelname '"'...
    ' AND batch=' num2str(META.batch) ...
    ' AND cellid="' XXX{1}.cellid '"'];
r=mysql(sql);
if length(r) == 1
    return
elseif length(r) > 1
    error('Duplicate values in DB found!');
    keyboard;
end

r_test = XXX{end}.score_test_corr;
if isnan(r_test)
    r_test = 0;
end

% Otherwise, generate a model plot and insert the results into the DB
plotpath = plot_model_summary();
[affected, id] = sqlinsert('NarfResults', ...
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

% fprintf('affected %d', affected);

end
