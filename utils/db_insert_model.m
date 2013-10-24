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

% -----------------------------------------------------------------------
% Require values for some and give defaults for others

if ~all(isfield(META, {'batch', 'modelname', 'modelpath', 'modelfile'})) || ...
   ~all(isfield(XXX{end}, {'cellid'}))
    error('Required for DB insertion: batch, modelname, modelpath, modelfile, and cellid');
end

if ~isfield(META,'git_commit')
    META.git_commit = 'unknown';
end

if ~isfield(XXX{end},'score_train_corr')
    XXX{end}.score_train_corr = 0.0;
end

if ~isfield(XXX{end},'score_test_corr')
    XXX{end}.score_test_corr = 0.0;
end

if ~isfield(XXX{end},'score')
    XXX{end}.score = 0.0;
end

if ~isfield(XXX{end},'sparsity')
    XXX{end}.sparsity= 0.0;
end
    
% -----------------------------------------------------------------------
dbopen(1); % make sure server connection not broken during a long fit
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

plotpath = plot_model_summary();

[affected, ~] = sqlinsert('NarfResults', ...
          'cellid',    XXX{1}.cellid,...
          'batch',     META.batch,...
          'r_fit',     nan2zero(XXX{end}.score_train_corr),...
          'r_test',    nan2zero(XXX{end}.score_test_corr),...
          'r_ceiling', nan2zero(XXX{end}.score_test_ceilingcorr),...
          'r_floor',   nan2zero(XXX{end}.score_test_floorcorr),...
          'score',     nan2zero(XXX{end}.score), ...
          'sparsity',  nan2zero(XXX{end}.sparsity), ...
          'modelname', META.modelname, ...
          'modelpath', META.modelpath, ...
          'modelfile', META.modelfile, ... 
          'githash',   META.git_commit, ...
          'figurefile', plotpath);
      
if affected ~= 1
    error('The number of affected sql entries was not exactly 1!');
end

end
