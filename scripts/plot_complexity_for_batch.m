function plot_complexity_for_batch(batch, cellids, modelnames)
% Makes a large complexity plot showing all valid models for the current
% batch

% Get the number of cells in the batch
cells = request_celldb_batch(batch);
N_cells = length(cells);

% Get a list of all the modelnames
sql = ['SELECT DISTINCT(modelname) FROM NarfResults WHERE batch=' num2str(batch) ''];
results = mysql(sql);
N = length(results);
models = {};
for ii = 1:N
    models{ii} = char(results(ii).modelname);
end

% Load the modelnames and add them to the data IFF they are complete
% in both r_ceil and n_parms
data = [];
modelnames = {};
jj = 1;
for ii = 1:N
    sql = ['SELECT * FROM NarfResults WHERE batch=' num2str(batch)];
    sql = [sql ' AND modelname="' models{ii} '"'];
    results = mysql(sql);
    if length(results) ~= N_cells || any(isempty([results.r_ceiling])) || any(isempty([results.n_parms]))        
        continue;
    else
        data(:, jj, 1) = [results.n_parms];
        data(:, jj, 2) = [results.r_test];        
        modelnames{jj} = models{ii};
        jj = jj + 1;
    end
end

plot_complexity(data, modelnames, 'r_ceiling');
