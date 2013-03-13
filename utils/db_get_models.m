function models = db_get_models(batch, cellid, matching_tokens, sort_field) 
% Returns all models which have ALL of the matching tokens in their
% modelnames, sorted by sortfield.
% ARGUMENTS:
%    batch            An integer number. 
%    cellid           Self explanatory. 
%    matching_tokens  A cell array of tokens found in module_groups()
%                     Example: {'firn', 'npnl'}
%                     An empty cell array matches everything (default)
%    sortfield        Which field to sort by (default is r_test)

constraints = {};

sql = ['SELECT * FROM NarfResults'];

if exist('batch', 'var')
    constraints{end+1} = ['batch=' num2str(batch)];
end

if exist('cellid', 'var')
    constraints{end+1} = ['cellid="' cellid '"'];
end

if exist('matching_tokens', 'var') && ~isempty(matching_tokens)
    for ii = 1 :length(matching_tokens)
        constraints{end+1} = ['modelname like "%' matching_tokens{ii} '%"'];
    end
end

if ~exist('sort_field', 'var')
    sort_field = 'r_test';
end
    
% Unwrap the constraints, if any
for ii = 1:length(constraints),
    if ii == 1
        sql = [sql ' WHERE ' constraints{ii}];
    else
        sql = [sql ' AND ' constraints{ii}];
    end
end

% Add the sort order
sql = [sql ' ORDER BY ' sort_field];

models=mysql(sql);
