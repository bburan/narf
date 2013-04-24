function models = db_get_models(batch, cellid, matching_tokens, sort_field) 
% models = db_get_models(batch, cellid, matching_tokens, sort_field) 
%
% Returns all models which have ALL of the matching tokens in their
% modelnames, sorted by sortfield. If you want to search for all cellids,
% please leave it empty ('')

% ARGUMENTS:
%    batch            An integer number. 
%    cellid           Using an empty string here matches all cellids.
%                     Default: ''
%    matching_tokens  A cell array of tokens found in module_groups()
%                     Example: {'firn', 'npnl'}
%                     Default: {}  (match everything)
%    sortfield        Which field to sort by.
%                     The default is r_test.
%
% RETURNS: 
%    models           A structarray returned directly from the mysql query.


constraints = {};

sql = 'SELECT * FROM NarfResults';

if exist('batch', 'var')
    constraints{end+1} = ['batch=' num2str(batch)];
end

if exist('cellid', 'var') && ~isempty(cellid)
    constraints{end+1} = ['cellid="' cellid '"'];
end

if exist('matching_tokens', 'var') && ~isempty(matching_tokens)
    for ii = 1 :length(matching_tokens)
        constraints{end+1} = ['modelname REGEXP "(^|[^[:alnum:]])' ...
                                matching_tokens{ii}  ...
                                '([^[:alnum:]]|$)"'];
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

models = mysql(sql);
