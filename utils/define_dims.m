function M2D = define_dims(M, dimension_names)
% Flatten an N dimensional matrix into a 2D matrix and provide a function
% which makes it easy to select only columns that match the criteria you
% provide. 
%
% Assumes that the 'time' dimension is the first dimension, so you don't
% need to explicitly list 'time' in your DIMENSION_NAMES cell array. 
% 
% The returned 2D matrix is basically a pseudoobject with these fields:
%   M2D.m          Original matrix
%   M2D.mat        Matrix compressed into 2D
%   M2D.selcols    Fn to selects columns by field name
%   M2D.length     Fn to get length of a particular field
%   M2D.foreach    Fn to get a cell array of columns idxes for every element
%
% TODO: An example of how to use this.
% 

dimension_sizes = size(M);

M2D = [];
M2D.m = M;
M2D.mat = M(:,:);                           
M2D.length = @get_length;
M2D.selcols = @select_cols_fn;
M2D.foreach = @foreach;

function len = get_length(dname)
    % Ensure that dname is a field
    if ~isempty(find(dimension_names, dname))
        error('Unknown field name: %s', dname);
    end
    dimension_sizes(find(dimension_names, dname));
end

function column_indexes = select_cols_fn(varargin)
    s = struct(varargin{:});
    
    % Check that all fields are in the dimension_names
%     for f = fieldnames(s)', f = f{1};
%         if ~isempty(find(dimension_names, f))
%             
%         else
%             error('Unknown field name: %s', f);
%         end
%     end
    
    sels = ones(prod(dimension_sizes), 1);
    for ii = 1:length(sels)
        sels(ii) = ii;
    end
    
    sels = reshape(sels, dimension_sizes);   
    mask = true(size(sels));

    for ii = 1:length(dimension_names);
        f = dimension_names{ii};
        n = dimension_sizes(ii+1);
        if isfield(s, f)
            tempmask = false(size(sels));
            if ii == 1
                tempmask(s.(f), :,:,:,:,:,:,:,:) = true;
            elseif ii == 2
                tempmask(:, s.(f),:,:,:,:,:,:,:) = true;
            elseif ii == 3
                tempmask(:,:,s.(f), :,:,:,:,:,:) = true;
            elseif ii == 4
                tempmask(:,:,:,s.(f), :,:,:,:,:) = true;
            elseif ii == 5
                tempmask(:,:,:,:,s.(f), :,:,:,:) = true;
            elseif ii == 6
                tempmask(:,:,:,:,:,s.(f), :,:,:) = true;
            elseif ii == 7
                tempmask(:,:,:,:,:,:,s.(f), :,:) = true;
            elseif ii == 8
                tempmask(:,:,:,:,:,:,:, s.(f),:) = true;
            elseif ii == 9
                tempmask(:,:,:,:,:,:,:,:, s.(f)) = true;
            else
                error('I only support 9 dimensions at most!');
            end
            mask = mask & tempmask;
        end
    end    
    column_indexes = sels(mask);
end

% Return a cell array of indexes which can be used in a for loop.
function column_indexes = foreach(dname)
    
    % Ensure that dname is a field
    if ~isempty(find(dimension_names, dname))
        error('Unknown field name: %s', dname);
    end
    
    % Build up a list of indexes
    column_indexes = {};
    idx = find(dimension_names, dname);
    jj = 1;
    for ii = 1:dimension_sizes(idx)
        column_indexes{jj} = select_cols_fn(dname, ii);
        jj = jj + 1;
    end
    
    sels = reshape(sels, dimension_sizes);   
    mask = true(size(sels));

end

end