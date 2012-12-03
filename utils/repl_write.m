function s = repl_write(obj)
% Prints obj in a readable manner. 
if ismatrix(obj) & isnumeric(obj) & any(size(obj) ~= 1)  % Matrices
    s = '[';
    [M, N] = size(obj);
    for ii = 1:M-1
        s = strcat(s, num2str(obj(ii,:)), '; ');
    end
    s = strcat(s, num2str(obj(M,:)), ']');
elseif isstr(obj)       % Single strings
    s = obj; 
elseif isnumeric(obj)   % Single numbers
    s = num2str(obj);
elseif islogical(obj) % Booleans
    if obj
        s = 'true'; 
    else
        s = 'false';
    end
elseif isa(obj, 'function_handle')
    s = ['@' func2str(obj)];
else
    error('Not sure how to print: %s', obj);
end