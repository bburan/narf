function s = repl_write(obj)
% Prints obj in a readable manner. 
if ismatrix(obj) & isnumeric(obj) & any(size(obj) ~= 1)  % Matrices
    s = strcat('[', num2str(obj), ']');
    s = regexprep(s, '\n', '; ');
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