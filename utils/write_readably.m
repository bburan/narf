function s = write_readably(obj)
% s = write_readably(obj)
%
% In languages like LISP, there is a concept that you should be able to 
% write a data structure such that it can be "eval'd" again by the program
% language interpreter. This is a hacky way of achieving the same thing in
% matlab for some (but certainly not all) non-recursive data structures. 
% 
% ARGUMENTS:
%    obj     A string, vector, matrix, cell array, or function handle.
%
% RETURNS:
%    s       A string representation of obj that can be converted into obj
%            using eval(s).
%
if ismatrix(obj) & isnumeric(obj) & any(size(obj) ~= 1)  % Matrices
     if isempty(obj)
         s = '[]';
     else       
        s = '[';
        [M, N] = size(obj);
        for ii = 1:M-1
            s = strcat(s, num2str(obj(ii,:)), '; ');
        end
        s = strcat(s, num2str(obj(M,:)), ']');
     end
elseif iscell(obj) 
    if isempty(obj)
         s = '{}';
    else         
        s = '{';
        for ii = 1:length(obj)
            s = strcat(s, write_readably(obj{ii}), ', ');
        end
        s = strcat(s, '}');   
    end
elseif isstr(obj)       % Single strings must be quoted
    s = ['''' obj '''']; 
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