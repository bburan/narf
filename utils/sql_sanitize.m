function string = sql_sanitize(string)
% A function to sanitize a string and make it safe for putting in a DB
% Escapes the characters \x00, \n, \r, \, ', " and \x1a. 

if ismatrix(string)
    s = '';
    for ii = 1:size(string, 1)
        s = cat(2, s, strtrim(string(ii,:))); %% FIXME: Right now this removes newlines
    end
    string = s;
end

string = regexprep(string, '''', '\\''');
string = regexprep(string, '"', '\\"');

string = regexprep(string, sprintf('\n'), '\\n');
string = regexprep(string, sprintf('\r'), '\\r');

% string = regexprep(string, sprintf(','), '\\,'); % Needed?
%string = regexprep(string, sprintf('\t'), '\\t'); % Needed?
%string = regexprep(string, sprintf('\'), '\\%'); % Needed?

% TODO: I'm not sure how to escape x00, x1a, although they can be used for
% SQL injection attacks apparently. 

end