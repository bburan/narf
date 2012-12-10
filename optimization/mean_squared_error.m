function z = mean_squared_error(w)
% Returns the mean squared error of lf_stim and raw_resp
global STACK XXX;

% Unpack the vector and set the stack up to reflect it
unpack_fittables(w);

% Recalculate the stack, starting at the needed point
start_depth = find_fit_start_depth(STACK);
recalc_stack(start_depth);

% Concatenate everything together
x = XXX{end};
V1 = [];
V2 = [];
for sf = fieldnames(x.dat)', sf = sf{1};
    [S, T] = size(x.dat.(sf).lf_stim);
    V1 = cat(1, V1, reshape(x.dat.(sf).raw_respavg',[],1));
    V2 = cat(1, V2, reshape(x.dat.(sf).lf_stim',[],1));
end

% Return the mean squared error
z = mean(abs(V1 - V2));