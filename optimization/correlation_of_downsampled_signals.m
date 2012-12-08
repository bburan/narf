function z = correlation_of_downsampled_signals(w)
% Returns the correlation of lf_stim and raw_resp
global STACK XXX;

% Unpack the vector and set the stack up to reflect it
unpack_fittables(w);

% Recalculate the stack, starting at the needed point
start_depth = find_fit_start_depth(STACK);
recalc_stack(start_depth);

% Compute correlation after concatenating everything together
x = XXX{end};
V1 = [];
V2 = [];
for sf = fieldnames(x.dat)', sf = sf{1};
    [S, T] = size(x.dat.(sf).lf_stim);
    V1 = cat(1, V1, reshape(x.dat.(sf).raw_respavg',[],1));
    V2 = cat(1, V2, reshape(x.dat.(sf).lf_stim',[],1));
end
R = corrcoef(V1,V2);
R(isnan(R)) = 0; % corrcoef returns NaNs if FIR had all zero coefficients
z = R(2,1)^2;  % Return r^2
