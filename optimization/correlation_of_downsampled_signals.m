function z = correlation_of_downsampled_signals(w)
% Returns the correlation of lf_stim and raw_resp
global STACK XXX;

% Unpack the vector and set the stack up to reflect it
unpack_fittables(w);

% -------------
% TODO: Move this block to a reusable function, since everybody will do
% this over and over again.

% Recalculate the stack, starting at the needed point
start_depth = find_fit_start_depth(STACK);

% If data already exists, invalidate it just in case
if start_depth < length(XXX)
    XXX = XXX(1:start_depth);  % Invalidate later data so it cannot be used
end

% If there is not enough data, begin the calculation from the top
if start_depth > length(XXX)
    start_depth = length(XXX);
end

% Now, do the recalculation of the data
for ii = start_depth:length(STACK);
    if ~STACK{ii}.isready_pred(STACK(1:ii), XXX(1:ii));
        error('Stack was not fully ready at depth %d', ii);
    end
    XXX{ii+1} = STACK{ii}.fn(STACK(1:ii), XXX(1:ii));
end
% ------------------

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
