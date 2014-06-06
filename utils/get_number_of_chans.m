function num_chans = get_number_of_chans(XXX, signal)
% get_number_of_chans(XXX)
%
% Returns the number of channels at the end of the stack for signal

x = XXX{end};
fns = fieldnames(x.dat);
sf = fns{1};
[T, S, num_chans] = size(x.dat.(sf).(signal));            