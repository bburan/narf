function m = pop_module()
% POP_MODULE()
%
% Pops the last module off of the current current global cell array STACK
% and returns it. 

global STACK XXX;
l = length(STACK);
m = STACK{l};
STACK = STACK(1:l-1);
XXX = XXX(1:l);
