function flatstack = get_flatstack()
global STACK;

flatstack={};
for ii = 1:length(STACK)
    for jj = 1:length(STACK{ii})
        flatstack{end+1} = STACK{ii}{jj};
    end
end
