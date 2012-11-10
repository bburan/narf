function c = pickcolor(n)
m = mod(n,7);
switch m
    case 0 
        c='k-';
    case 1 
        c='b-';
    case 2 
        c='g-';
    case 3 
        c='r-';
    case 4 
        c='c-';
    case 5
        c='m-';
    case 6
        c='y-';
end
