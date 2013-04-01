function c = pickcolor(n)
% c = pickcolor(n)
%
% Returns a color and linestyle for a particular index. 
% Useful when 'hold on' is set and you need to plot different lines with
% repeated calls to PLOT() but want them to be the proper colors.
%
% Essentially, this function is a very clumsy way of accomplishing this:
%
% set(gca, 'ColorOrder', [0 0 1; 0 1 0; 1 0 0; 0 1 1; 1 0 1; 1 1 0;], ...
%         'LineStyleOrder',{'-','--',':'},'NextPlot','ReplaceChildren');

m = mod(n,19);
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
    case 7
        c='b--';
    case 8 
        c='g--';
    case 9 
        c='r--';
    case 10 
        c='c--';
    case 11
        c='m--';
    case 12
        c='y--';
    case 13
        c='b:';
    case 14
        c='g:';
    case 15
        c='r:';
    case 16 
        c='c:';
    case 17
        c='m:';
    case 18
        c='y:';      
end