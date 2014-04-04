% A demo of nonlinearities

x = [0:0.001:1]';

d = {};

% Single parameter models
d{1} = [];
d{1}.lin   = {[1 0 0], x};
d{1}.logn =  {[0 1 0], nl_log([-.3], x)};
d{1}.rootn = {[0 0 1], nl_root([2], x)};
d{1}.dlog =  {[0 1 1], nl_dlog([-.3], x)};

% Two parameter models
d{2}.poly1 =   {[1 0 0], polyval([1 0.2], x)};
d{2}.lognb =   {[0 1 0], nl_log([-.3 .2], x)};
d{2}.rootnb =  {[0 0 1], nl_root([2 .1], x)};
d{2}.dlogb  =  {[0 1 1], nl_dlog([-.3 .2], x)};
d{2}.invlin =  {[1 0 1], nl_invlin([.2 2], x)};
d{2}.exp     = {[0 0 0], nl_exponential([.5 3], x)};

% Three parameter models
d{3}.poly2 =   {[1 0 0], polyval([-.4 1 0.2], x)};
d{3}.lognbz =  {[0 1 0], nl_log([-.3 .2 .1], x)};
d{3}.rootnbz = {[0 0 1], nl_root([2 .2 .1], x)};
d{3}.dlogbz =  {[0 1 1], nl_dlog([-.3 .2 .1], x)};
d{3}.invlinz = {[1 0 1], nl_invlin([.2 2 .1], x)};
d{3}.zthresh = {[0 0 0], nl_zerothresh([.1 2 .2], x)};

% Four parameter models
d{4}.poly3 = {[1 0 0], polyval([1 -1 1 0.2], x)};
d{4}.sig   = {[0 1 0], nl_sigmoid([.5 .1 .5 .2], x)};
d{4}.dexp =  {[0 0 1], nl_dexp([.2 1 .5 10], x)};
d{4}.zsoft   = {[0 1 1], nl_softzero([.2 20 0.5 .2], x)};

% Five parameter models
d{5}.sigrich   = {[1 0 0], nl_richards([.1 1 .5 30 1], x)};
d{5}.sigerf    = {[0 1 0], nl_sig_erf([.1 .5 .5 0.2 0.1], x)};
d{5}.sigell    = {[0 0 1], nl_sig_elliot([.1 1 .5 5 10], x)}; 
d{5}.siglog    = {[0 0 0], nl_sig_logistic([.1 1 .5 15 20], x)};
d{5}.sigumbel  = {[1 0 1], nl_sig_gumbel([.1 1 .5 0.2 0.1], x)}; 
d{5}.sigumber  = {[0 1 1], nl_sig_gumber([.1 1 .5 0.2 0.1], x)}; 
d{5}.sigcauchy = {[.5 .5 1], nl_sig_cauchy([.1 1 .5 5 10], x)}; 

figure();
N = length(d);
for jj = 1:length(d)    
    names = fieldnames(d{jj});
    leg = {};
    subplot(1,N,jj)
    hold on;
    for ii = 1:length(names)
        name = names{ii};
        plot(x, d{jj}.(name){2}, 'Linewidth', 2, 'Linestyle', '-', 'color', d{jj}.(name){1});
        leg{end+1} = name;
        axis([0 1 0 1]);        
    end
    legend(leg, 'location', 'SouthOutside');
    xlabel('Input');
    ylabel('Output');
    title(sprintf('%d Parameter Curves', jj));    
    axis square;
    hold off;
end
