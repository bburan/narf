function verify_model_polarity()
% verify_model_polarity()
%
% (Requires that the STACK and XXX data structures be fully initialized.)
%
% For each FIR filter in the STACK, 
%     Flips the sign of every element in that filter iff its output is
%     creating a scatter plot with a negative slope.
%
% No arguments or return values.

global STACK XXX;

[~, firmod_idxs] = find_modules(STACK, 'fir_filter');

for kk = 1:length(firmod_idxs)
    idx = firmod_idxs{kk};
    
    firmod = STACK{idx}{1};
    x = XXX{idx+1};
        
    % Flip the sign of the coefs if there is a negative slope
    V1 = [];
    V2 = [];
    for ii = 1:length(x.training_set),
        sf = x.training_set{ii};
        V1 = cat(1, V1, x.dat.(sf).(firmod.output)(:));
        V2 = cat(1, V2, x.dat.(sf).respavg(:));
    end
   
    [~, idxs, ~] = unique(V1, 'first'); % Remove duplicates
    D = [V1(idxs) V2(idxs)]; 
    D = excise(D);
    D = sortrows(D);
    D = conv_fn(D, 1, @nanmean, ceil(size(D, 1)/100), 0);
    xs = D(:,1);
    ys = D(:,2);
    
%     % Maybe filtering will help remove outliers?
% 	winsize = ceil(length(ys) / 4);
%     if mod(winsize, 2) == 1
%         winsize = winsize + 1;
%     end
%     gf = gausswin(winsize);
%     gf = gf / sum(gf);
%     %x = filter(gf, 1, xs);
%     y = filter(gf, 1, ys);    
    idx10 = ceil(length(ys) / 20);
    y1 = mean(ys(1:idx10));
    y2 = mean(ys(end-idx10:end));
    if (y1 > y2)
        fprintf('verify_model_polarity.m: Detected negative slope. Flipping.\n');
        for aa = 1:length(STACK{idx})
            STACK{idx}{aa}.coefs = - STACK{idx}{aa}.coefs;
        end
    end
    
    % Possible 4th attempt: Find the index of the peak point
    % If it is to the 'left' of the center of the range, flip?
    
    % Debug crap
    figure; plot(xs, ys, 'b.');
    figure; imagesc(STACK{idx}{1}.coefs);
    ca = caxis;
    lim = max(abs(ca));
    caxis([-lim, +lim]);
    keyboard
    
    % OLD way wasn't sufficiently reliable: just linear isn't enough?
    % P = polyfit(D(:,1),D(:,2), 1);
    % if (P(1) < 0)
    %     for aa = 1:length(STACK{idx})
    %         STACK{idx}{aa}.coefs = - STACK{idx}{aa}.coefs;
    %     end
    % end
    
    % Another way: try to fit a sigmoid to the curve 
%     xs = D(:,1);
%     ys = D(:,2);
%     phi_init1 = [mean(xs), var(xs), 0.5*(max(xs)-min(xs)), 0]; 
%     phi_init2 = [- mean(xs), var(xs), 0.5*(max(xs)-min(xs)), 0];    
%     opts = optimset('Display','off');
%     [phi_best1, rn1] = lsqcurvefit(@nl_sigmoid, ...
%                            phi_init1, xs, ys, [], [], opts);
%     [phi_best2, rn2] = lsqcurvefit(@nl_sigmoid, ...
%                            phi_init2, xs, ys, [], [], opts);
%     if (rn2 < rn1)
%         for aa = 1:length(STACK{idx})
%             STACK{idx}{aa}.coefs = - STACK{idx}{aa}.coefs;
%         end
%     end   
    
    
    calc_xxx(idx); 
    
end
