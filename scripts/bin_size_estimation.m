function bin_size_estimization(batch, cellids, modelnames)

global STACK XXX META MODULES;

if length(modelnames) > 1
    error('This script only works when given a single model and one or more cellids.');
end

freqs = logspace(0, log10(1000), 50); % 50 points between 1 and 1000Hz
stimfiles = {};
c = [];

    function results = estimate_bin_size(~)
        results = [];
        
        for ii = 1:length(freqs)
            freq = freqs(ii);
            STACK{1}{1}.raw_stim_fs = 100;
            STACK{1}{1}.raw_resp_fs = freq;    
            calc_xxx(1,2);
        
            sfs = fieldnames(XXX{2}.dat);
            for jj = 1:length(sfs),
                sf = sfs{jj};
                r = XXX{2}.dat.(sf).resp;
                reps = size(r, 3);        
                stimfiles{end+1} = [sf sprintf('(%d reps)', reps)];        
                ki = nansum(reshape(r, [], size(r, 3)), 2);
                n = length(ki);
                delta = 1 / freq;
                m = mean(ki);
                v = sum((ki - m).^2) / n;
                score(ii, jj) = (2*m - v) / delta^2;    
            end    
        end        
        c = cat(2, c, score);    
    end
        
nullfn = @(x) 0;
[~, ~, ~] = load_model_batch(batch, cellids, modelnames, ...
            nullfn, @estimate_bin_size);

% 
% %STACK = {};
% XXX = {};
% XXX{1}.cellid = 'sti019c-d1';
% XXX{1}.training_set = {'sti019c06_p_VOC'};
% XXX{1}.test_set = {'sti019c05_p_VOC'};
% XXX{1}.filecodes = {};

%append_module(MODULES.load_stim_resps_from_baphy)
% 
% freqs = logspace(0, log10(1000), 50); % 50 points between 1 and 250Hz
% stimfiles = {};
% c = [];
% for ii = 1:length(freqs)
%     freq = freqs(ii);
%     STACK{1}{1}.raw_stim_fs = 100;
%     STACK{1}{1}.raw_resp_fs = freq;    
%     calc_xxx(1,2);
%     
%     sfs = fieldnames(XXX{2}.dat);
%     for jj = 1:length(sfs),
%         sf = sfs{jj};
%         r = XXX{2}.dat.(sf).resp;
%         reps = size(r, 3);        
%         stimfiles{end+1} = [sf sprintf('(%d reps)', reps)];        
%         ki = nansum(reshape(r, [], size(r, 3)), 2);
%         n = length(ki);
%         delta = 1 / freq;
%         m = mean(ki);
%         v = sum((ki - m).^2) / n;
%         c(ii, jj) = (2*m - v) / delta^2;    
%     end    
% end


% Normalize C since we just care about the minimum
c_min = repmat(min(c), [size(c, 1), 1]);
c_max = repmat(max(c), [size(c, 1), 1]);
c_scale = c_max - c_min;
c = (c - c_min) ./ c_scale;

figure;
%plot(freqs, c);
semilogx(freqs, c);
title('Choosing Optimal Bin Size (Lower is Better)')
xlabel('Sampling Frequency [Hz]');
ylabel('Normalized Mean Integrated Squared Error (NMISE)');
l = legend(stimfiles, 'Interpreter', 'none');
end
