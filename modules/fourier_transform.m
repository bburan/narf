function m = fourier_transform(args)
% Fourier Transform NPNL. Returns many channels, logarithmically spaced
% along the frequency axis. 
%
% Returns a function module 'm' which implements the MODULE interface.
% See documentation for more information on how modules are typically used.

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @fourier_transform;
m.name = 'fourier_transform';
m.fn = @do_fourier_transform;
m.pretty_name = 'Fourier Transform';
m.editable_fields = {'n_output_chans', 'min_freq', 'input', 'input_time', ...
                     'input_freq', 'output', 'output_time', 'output_freq'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.n_output_chans = 128; 
m.min_freq = 200;
m.max_freq = 20000;
m.input = 'stim';
m.input_freq = 50000;
m.input_time = 'stim_time';
m.output = 'stim';
m.output_freq = 200;
m.output_time = 'stim_time';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_fft_spectrogram;
m.plot_fns{1}.pretty_name = 'Output Spectrogram';

function x = do_fourier_transform(mdl, x, stack, xxx)
    sfs = fieldnames(x.dat);    
    for ii = 1:length(sfs)
        sf = sfs{ii};        
        
        sig = x.dat.(sf).(mdl.input);
        [T,S,C] = size(sig);
        
        if C~=1
            error('I can only do an FFT on one channel!');
        end
               
        blocksize = ceil(mdl.input_freq / mdl.output_freq);
        overlap = blocksize/2;
        NFFT = 2^nextpow2(blocksize+2*overlap);
        
        tmp = zeros(ceil(T/blocksize), S, mdl.n_output_chans);       
        
        for s = 1:S          
            trial = sig(:,s,1);
            [~, F, TT, P] = spectrogram(trial, blocksize, overlap, NFFT, mdl.input_freq);
            mask = (F>mdl.min_freq) & (F<mdl.max_freq);
            region = log10(abs(P(mask,:)).^2);
            out = imresize(region, [mdl.n_output_chans, ceil(T/blocksize)], ...
                'Method', 'bicubic', 'Dither', true, 'Antialias', true);
            tmp(:,s,:) = out';
        end
        
        x.dat.(sf).(mdl.output) = tmp;
        x.dat.(sf).(mdl.output_time) = (1/mdl.output_freq) * [0:size(tmp,1)-1]';
        x.dat.(sf).freqmapthingy = F(mask);
    end
end

function do_plot_fft_spectrogram(sel, stack, xxx)
    mdl = stack{end}{1};
    x = xxx{end};                    

    fratio = (mdl.max_freq/mdl.min_freq)^(1/mdl.n_output_chans);
    logffrqs = mdl.min_freq .* fratio.^[0:(mdl.n_output_chans-1)];  % Bin frequencies, log scale
    
    X = squeeze(x.dat.(sel.stimfile).(mdl.input)(:, sel.chan_idx, 1));
    y = squeeze(x.dat.(sel.stimfile).(mdl.output)(:, sel.chan_idx, :))';
    imagesc([0 length(X)/mdl.input_freq],[1 mdl.n_output_chans], y);
    axis xy;
    %caxis([-20,40]);  % TODO: use a 'smarter' caxis here
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    axis tight;  
    yt = get(gca,'YTick');
    for i = 1:length(yt)
       %ytl{i} = sprintf('%.0f',   logffrqs(yt(i)));
       ytl{i} = sprintf('%.0f',  x.dat.(sel.stimfile).freqmapthingy(yt(i)));
    end
    set(gca,'YTickLabel',ytl);
end


end