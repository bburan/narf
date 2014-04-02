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
               
        WIN = ceil(mdl.input_freq / mdl.output_freq); % Window Size        
        N = 2^nextpow2(WIN);  % Number of FFT points
        NOV = WIN/2;        % Overlap points       
        SR = mdl.input_freq;  % Source sample rate
        FMIN = mdl.min_freq;  % Minimum frequency to display
        FMAX = SR/2;
        NCHANS = mdl.n_output_chans;
        fratio = (FMAX/FMIN)^(1/NCHANS); % Ratio between adjacent freqs in log-f axis
            
        % Freqs corresponding to each bin in FFT
        fftfrqs = [0:(N/2)]*(SR/N);
        nfftbins = N/2+1; 
        
        logffrqs = FMIN .* fratio.^[0:(NCHANS-1)];  % Bin frequencies, log scale
        logfbws = logffrqs * (fratio - 1); % Bandwidths of each bin in log F
        logfbws = max(logfbws, SR/N);      % Bandwidth cannot be less than FFT binwidth    

        % Control overlap between adjacent bands
        ovfctr = 0.5475;   % Adjusted by hand to make sum(mx'*mx) close to 1.0
                
        % Weighting matrix mapping energy in FFT bins to logF bins
        % is a set of Gaussian profiles depending on the difference in
        % frequencies, scaled by the bandwidth of that bin        
        freqdiff = ( repmat(logffrqs',1,nfftbins) - repmat(fftfrqs,NCHANS,1) )./repmat(ovfctr*logfbws',1,nfftbins);
        mx = exp( -0.5*freqdiff.^2 );
        % Normalize rows by sqrt(E), so multiplying by mx' gets approx orig spec back
        mx = mx ./ repmat(sqrt(2*sum(mx.^2,2)), 1, nfftbins);
                  
        % err = sum(mx'*mx)
        % Try to normalize MX a little better
        mx = mx .* repmat(1./(sqrt(sum(mx'*mx))), size(mx,1), 1);
        
        tmp = [];
        
        for s = 1:S                               
            trial = sig(:,s,1);
            XX = spectrogram(trial,WIN,NOV,N,SR);            
            y = sqrt( mx * (abs(XX).^2) );            
            tmp(:,s,:) = y';
        end
        
        x.dat.(sf).(mdl.output) = tmp;
        x.dat.(sf).(mdl.output_time) = (1/mdl.output_freq) * [0:size(tmp,1)-1]';
    end
end

function do_plot_fft_spectrogram(sel, stack, xxx)
    mdl = stack{end}{1};
    x = xxx{end};                    

    fratio = (mdl.input_freq/mdl.min_freq)^(1/mdl.n_output_chans);
    logffrqs = mdl.min_freq .* fratio.^[0:(mdl.n_output_chans-1)];  % Bin frequencies, log scale
    
    X = squeeze(x.dat.(sel.stimfile).(mdl.input)(:, sel.chan_idx, 1));
    y = squeeze(x.dat.(sel.stimfile).(mdl.output)(:, sel.chan_idx, :))';
    imagesc([0 length(X)/mdl.input_freq],[1 mdl.n_output_chans], 20*log10(y));
    axis xy;
    caxis([-20,40]);  % TODO: use a 'smarter' caxis here
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    axis tight;  
    yt = get(gca,'YTick');
    for i = 1:length(yt)
       ytl{i} = sprintf('%.0f',logffrqs(yt(i)));
    end
    set(gca,'YTickLabel',ytl);
end


end