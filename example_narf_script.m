% Database query: select just one very good cellfile for which to generate plots
[cfd, cellids, cellfileids] = dbgetscellfile('cellid','dai008a-c1','rawid', 64662);
index = 1;
channel = 1;

% ------------------------------------------------------------------------
% THINGS ABOVE THIS LINE SELECT DATA, THINGS BELOW IT DEFINE THE MODEL

% Create a map from a parameter struct to a vector theta
p2theta = @(p) [p.prefiltercenterfs,    % 1
                p.linearcoefficients,   % depends on length 
                ];
            
% Update user-readable struct from machine-usable parameter vector theta
function p = updatep(p, theta)  
    p.stimfs = 20000;      % Stimulus frequency
    p.respfs = 200;        % Response histogram bin frequency (200Hz -> 5ms)   
    
    % A function to create the prefilter:
    p.prefilterfn = @(S, Sfs) gammatone.gammatone(S, Sfs, theta(1), p.prefilterphasealign); 
    
    % A function to create the (non)linear model
    p.modelfn     = @(theta) firfilter.firfilter();
    
end

% Define the initial values
p=[];

p.stimfs = 20000;      % Stimulus frequency
p.respfs = 200;        % Response histogram bin frequency (200Hz -> 5ms) 

% Prefilter
p.prefiltercenterfs   = 8000;    % Center frequency [Hz]
p.prefilterbandwidth  = 1000;    % Bandpass width [Hz] 
p.prefilterphasealign = false;   % Align the gamma filters to be in phase?
p.prefiltersmoothfs   =  200;    % Smoothing frequency [Hz]

% Goodness of fit function
function r = jackcorr(s1, s2)
    r=0.012345;
end
p.evaluatefn = @r;

% Convert that into a vector
% TODO


% Define the regularization parameters: limits, priors, and penalties 
% TODO

% ------------------------------------------------------------------------
% THERE SHOULD BE NO NEED TO EDIT BELOW THIS LINE IF THE FRAMEWORK IS GOOD

% Prepare DB options for requesting the stimulus and response files
options = []; 
options.includeprestim = 0;  % Include pre and post stimulation results
options.unit     = cfd(index).unit;
options.channel  = channel;
options.rasterfs = cfg.respfs; 

% Load and raster the auditory stimulus
fprintf('Loading stimulus file: %s%s\n', cfd(index).stimpath, cfd(index).stimfile);
stimfile = [cfd(index).stimpath cfd(index).stimfile];
stim     = loadstimfrombaphy(stimfile, [], [], 'wav', cfg.stimfs, 1, 0, options.includeprestim);
stimenv  = loadstimfrombaphy(stimfile, [], [], 'envelope', options.rasterfs, 1, 0, options.includeprestim);

% Load and raster the neural spike response
fsprintf('Loading response file: %s/%s\n', cfd(index).path, cfd(index).respfile);
respfile = [cfd(index).path cfd(index).respfile]; 
fsprintf('Rastering response file...\n');
[resp, tags] = loadspikeraster(respfile, options); 

% Create memoized versions of functions so that results are cached automatically
% TODO
[N_vals, N_reps, N_stims] = size(stim);

% ------------------------------------------------------------------------
% Begin the regression procedure

% TODO: This should also be iterated in the context of the best model
fprintf('Computing the Expectation-Maximizing spike rate average\n');
[Yhat, Ytop, Ybot] = KalmanPSTH(resp);

% initialize the optimization memory vector
psi = [ ];

% Use Matlab's impulse estimation tools to infer the impulse response

% Loop until the termination condition Omega is satisfied
%while (Omega(Y, Yhat, psi) == false)
    
    % Run the prefilter on the stimulus      
    [Z, Zfs] = p.prefilterfn(stim, stimfs);
    
    % Run Kalman filter on the prefiltered stimulus 
    Y = p.prefilterfn(Z, Zfs);
    
    % Check the correlation coefficient
    %R=corrcoef(resp', stimenv);
    %disp(sprintf('Correlation psth-env:          %f', R(2,1)));
    
    % Update theta
%end

% Save the best model params somewhere safe

% Plot the stimulus, response, and best-fitting model
fprintf('Generating stimulus/response plot...\n');                   
k=1;

t_hop = 1/resp_fs;    % Time between gamma window centers [seconds]
N_gfs = 64;           % Number of gamma tone filter channels
f_min = 100;          % Lowest frequency [Hz]
f_max = stim_fs/2;    % Highest (nyquist) frequency [Hz]
N_xticks = 10;        % Number of ticks to display on the x axis
N_yticks = 5;         % Number of ticks to display on the y axis
fh = plotutil.stimrespplot(stim(1,:,k), stim_fs, sum(resp(:,:,k),2), options.rasterfs, stimenv(1,:,k));

% Plot the correlation of the model, and error residuals

% TODO: Any division stuff should avoid divide by zeros (smoothfs, for
% example)