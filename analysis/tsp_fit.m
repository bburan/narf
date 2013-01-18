function tsp_fit(cellid)

% Fit a model to five stimulus files to initialize model params.
% Fit a specific model parameter again using each file individually.
% Save the fitted params and plot the results.

if ~exist('dbgetscellfile','file'), baphy_set_path; end
if ~exist('fir_filter','file'), narf_set_path; end
global NARF_PATH STACK XXX;

close all

files_to_analyze = {...
    % OLD SINGLE-SPEAKER RECORDINGS
    
    
    %{'por020a-c1',  'por020a04_p_SPN',	'por020a08_a_TSP',  'por020a09_p_SPN',	'por020a10_a_TSP', 	'por020a11_p_SPN'};
    %{'por020a-d1',  'por020a04_p_SPN',	'por020a08_a_TSP',  'por020a09_p_SPN',	'por020a10_a_TSP', 	'por020a11_p_SPN'}; % Crummy STRF. Not good fit!
    %{'por019b-a1',  'por019b02_p_SPN',	'por019b05_a_TSP',  'por019b07_p_SPN', 	'por019b08_a_TSP', 	'por019b09_p_SPN'}; % Crash
    %{'por018b-a1',  'por018b04_p_SPN',	'por018b06_a_TSP',  'por018b07_p_SPN',	'por018b08_a_TSP',	'por018b09_p_SPN'}; % Crash 
    %{'por018b-c1',  'por018b04_p_SPN',	'por018b06_a_TSP',  'por018b07_p_SPN',	'por018b08_a_TSP',	'por018b09_p_SPN'}; % Stable one. Crash.
    %{'por018b-d1',  'por018b04_p_SPN',	'por018b06_a_TSP',  'por018b07_p_SPN',	'por018b08_a_TSP',	'por018b09_p_SPN'}; % Crash
    %{'por018c-a1',  'por018c05_p_SPN',	'por018c06_a_TSP',  'por018c08_p_SPN',	'por018c09_p_TSP',	'por018c10_p_SPN'};
    %{'por018c-a2',  'por018c05_p_SPN',	'por018c06_a_TSP',  'por018c08_p_SPN',	'por018c09_p_TSP',	'por018c10_p_SPN'};
    %{'por018c-c1',  'por018c05_p_SPN',	'por018c06_a_TSP',  'por018c08_p_SPN',	'por018c09_p_TSP',	'por018c10_p_SPN'}; % Crash
    
    
    % TWO-SPEAKER RECORDINGS
    
    {'por022a-a1',  'por022a12_p_SPN', 	'por022a13_a_TSP',	'por022a14_p_SPN'};  
    {'por022a-c1', 'por022a08_p_SPN'	'por022a09_a_TSP', 	'por022a12_p_SPN',	'por022a13_a_TSP',	'por022a14_p_SPN'};
    {'por022b-a1', 'por022b08_p_SPN', 	'por022b10_a_TSP', 	'por022b11_p_SPN'}; 
    {'por022b-a2', 'por022b08_p_SPN',	'por022b10_a_TSP', 	'por022b11_p_SPN'}; 
    {'por023a-a1', 'por023a06_p_SPN', 	'por023a07_a_TSP',  'por023a09_p_TSP',	'por023a10_p_SPN'};                     
    %%{'por023a-b1', 'por023a06_p_SPN', 'por023a07_a_TSP', 	'por023a09_p_TSP',	'por023a10_p_SPN'};                                         
    {'por023b-a1', 'por023b12_p_SPN', 'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'}; 
    {'por023b-b1', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'}; 
    {'por023b-d1', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN'}; 
    {'por023b-d2', 'por023b12_p_SPN', 	'por023b13_a_TSP', 	'por023b15_p_SPN',	'por023b16_a_TSP', 	'por023b18_p_SPN'}; 
    {'por024a-c1', 'por024a09_a_SPN',  'por024a13_a_TSP',  'por024a14_a_TSP',  'por024a15_p_SPN'};                      
    {'por025a-b1', 'por025a09_p_SPN', 	'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'};
    %%{'por025a-c1', 'por025a09_p_SPN', 'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'}; 
    %%{'por025a-d1', 'por025a09_p_SPN', 'por025a10_a_TSP', 	'por025a11_p_SPN', 	'por025a12_a_TSP', 	'por025a13_p_SPN'};
    {'por025c-c2', 'por025c04_p_SPN',	'por025c05_a_TSP', 	'por025c06_p_SPN', 	'por025c07_a_TSP', 	'por025c08_a_SPN'}; % crash
    %%{'por026a-b1', 'por026a08_p_SPN', 'por026a09_a_TSP', 	'por026a10_p_SPN', 	'por026a11_a_TSP',	'por026a12_p_SPN'};
    %%{'por026a-d1', 'por026a08_p_SPN', 'por026a09_a_TSP', 	'por026a10_p_SPN'};                                         % Better w/o these two 'por026a11_a_TSP','por026a12_p_SPN'
    {'por027a-a1', 'por027a06_p_SPN',	'por027a07_a_TSP',	'por027a08_p_SPN',	'por027a09_a_TSP',	'por027a13_p_SPN'}; 
    {'por026b-b1', 'por026b11_p_SPN',  'por026b12_a_TSP',  'por026b13_p_SPN'};  
    {'por027a-a1', 'por027a06_p_SPN',	'por027a07_a_TSP',	'por027a08_p_SPN'};
    %%{'por027a-b1', 'por027a06_p_SPN',	'por027a07_a_TSP',	'por027a08_p_SPN'}; 
    
    %{'por027a-c1', 'por027a06_p_SPN',	'por027a07_a_TSP',	'por027a08_p_SPN',	'por027a09_a_TSP',	'por027a13_p_SPN'};
    %%{'por027b-b1', 'por027b06_p_SPN', 'por027b08_a_TSP',	'por027b11_p_SPN',	'por027b14_a_TSP',	'por027b16_p_SPN'}; 
    {'por027b-b2', 'por027b06_p_SPN', 	'por027b08_a_TSP',	'por027b11_p_SPN',	'por027b14_a_TSP',	'por027b16_p_SPN'};                                         
    {'por027b-c1', 'por027b06_p_SPN', 	'por027b08_a_TSP',	'por027b11_p_SPN',	'por027b14_a_TSP',	'por027b16_p_SPN'}; 
    {'por028b-a1', 'por028b02_p_SPN',	'por028b05_a_TSP',	'por028b06_p_SPN',	'por028b10_a_TSP',	'por028b12_p_SPN'}; 
    {'por028b-b1', 'por028b02_p_SPN',	'por028b05_a_TSP',	'por028b06_p_SPN',	'por028b10_a_TSP',	'por028b12_p_SPN'}; 
    {'por028b-b2', 'por028b02_p_SPN',	'por028b05_a_TSP',	'por028b06_p_SPN',	'por028b10_a_TSP',	'por028b12_p_SPN'}; 
    {'por028b-c1', 'por028b02_p_SPN',	'por028b05_a_TSP',	'por028b06_p_SPN',	'por028b10_a_TSP',	'por028b12_p_SPN'};
    {'por028b-d1', 'por028b02_p_SPN',	'por028b05_a_TSP',	'por028b06_p_SPN',	'por028b10_a_TSP',	'por028b12_p_SPN'}; 
    };

% -------------------------------------------------------------------------
% Define the model
mdls = scan_directory_for_modules([NARF_PATH filesep 'modules']);

raster_fs = 100;
filter_length = 10;
n_channels = 2;

% -------------------------------------------------------------------------
ci=0;
for ii=1:length(files_to_analyze),
    if strcmpi(cellid,files_to_analyze{ii}{1}),
        ci=ii;
    end
end

if ~ci,
    disp('cellid not found in list');
end

savepath = [NARF_PATH filesep 'saved_models' filesep cellid];
if ~exist(savepath,'dir'),
    mkdir(fileparts(savepath),cellid);
end
filecount=length(files_to_analyze{ci});
respfiles = {files_to_analyze{ci}{2:filecount}};


% Re-initialize the data to evaluate
XXX = {};
XXX{1}.cellid = cellid;  
XXX{1}.training_set = respfiles;
XXX{1}.test_set = {};
    
STACK = {};
append_module(mdls.load_stim_resps_from_baphy.mdl(...
    struct('raw_resp_fs', raster_fs, ...
           'raw_stim_fs', raster_fs,...
           'include_prestim',0,...
           'stimulus_format','envelope', ...
           'stimulus_channel_count', n_channels)));

%append_module(mdls.depression_filter_bank.mdl(struct('num_channels',1,...
%                                                  'strength',0.2,...
%                                                  'tau',50)));
%append_module(mdls.normalize_channels);
append_module(mdls.normalize_channels);
append_module(mdls.fir_filter.mdl(struct('num_coefs',10,...
                                         'fit_fields','coefs')));
append_module(mdls.correlation);

%filename = sprintf('%s/%s_all_xccorefet.mat', savepath, cellid);
%save_model_stack(filename, STACK, XXX);

% re-fit fir filter using least-square
firmodidx = 3; % find_module(STACK, 'fir_filter');
STACK{firmodidx}.fit_fields = {'coefs','baseline'};
fit_with_lsqcurvefit();

% re-fit fir filter using least-square
STACK{firmodidx}.fit_fields = {'coefs'};
fit_with_lsqcurvefit();

% add output NL to stack
sf=respfiles{1};
resp=XXX{5}.dat.(sf).respavg;
linpredmin=min(XXX{5}.dat.(sf).stim(:));
linpredmax=max(XXX{5}.dat.(sf).stim(:));

phi1=mean(std(XXX{4}.dat.(sf).stim(:)))-std(XXX{4}.dat.(sf).stim(:))*2;
phi2=0.001;
phi3=max(resp(:))/6;
%STACK{4} = mdls.nonlinearity.mdl(struct('phi', [phi1 phi2 phi3 0], ...
%                                        'nlfn', @sigmoidal));
%STACK{3}.fit_fields = {};
%STACK{4}.fit_fields = {'phi'};
%fit_with_lsqcurvefit();
%phi_init=STACK{4}.phi;

STACK{5} = mdls.nonparm_nonlinearity.auto_init(STACK,XXX);
STACK{6} = mdls.correlation;

filename = sprintf('%s/%s_all.mat', savepath, cellid);
save_model_stack(filename, STACK, XXX);

% Now train the model again on each respfile individually, holding the
% depression and FIR coefficients free and letting the nonlinearity slide
nlbins=200;
nlmatrix=zeros(nlbins,filecount-1);
pp=zeros(STACK{5}.bincount,filecount-1);
rr=zeros(STACK{5}.bincount,filecount-1);
rre=zeros(STACK{5}.bincount,filecount-1);
for rfi = 2:filecount,
    rf = files_to_analyze{ci}{rfi};
    XXX = {};
    XXX{1} = [];
    XXX{1}.training_set = {rf};
    XXX{1}.cellid = cellid;  
    XXX{1}.test_set = {};
    recalc_xxx(1);
    STACK{5} = mdls.nonparm_nonlinearity.auto_init(STACK(1:3),XXX(1:4));
    
    %STACK{4}.phi = phi_init;  % Reset phi each time to avoid minima
    %fit_with_lsqcurvefit();
    
    filename = sprintf('%s/%s_%s.mat', savepath, cellid, rf);
    save_model_stack(filename, STACK, XXX);
    
    % save some results for plotting
    x=linspace(linpredmin,linpredmax,nlbins)';
    nlmatrix(:,rfi-1)=raw_nl(STACK{5}.phi,x);
    pp(:,rfi-1)=STACK{5}.phi{1};
    rr(:,rfi-1)=STACK{5}.phi{2};
    rre(:,rfi-1)=STACK{5}.outbinserr;
end

filename = sprintf('%s/%s_all.mat', savepath, cellid);
load_model_stack(filename);
open_narf_gui();
set(gcf,'name',[cellid ' : ' basename(filename)]);

tsp_psth(cellid);

subplot(6,2,7);
h=STACK{3}.coefs;
mm=max(abs(h(:)));
imagesc(h,[-mm mm]);
axis image xy;
title('Linear filter');

ha=subplot(3,2,5);
tsp_gain_comp(cellid,ha);

subplot(3,2,6);
colorset=[0 0 1; 1 0 0; 0.2 0.2 0.8; 0 1 0; 0.4 0.4 0.6];
ht=zeros(filecount-1,1);
lstr={};
for ff=1:filecount-1,
    ht(ff)=errorbar(pp(:,ff),rr(:,ff),rre(:,ff),'Color',colorset(ff,:));
    hold on
    lstr{ff}=respfiles{ff}(8:(end-4));
end
hl=legend(ht,lstr,'Location','Northwest','FontSize',6);
set(hl,'interpreter','none');
set(gcf,'name',[cellid ' : gain per file']);
xlabel('Avg STRF output');
ylabel('File-specific response');

fullpage portrait

fprintf('print command:  print -f2 -djpeg -r150 /auto/users/svd/data/batch241/%s.jpg\n',cellid);

%keyboard