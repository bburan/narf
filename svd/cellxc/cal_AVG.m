function [stim_avg, Avg_psth, psth, errFlg] = cal_AVG(DS, nband, lin_flag, sil_window)
%
% [stim_avg, Avg_psth, psth] = cal_AVG(DS, lin_flag, sil_window)
%     -- Calculate the average stimulus value across time.
%     -- Calculate average of psth over all trials 
%     -- Calculate the psth of response file over all trials.
%   Input:
%      DS(required input):  the cell of each data struct 
%                               that contains four fields:
%               stimfiles  - stimulus file name
%               respfiles  - response file name
%               nlen    - length of time domain
%               ntrials    - num of trials
%              e.g. DS{1} = struct('stimfiles', 'stim1.dat', 'respfiles', 
%                    'resp1.dat', 'nlen', 1723, 'ntrials', 20);
%      lin_flag: the flag to show whether we need take log on data
%              e.g. lin_flag = 0 if we need take log, 1(default) if otherwise 
%      sil_window: the interval for not counting when preprocessing data
%              e.g. sil_window = 0(default) 
%   Output:
%      stim_avg:  average stimulus value which is only function of space.
%      Avg_psth: average psth which is the number
%      psth: the cell of psth for one data pair which is only func of time 
%      
%             STRFPAK: STRF Estimation Software
% Copyright ©2003. The Regents of the University of California (Regents).
% All Rights Reserved.
% Created by Theunissen Lab and Gallant Lab, Department of Psychology, Un
% -iversity of California, Berkeley.
%
% Permission to use, copy, and modify this software and its documentation
% for educational, research, and not-for-profit purposes, without fee and
% without a signed licensing agreement, is hereby granted, provided that
% the above copyright notice, this paragraph and the following two paragr
% -aphs appear in all copies and modifications. Contact The Office of Tec
% -hnology Licensing, UC Berkeley, 2150 Shattuck Avenue, Suite 510,
% Berkeley, CA 94720-1620, (510) 643-7201, for commercial licensing
% opportunities.
%
%IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
%SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
%ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
%REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
%REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT
%LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
%PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY,
%PROVIDED HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PRO
%-VIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

% created by JXZ
% Sep. 25, 2002
%
%


% ========================================================
% check whether we have valid required input
% ========================================================
errFlg = 0;
if isempty(DS) 
    errordlg('ERROR: Please enter non-empty data filename',...
        'Input Data Error', 'modal')
    errFlg = 1;
    return
end

if ~exist('nband')
    global NBAND
    if isempty(NBAND)
        errordlg('You need assign global variable NBAND first.',...
                 'Global Variable Error', 'modal');
        errFlg = 1;    
        return
    end
    nband = NBAND;
end

% ========================================================
% if user does not provide other input, we set default.
% ========================================================
% lin_flag refers to linear data if lin_flag = 1.
% otherwise, we take logrithm on data if lin_flag = 0.
if ~exist('lin_flag')
    lin_flag = 1;
end

if ~exist('sil_window')
    sil_window = 0;
end

% ========================================================
% initialize output and declare local variables 
% ========================================================
% find out the total number of data files
ndata_files = length(DS);

stim_avg = zeros(nband, 1);
count_avg = 0;
tot_trials = 0;
psth = {};
Avg_psth = 0;

% ========================================================
% calculate the output over all the data file 
% ========================================================
for n = 1:ndata_files
   % load stimulus files
   stim_env = Check_And_Load(DS{n}.stimfiles);

   % load response files
   psth_rec = Check_And_Load(DS{n}.respfiles);
   
   % take logrithm of data based on lin_flag 
   if lin_flag == 0
	stim_env = log(stim_env + 1.0);
   end

   % calculate stim avg
   %
   % Before Do calculation, we want to check if we got the correct input
   tempXsize = size(stim_env,1);
   if tempXsize ~= nband
       errordlg(['Data Error: Please check your input data by clicking ',...
     '"Get Files" Button in the main window: The first data file need ',... 
     'to be stimuli and the second data file need to be its corresponding',... 
     ' response file. If you made a mistake, please type "clear all" ',...
     ' or hit "reset" button first and then choose input data again.'],...
     'Input Data Error', 'modal');
        errFlg = 1;
        return
    end
    
   stim_avg = stim_avg + sum(stim_env*DS{n}.ntrials, 2);
   count_avg = count_avg +(DS{n}.nlen + 2*sil_window)*DS{n}.ntrials;  
    
    % then calculate response_avg
   if DS{n}.ntrials > 1
    psth{n} = mean(psth_rec);
   else
    psth{n} = psth_rec;
   end

    tot_trials = tot_trials + DS{n}.nlen + sil_window;
    
    % calculate the total spike/response avg.
    Avg_psth = Avg_psth + sum(psth{n});
    
    % clear workspace
    clear stim_env
    clear psth_rec
end

% ========================================================
% save the stim_avg into the data file
% ========================================================
currentPath = pwd;
global outputPath
if ~isempty(outputPath)
    cd (outputPath);
else
    disp('Saving output to Output Dir.');
    stat = mkdir('Output');
    cd('Output');
    outputPath = pwd;
end

stim_avg = stim_avg/count_avg;
Avg_psth = Avg_psth / tot_trials;

save('stim_avg.mat', 'stim_avg');

cd(currentPath);
% ========================================================
% END OF CAL_AVG
% ========================================================
