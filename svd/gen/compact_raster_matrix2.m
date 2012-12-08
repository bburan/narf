function[compact_mat] = compact_raster_matrix2(raw_mat, bins_per_frame)
%
%  Required inputs:
%       raw_mat         : A matrix of responses (WITHOUT the psth as the leading row).
%                         time must run along rows while trials runs along columns
%
%       num_segments    : The number of 5 second movie segments
% 
%       bins_per_frame  : The number of time bins per 72.283 Hz movie frame
%
%  Outputs:
%       compact_mat     :  A matrix of responses. Time runs along rows while trials runs along columns.
%                          All -1's have been removed from this data.
%                          The first 13 movie frames have been clipped.
% 
%
% Program Outline:
%       
% 1)      Go through segments and compactify, starting at one
% 2)      Trim header
% 3)      Trim excess trials
% 4)      Truncate trailing -1's
%
%  NB: There is a  discrepancy (two ~14-ms frames) at the end of the first segment/beginning of the second in
%  the psth calculated from the compact_raster_matrix and the psth that is output from d2psth code (or JG's
%  makepsth /psth code).
%
%       
%  This arises because there is a two frame overlap in the data from the two movie segments (i.e. every
%  movie trial contains these two frame). The compact_raster_matrix throws away the leading two frames from
%  the second movie segments. The above mentioned psth programs take these frames into account.

%
% Notes on compacting movies: For multiple segments:
% Each segment appears to be ~4994.3 ms long.
%  
% They are presented in an overlapping fashion: The overlap appears to be 207.5 ms 
% i.e.:
%  X
%  X
%  X   Y
%      Y
%      Y   Z
%          Z
%          Z
%
% Compact raster returns:
%  
%  X
%  X
%  X
%  Y
%  Y
%  Z
%  Z
%
%  My program uses a robust and 'non-analytic' solution to handling overlaps.
%  (The reason for this is the occasional appearance of trials that differ in length by 1 bin under
%  various circumstances of data binning and data set. The origin of these deviations is not 
%  understood).
%
%  The algorithm works in the following manner: 
%       Start at segment #1:
%       Look at first column.
%       Is there data?
%       If no, go to next column.
%       If yes, find the end of the data segment.
%       Extract and store the data.
%       Is the end of the trial data shorter than the previous smallest amount of data for this
%       trial? If so then move the prev_seg_stop to the new shorter length.
%       On the next segment the data start point will reflect this fact (the 'ragged' edge of the
%       previous segment will be over-written but this is ok).
%
%  Graphical example:
% raw_mat:
%  X    X
%  X    X
%  X  Y X Y 
%  X  Y   Y 
%     Y   Y
%     Y   Y
%
%
% output:
%
%  X
%  X 
%  X
%  Y
%  Y
%  y
%   

% Step 0: Setup
% HARD CODE OF _ROUGHLY_HOW MANY BINS THERE ARE IN A 5-sec SEGMENT:
% The 'non-analytic' algorithm does not require an exact knowledge of how many bins there
% are...
% NB: This depends upon binning resolution
% Given in units of rows
seg_length = 361 * bins_per_frame;



% switch statement made obsolete by robust algorithm
%switch bins_per_frame,
%       
%        case 1,
%                seg_length = 361;
%
%        case 2,
%                seg_length = 721;       % 360.5 frames
%
%        case 4,
%                seg_length = 1441;      % 360.25 frames
%
%        case 7,
%                seg_length = 2521;      % 360.1429 frames
%
%        case 8,
%                seg_length = 2880;      % 360.125 frames
%
%        case 14,
%                seg_length = 5040;      % 360 frames
%        
%        otherwise,
%
%                disp('Your bins-per-frame is not supported !');
%                disp('You can inspect the code, guess a seg_length and manually set it.');
%                
%
%end


% HARD CODE FOR LENGTH TO TRIM OFF FRONT OF DATA:
% Given in units of rows
head_trim = ceil(13 * bins_per_frame);

num_rows = size(raw_mat, 1);
num_cols = size(raw_mat, 2);



% Figure out the number of segments
num_segments = round( num_rows / seg_length );



% Step 1: Compactify segments
compact_mat = zeros(size(raw_mat));
min_trials = inf;
prev_seg_stop = 0;

for i = 1:num_segments,

       % Draw a sample point near the middle of the segment. Due to the overlap in presentations this
       % scheme may fail if there are a gigantic number of segments. I don't expect this to happen
           
       sample_point = (seg_length * (i-1)) + round( (seg_length / 2) );

       % Make sure sample point is integer:
       sample_point = round(sample_point);


       if (sample_point < prev_seg_stop),
                disp('PANIC: Pathology in compact_raster_matrix.m !');
                disp('This condition should only be reached if there are an absurd # of data segments !');
                keyboard
       end

       trial_index = 0; 
       lowest_seg_stop = inf;

       for j = 1:num_cols,

                seg_start = prev_seg_stop + 1;

                % Take sample from 'middle' of segment to decide whether the trial was run
                sample = raw_mat(sample_point, j);

                if (sample ~= -1),
                        % Then there is data in this here trial
                        
                        seg_stop = find_seg_stop(raw_mat, sample_point, j);

                        trial_index = trial_index + 1;
                        compact_mat(seg_start:seg_stop,trial_index) = raw_mat(seg_start:seg_stop,j);

                        if (seg_stop < lowest_seg_stop),
                                lowest_seg_stop = seg_stop;
                        end        % if test for lowest_seg_stop 
                end                % if test for data in segment          

        end                        % j loop over columns

        % Set prev_seg_stop
        prev_seg_stop = lowest_seg_stop;
        if (prev_seg_stop == inf),
                disp('PANIC: Pathology in compact_raster_matrix.m !');
                disp('This condition should only be reached if there is a trial which was never presented !');
                keyboard
        end
        
        % check whether current trials are the fewest trials for any segment
        if (trial_index < min_trials),
                min_trials = trial_index;
        end
        

end             % i loop over segments


% Finally trim away any part of the matrix that wasn't filled (i.e. any part that has trailing -1's

compact_mat = compact_mat(1:prev_seg_stop,:);


% Step 2: Trim header I.E. FIRST 13 FRAMES
num_rows = size(compact_mat, 1);

compact_mat = compact_mat( (head_trim+1):num_rows, :);



 
% Step 3: Trim excess trials to make compact

compact_mat = compact_mat(:,1:min_trials);


% Step 4: Truncate off any trailing -1's that happen to exist
% I suspect that the new algorithm makes this step irrelevant.

cut_point = size(compact_mat, 1);

for i = 1:min_trials,
        
        neg_point = min( find( compact_mat(:,i) == -1 ));
        
        if (neg_point < cut_point),
                cut_point = neg_point - 1;
        end

end             % i loop
       

compact_mat = compact_mat(1:cut_point,:);

        







function[seg_stop] = find_seg_stop(raw_mat, sample_point, j)
%
%
% This is hyper-conservative... I start at the sample point and plod through looking
% for the first -1.
% 

num_rows = size(raw_mat, 1);
stop_flag = 0;       
counter = 1;
 
while (stop_flag == 0),

        test_point = sample_point + counter;

        test = raw_mat(test_point, j);

        if (test_point == num_rows),            % Stop if you are at the end of file
                stop_flag = 1;
                seg_stop = test_point;
        end

        if (test == -1),                        % Also stop if you have reached a -1
                stop_flag = 1;
                seg_stop = test_point - 1;      % Back up one bin since you have gone too far        
        end

        counter = counter + 1;

end     % while loop