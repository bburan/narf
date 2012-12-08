%to use call:
%
%cleaned = CleanData(rawdata);     % to clean the data without plotting
%cleaned = CleanData(rawdata,0);   % to clean the data without plotting
%cleaned = CleanData(rawdata,1);   % to clean the data with plotting
%
%rawdata is the original data of channels and time, in a 2-D matlab array.
%cleaned is a 2-D array, [time,channels+2] of cleaned data. The last
%two columns returned are the first and second principal component vectors.
%
%CleanData.m will clean the data array using a principal
%component analysis to find common signal across channels
%and remove it.  Edit the initalized value of the global variables in this
%file appropriately to optimize cleaning of your data.

%***>NOTE: Low sampling rate relative to data bandwidth and multiplexing 
%   offsets between channels may limit performance. See:
%   "Signal to noise ratio improvement in multiple electrode recording", 
%   Musial, Baker, Gerstein, King and Keating, 
%   Journal of Neuroscience Methods 115: 29-43, 2002.  		<****

%Questions?  jeff@mulab.physiol.upenn.edu
%	   george@mulab.physiol.upenn.edu
