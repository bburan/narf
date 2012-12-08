function r = gzipped(fname)
%function r = gzipped(fname)
%
% is a file gzipped?
%
%Fri Apr 13 16:46:02 2001 mazer 

r = strcmp(fname((end-2):end), '.gz');
