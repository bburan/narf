function load_model(filepath)
% load_model(filepath)
%
% Loads the NARF model found at path FILEPATH into the global variables 
% STACK, XXX, and META. Their contents are overwritten if they exist.
%
% ARGUMENTS
%    filepath    The absolute path to where the NARF .mat file is found
% 
% RETURNS: Nothing

global STACK XXX META;
vars = load(filepath, 'stack', 'xxx', 'meta');

XXX = {};

% Swap in the new globals
XXX{1} = vars.xxx;
META = vars.meta;
STACK = vars.stack;  
