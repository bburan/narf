function narf_set_path()
    global NARF_PATH NARF_MODULES_PATH NARF_SAVED_MODELS_PATH;
    NARF_PATH = fileparts(which('narf_set_path')); 
    NARF_MODULES_PATH = [NARF_PATH filesep 'modules'];
    NARF_SAVED_MODELS_PATH = [NARF_PATH filesep 'saved_models'];
    addpath(NARF_PATH, ...
            NARF_MODULES_PATH, ...
            [NARF_PATH filesep 'analysis'], ...
            [NARF_PATH filesep 'optimization'], ...
            [NARF_PATH filesep 'svd'], ...    
            [NARF_PATH filesep 'svd/autils'], ...
            [NARF_PATH filesep 'svd/cellxc'], ...
            [NARF_PATH filesep 'svd/db'], ...
            [NARF_PATH filesep 'svd/gen'], ...
            [NARF_PATH filesep 'utils']);
end