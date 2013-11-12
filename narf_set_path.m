function narf_set_path()   
    % Bad side effect of javaaddpath: it clears all variables!
    if isempty(which('TableSorter'))
        javaaddpath('/home/ivar/matlab/narf/libs/TableSorter.jar');
    end
    if isempty(which('TableColumnAdjuster'))
        javaaddpath('/home/ivar/matlab/narf/libs/TableColumnAdjuster.jar');
    end
    
    global NARF_PATH NARF_MODULES_PATH NARF_SAVED_MODELS_PATH ...
           NARF_SAVED_IMAGES_PATH NARF_SCRIPTS_PATH NARF_LIBS_PATH...
           NARF_KEYWORDS_PATH NARF_FITTERS_PATH NARF_MEMOIZATION_PATH;
    NARF_PATH = fileparts(which('narf_set_path'));   
    NARF_FITTERS_PATH = [NARF_PATH filesep 'fitters'];
    NARF_KEYWORDS_PATH = [NARF_PATH filesep 'keywords'];
    NARF_MODULES_PATH = [NARF_PATH filesep 'modules'];   
    NARF_SAVED_MODELS_PATH   = '/auto/data/code/saved_models';
    NARF_SAVED_IMAGES_PATH   = '/auto/data/code/saved_images';
    NARF_MEMOIZATION_PATH    = '/auto/data/code/memoization';
    NARF_SCRIPTS_PATH = [NARF_PATH filesep 'scripts'];
    NARF_LIBS_PATH = [NARF_PATH filesep 'libs'];
    
    warning off MATLAB:dispatcher:nameConflict;
    addpath(NARF_MODULES_PATH, ...
            NARF_FITTERS_PATH, ...
            NARF_LIBS_PATH, ...        
            [NARF_PATH filesep 'methods'], ...
            [NARF_PATH filesep 'queue'], ...    
            [NARF_PATH filesep 'svd'], ...    
            [NARF_PATH filesep 'svd' filesep 'autils'], ...
            [NARF_PATH filesep 'svd' filesep 'cellxc'], ...
            [NARF_PATH filesep 'svd' filesep 'gen'], ...
            [NARF_PATH filesep 'utils']);

    addpath(NARF_PATH);
    
    global BAPHYHOME
    if isempty(BAPHYHOME),
        baphy_set_path
    end
    warning on MATLAB:dispatcher:nameConflict;        
end
