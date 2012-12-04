function narf_set_path()
    NARF_PATH = '/home/ivar/matlab/narf/';
    addpath([NARF_PATH filesep 'utils'], ...
            [NARF_PATH filesep 'modules']);
end