function MF = memoize(F)
% Memoizes function F. Doesn't work recursively.

global NARF_MEMOIZATION_PATH;

function_name = func2str(F);
memoization_dir = [NARF_MEMOIZATION_PATH '/' function_name '/'];
[s,w]=unix(['chmod 777 ' memoization_dir]);

    function varargout = INEEDAGENSYM (varargin) 
        %save('/tmp/ineedagensymdump', 'varargin', '-mat');
        %[s, h] = unix('md5sum /tmp/ineedagensymdump.mat | cut -d " " -f1');
        %h = h(1:end-1); % Remove newline
        h = DataHash(varargin);       % Calc MD5 sum hash of arguments
        thedir = [memoization_dir h(1) '/' h(2) '/'];
        f = [thedir h];
        
        % If matching MD5 hash exists, load result immediately
        if exist(thedir, 'dir') && exist(f, 'file')
            vars = load([thedir h], '-mat');
            % FIXME: Check for a very unlikely hash collision
            % I wanted to do isequal(vars.varargin, varargin), 
            % to check for a hash collision more explicitly, but for some
            % reason isequal fails even when the hashes are equal. 
            % Since I see no way to resolve this, we will just have to
            % accept very infrequent hash collisions. 
            fprintf('Using memoized result [%s]\n', ...
                [thedir h]);
            varargout = vars.varargout;
            return
        end 
        % If not, run function and cache result if results take longer than
        % 5 minutes to run. 
        t_start = tic;
        nouts = nargout(F);
        if (nouts == 1)
            [arg1] = F(varargin{:});
            varargout = {arg1};
        elseif (nouts == 2)
            [arg1 arg2] = F(varargin{:});
            varargout = {arg1 arg2};
        elseif (nouts == 3)
            [arg1 arg2 arg3] = F(varargin{:});
            varargout = {arg1 arg2 arg3};
        elseif (nouts == 4)
            [arg1 arg2 arg3 arg4] = F(varargin{:});
            varargout = {arg1 arg2 arg3 arg4};
        else
            error('Only supports 4 output args right now because MATLAB has crappy syntax, wonky destructuring binds, and I cannot figure out how to use varargout.');
        end
        t_elapsed = toc(t_start);   
        
        if t_elapsed > 300 % Ensure that we only save if it was longer than 5 minutes
            if ~exist(thedir, 'dir')
                mkdir(thedir);
                unix(['chmod 777 ' thedir]);
            end
            save(f, 'varargin', 'varargout', '-mat');
            unix(['chmod 777 ' f]);
        end
    end

MF = @(varargin) INEEDAGENSYM(varargin{:});    

end
