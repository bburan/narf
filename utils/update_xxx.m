function update_xxx(start_depth, end_depth)
% update_xxx(start_depth, end_depth)
%
% Updates the XXX data structure, overwriting any previous data as needed.
%
% A very important function, update_xxx is responsible for actually 
% computing the prediction made by a model.
%
% Note that START_DEPTH is not actually just the STACK index, but actually
% the "depth" of the module, it's position relative to the start of the
% STACK structure. That is, if you had STACK{1}{1} = A, STACK{2}{1} = B, 
% STACK{2}{2} = C, and STACK{3,1} = D, then the "Depth" of D is 4. 
%
% ARGUMENTS:
%    start_depth    The depth of the module to begin recomputation from.
%    end_depth      Defaults to the last possible module.
% 
% RETURNS: Nothing. But it usually has a lot of side effects!

global XXX STACK FLATXXX;

if start_depth < 1
    error('Start depth less than 1 not allowed.');
end

% Make a flattened list of modules from STACK
flatstack = {};
xxxindexes = [];
for ii = 1:length(STACK)
    for jj = 1:length(STACK{ii})
        flatstack{end+1} = STACK{ii}{jj};      
        xxxindexes(end+1) = ii; 
    end
end

if ~exist('end_depth', 'var')
    end_depth = length(flatstack);
end

% If XXX is not fully initialized, we need to compute earlier
if length(XXX)==1,
   FLATXXX=XXX;
end
if xxxindexes(start_depth) > length(XXX)
    start_depth = find(xxxindexes, length(XXX), 'first');
    FLATXXX=FLATXXX(1:start_depth);
end
if length(FLATXXX)<start_depth,
   if isempty(FLATXXX),
      FLATXXX=XXX(1);
      start_depth=1;
   else
      start_depth=length(FLATXXX)-1;
   end
end

%%START SVD REMOVED 2014-06-14
% Remove references to data 
%if xxxindexes(start_depth) < length(XXX)
%    XXX = XXX(1:xxxindexes(start_depth)); 
%end

%fprintf('startdepth is: %d\n', start_depth);

% Build up a flattened XXX structure
%flatxxx = {};  % What goes IN to each module
%for ii = 1:start_depth
%    flatxxx{ii} = XXX{xxxindexes(ii)};     
%end
%%END SVD REMOVED 2014-06-14

% Work through the flattened stack, skipping computations that can be
% skipped whenever
modified = {}; % A list of modified signals
for ii = start_depth:end_depth,    
    % If the affected list is empty, we must compute because this is the first
    % module that we are updating. After that, we may skip the module if it
    % has requirements that need not be updated.    
    
    if ~isfield(flatstack{ii}, 'required') || ~isfield(flatstack{ii}, 'modifies') 
        fprintf('WARNING: NARF now requires that modules declare their dependencies (m.required) and which signals they alter (m.modifies).\n'); 
        fprintf('If you are seeing this in any other circumstance than when loading an old model, you probably have subtle errors regarding conditional evaluation!\n');
%%START SVD REMOVED 2014-06-14
        %flatxxx{ii+1} = flatstack{ii}.fn(flatstack{ii}, flatxxx{ii});        
%%END SVD REMOVED 2014-06-14
        FLATXXX{ii+1} = flatstack{ii}.fn(flatstack{ii}, FLATXXX{ii});        
    else    
%%START SVD REMOVED 2014-06-14
%        if ~isempty(modified) && ...
%                isempty(intersect(flatstack{ii}.required, modified)) && ...
%                length(XXX) > (xxxindexes(ii) + 1)
            % No computation required. Pass through existing signals.  
            %flatxxx{ii+1} = merge_structs(XXX{xxxindexes(ii)+1}, flatxxx{ii});
%             flatxxx{ii+1}=XXX{xxxindexes(ii)+1};
%             fns = fieldnames(flatxxx{ii+1}.dat);
%             for idx = 1:length(fns);
%                 for jj=1:length(modified),
%                 flatxxx{ii+1}.dat.(fns{idx}).(modified{jj})= ...
%                     flatxxx{ii}.dat.(fns{idx}).(modified{jj});
%                 end
%             end
%%END SVD REMOVED 2014-06-14
        if ~isempty(modified) && ...
                isempty(intersect(flatstack{ii}.required, modified)) &&...
                length(FLATXXX)>ii,
            % No computation required. Pass through existing signals.  
            fns = fieldnames(FLATXXX{ii+1}.dat);
            for jj=1:length(modified),
               if isfield(FLATXXX{ii},modified{jj}),
                  FLATXXX{ii+1}.(modified{jj})= ...
                     FLATXXX{ii}.(modified{jj});
               else
                  for idx = 1:length(fns)
                     FLATXXX{ii+1}.dat.(fns{idx}).(modified{jj})= ...
                        FLATXXX{ii}.dat.(fns{idx}).(modified{jj});
                  end
               end
            end
        else
            % Uncomment this to help debug your tree traversals
            %if isfield(flatstack{ii}, 'output')
            %    fprintf('updating[%d]: %s\n', ii, flatstack{ii}.output);
            %    disp(modified);
            %end
%%START SVD REMOVED 2014-06-14
            %flatxxx{ii+1} = flatstack{ii}.fn(flatstack{ii}, flatxxx{ii});
%%END SVD REMOVED 2014-06-14
            FLATXXX{ii+1} = flatstack{ii}.fn(flatstack{ii}, FLATXXX{ii});
            modified = union(flatstack{ii}.modifies , modified);            
        end
    end    
end

% Now update the main structure
for ii = start_depth:end_depth,    
    XXX{xxxindexes(ii)+1} = FLATXXX{ii+1};
end

