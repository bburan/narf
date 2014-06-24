function split_stack(startidx,endidx)

global STACK META XXX MODULES

if ~exist('startidx','var'),
   startidx=2;
end

if ~exist('endidx','var'),
   error('currently you need to specify endidx');
end

oldstack=STACK;

STACK=STACK(1:(startidx-1));
XXX=XXX(1:2);

append_module(MODULES.split_signal);
split_mod=STACK{end}{1};
split_count=length(split_mod.outputs);
for ii=startidx:endidx,
    mdl=oldstack{ii}{1};
    if isfield(mdl,'input'),
        oldinput=mdl.input;
    end
    if isfield(mdl,'input_stim'),
        oldinput=mdl.input_stim;
    end
    oldinputidx=find(strcmp(mdl.required,oldinput));
    for jj=1:split_count,
        newinput=split_mod.outputs{jj};
        if isfield(mdl,'input'),
            mdl.input=newinput;
        end
        if isfield(mdl,'input_stim'),
            mdl.input_stim=newinput;
        end
        if isfield(mdl,'output'),
            mdl.output=newinput;
        end
        if isfield(mdl,'output_stim'),
            mdl.output_stim=newinput;
        end
        
        mdl.modifies={newinput};
        mdl.required{oldinputidx}=newinput;
        STACK{ii+1}{jj}=mdl;
    end
end
update_xxx(startidx);

append_module(MODULES.unify_signal);
for ii=(endidx+1):length(oldstack),
    STACK{end+1}=oldstack{ii};
end

