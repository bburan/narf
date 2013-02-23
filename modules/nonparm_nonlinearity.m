function m = nonparm_nonlinearity(args)
% Applies a nonparametric static nonlinear function to the input. 
%
%
% You may of course define your own functions. 

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @nonparm_nonlinearity;
m.name = 'nonparm_nonlinearity';
m.fn = @do_np_nonlinearity;
m.pretty_name = 'Nonparametric Nonlinearity';
m.editable_fields = {'input_stim', 'input_resp', 'time', 'output', 'bincount'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.get_parms = @init_nonparm_nonlinearity;
m.input_stim = 'stim';
m.input_resp = 'respavg';
m.time = 'stim_time';
m.output = 'stim';
m.bincount = 20;
m.phi={zeros(m.bincount,1),zeros(m.bincount,1)};
m.outbinserr=zeros(m.bincount,1);

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_smooth_scatter_and_nonlinearity; 
m.plot_fns{1}.pretty_name = 'Stim/Resp Smooth Scatter';

m.plot_fns{2}.fn = @(stack, xxx) do_plot_signal(stack, xxx, stack{end}.time, stack{end}.output);
m.plot_fns{2}.pretty_name = 'Output vs Time';

%m.plot_fns{3}.fn = @(stack, xxx) do_plot_nonlinearity(stack, xxx, stack{end}.input_stim, @(x) stack{end}.nlfn(stack{end}.phi, x), false);
%m.plot_fns{3}.pretty_name = 'Nonlinearity';

%m.plot_fns{4}.fn = @(stack, xxx) do_plot_nonlinearity(stack, xxx, stack{end}.input_stim, @(x) stack{end}.nlfn(stack{end}.phi, x), true);
%m.plot_fns{4}.pretty_name = 'Nonlinearity + Histogram';

m.plot_fns{3}.fn = @do_plot_scatter_and_nonlinearity; 
m.plot_fns{3}.pretty_name = 'Stim/Resp Scatter';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

function [phi,outbinserr] = init_nonparm_nonlinearity(stack, xxx)
    % calculate nonparm_nonlinearity parameters
    
    mdl = stack{end};
    x = xxx{end};
    
    % find out parameters:
    pred=[]; resp=[];
    for sf = x.training_set,
        sf=sf{1};
        pred=cat(1,pred,x.dat.(sf).(mdl.input_stim)(:));
        resp=cat(1,resp,x.dat.(sf).(mdl.input_resp)(:));
    end
    
    keepidx=find(~isnan(resp));
    pred=pred(keepidx);
    resp=resp(keepidx);
    
    bincount=mdl.bincount;
    pp=zeros(bincount,1);
    rr=zeros(bincount,1);
    rre=zeros(bincount,1);
    
    if nansum(resp)>0,
        [ss,si1]=sort(pred);
        tb=bincount;
        b=[];
        
        while length(b)<bincount && tb<length(si1),
            tb=tb+1;
            edges1=round(linspace(1,length(si1)+1,tb));
            [b,ui,uj]=unique(ss(edges1(1:(end-1)))');
        end
        edgeend=edges1(end);
        edges1=edges1(ui);
        
        if length(b)<bincount,
            b=[b repmat(b(end),[1 bincount-length(b)])];
            edges1=[edges1 repmat(edges1(end),[1 bincount-length(edges1)])];
        end
        edges1=[edges1 edgeend];
        
        for bb=1:bincount,
            pp(bb)=mean(pred(si1(edges1(bb):(edges1(bb+1)-1))));
            rr(bb)=mean(resp(si1(edges1(bb):(edges1(bb+1)-1))));
            nn=sqrt(edges1(bb+1)-edges1(bb));
            if edges1(bb+1)>edges1(bb),
                rre(bb)=std(resp(si1(edges1(bb):(edges1(bb+1)-1))))./...
                    (nn+(nn==1));
            end
        end
        pp(isnan(pp))=0;
        rr(isnan(rr))=0;
        
        rr(:)=gsmooth(rr(:),1);
        % [pp(:,dd),rr(:,dd)]
    end
    
    phi={pp,rr};
    outbinserr=rre;
    
 end

function x = do_np_nonlinearity(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    [phi,outbinserr] = init_nonparm_nonlinearity(stack, xxx);
    
    for sf = fieldnames(x.dat)', sf=sf{1};
        [T, S, C] = size(x.dat.(sf).(mdl.input_stim));
        y = zeros(T, S, C);
        
        % TODO: If a scalar-valued function, use this
        % y = arrayfun(@(in) mdl.nlfn(mdl.phi, in), x.dat.(sf).(mdl.input_stim));
        
        % Otherwise use the much faster vector valued functions
        y = raw_nl(phi, x.dat.(sf).(mdl.input_stim)(:));      
        
        % TODO: Find a better solution than this hacky way of zeroing nans
        % so that optimization continue in the presence of singularities
        y(isnan(y)) = 0;
        
        x.dat.(sf).(mdl.output) = reshape(y,[T,S,C]);
    end
end

function do_plot_scatter_and_nonlinearity(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    hold on;
    do_plot_scatter(stack, xxx(1:end-1), mdl.input_stim, mdl.input_resp);
    xlims = xlim();
    % xs = linspace(xlims(1), xlims(2), 100);
    % plot(xs, raw_nl(mdl.phi, xs));
    [phi,outbinserr] = init_nonparm_nonlinearity(stack, xxx(1:end-1));
    errorbar(phi{1},phi{2},outbinserr);
    hold off
end

function do_plot_smooth_scatter_and_nonlinearity(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    hold on;
    do_plot_avg_scatter(stack, xxx(1:end-1), mdl.input_stim, mdl.input_resp);
    % xlims = xlim();
    % xs = linspace(xlims(1), xlims(2), 100)';
    % plot(xs, raw_nl(mdl.phi, xs));
    [phi,outbinserr] = init_nonparm_nonlinearity(stack, xxx(1:end-1));
    errorbar(phi{1},phi{2},outbinserr);
    hold off
end

end