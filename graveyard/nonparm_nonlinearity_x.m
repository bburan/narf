function m = nonparm_nonlinearity_x(args)
% Same as nonparametric_nonlinearity, but with per-file nonlinearities. 

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @nonparm_nonlinearity_x;
m.name = 'nonparm_nonlinearity_x';
m.fn = @do_np_nonlinearity_x;
m.pretty_name = 'Nonparametric Nonlinearity';
m.editable_fields = {'input_stim', 'input_resp', 'time', 'output', 'bincount'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.get_parms = @init_nonparm_nonlinearity_x;
m.input_stim = 'stim';
m.input_resp = 'respavg';
m.time = 'stim_time';
m.output = 'stim';
m.bincount = 20;
m.phi={zeros(m.bincount,1),zeros(m.bincount,1)};
m.outbinserr=zeros(m.bincount,1);

% Optional fields
m.plot_fns = {};
m.plot_fns{1}.fn = @do_plot_smooth_scatter_and_nonlinearity_x; 
m.plot_fns{1}.pretty_name = 'Stim/Resp Smooth Scatter';

m.plot_fns{2}.fn = @(stack, xxx) do_plot_signal(stack, xxx, stack{end}.time, stack{end}.output);
m.plot_fns{2}.pretty_name = 'Output vs Time';

% m.plot_fns{3}.fn = @(stack, xxx) do_plot_nonlinearity_x(stack, xxx, stack{end}.input_stim, @(x) stack{end}.nlfn(stack{end}.phi, x), false);
% m.plot_fns{3}.pretty_name = 'Nonlinearity';

% m.plot_fns{4}.fn = @(stack, xxx) do_plot_nonlinearity_x(stack, xxx, stack{end}.input_stim, @(x) stack{end}.nlfn(stack{end}.phi, x), true);
% m.plot_fns{4}.pretty_name = 'Nonlinearity + Histogram';
% 
% m.plot_fns{3}.fn = @do_plot_scatter_and_nonlinearity_x; 
% m.plot_fns{3}.pretty_name = 'Stim/Resp Scatter';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

function [ppp, rrr, outbinserr] = init_nonparm_nonlinearity_x(stack, xxx)
    % calculate nonparm_nonlinearity_x parameters
    
    mdl = stack{end};
    x = xxx{end};
    
    % For each file in the training set
    if length(x.training_set) ~= length(x.test_set)
        error('Training and Test sets must be equal to use NPNLX!');
    end
        
    ppp = {};
    rrr = {};
    outbinserr = {};
    
    % For each training/test file combination, have a unique nonlinearity
    for ii = 1:length(x.training_set),
        sf = x.training_set{ii};

        pred=x.dat.(sf).(mdl.input_stim)(:);
        resp=x.dat.(sf).(mdl.input_resp)(:);
           
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
    
        ppp{ii} = pp;
        rrr{ii} = rr;
        outbinserr{ii}=rre;
    end
    
 end

function x = do_np_nonlinearity_x(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    
    [ppp, rrr, outbinserr] = init_nonparm_nonlinearity_x(stack, xxx);
    
    function ret = use_nth_nl(n, sf)
        [T, S, C] = size(x.dat.(sf).(mdl.input_stim));
        y = zeros(T, S, C);
        y = raw_nl({ppp{ii}, rrr{ii}}, x.dat.(sf).(mdl.input_stim)(:));      
        y(isnan(y)) = 0;
        ret = reshape(y,[T,S,C]);
    end
    
    % Training sets only
    for ii = 1:length(x.training_set),        
        sf = x.training_set{ii};
        x.dat.(sf).(mdl.output) = use_nth_nl(ii, sf);
    end
       
    % Test sets only    
    for ii = 1:length(x.test_set),        
        sf = x.test_set{ii};
        x.dat.(sf).(mdl.output) = use_nth_nl(ii, sf);
    end 
    
end

% function do_plot_scatter_and_nonlinearity_x(stack, xxx)
%     mdl = stack{end};
%     x = xxx{end};
%     
%     hold on;
%     do_plot_scatter(stack, xxx(1:end-1), mdl.input_stim, mdl.input_resp);
%     xlims = xlim();
%     % xs = linspace(xlims(1), xlims(2), 100);
%     % plot(xs, raw_nl(mdl.phi, xs));
%     [phi,outbinserr] = init_nonparm_nonlinearity_x(stack, xxx(1:end-1));
%     errorbar(phi{1},phi{2},outbinserr);
%     hold off
% end

function do_plot_smooth_scatter_and_nonlinearity_x(stack, xxx)
    mdl = stack{end};
    x = xxx{end};
    hold on;
    do_plot_avg_scatter(stack, xxx(1:end-1), mdl.input_stim, mdl.input_resp);
    % xlims = xlim();
    % xs = linspace(xlims(1), xlims(2), 100)';
    % plot(xs, raw_nl(mdl.phi, xs));
    [ppp, rrr, outbinserr] = init_nonparm_nonlinearity_x(stack, xxx(1:end-1));
    
    [sf, ~, ~] = get_baphy_plot_controls(stack);
    for ii = 1:length(x.training_set)
        errorbar(ppp{ii}, rrr{ii}, outbinserr{ii}, pickcolor(ii));
    end
    h = legend(sf, x.training_set{:},'Location','EastOutside');
    set(h,'Interpreter','none');
    hold off
    minrange=cat(1,rrr{:})-cat(1,outbinserr{:})*2;
    maxrange=cat(1,rrr{:})+cat(1,outbinserr{:})*2;
    aa=axis;
    axis([aa(1:2) min(minrange) max(maxrange)]);
end

end