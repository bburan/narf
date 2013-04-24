function m = nonparm_nonlinearity(args)
% Applies a nonparametric static nonlinear function to the input. 

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

% Optional fields
m.plot_fns = {};
m.auto_plot = @do_plot_smooth_scatter_and_nonlinearity;
m.plot_fns{1}.fn = @do_plot_smooth_scatter_and_nonlinearity; 
m.plot_fns{1}.pretty_name = 'Stim/Resp Smooth Scatter';
m.plot_fns{2}.fn = @do_plot_scatter_and_nonlinearity; 
m.plot_fns{2}.pretty_name = 'Stim/Resp Scatter';
m.plot_fns{3}.fn = @do_plot_all_default_outputs;
m.plot_fns{3}.pretty_name = 'Output Channels (All)';
m.plot_fns{4}.fn = @do_plot_single_default_output;
m.plot_fns{4}.pretty_name = 'Output Channel (Single)';


% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% ------------------------------------------------------------------------
% Methods

function [phi,outbinserr] = init_nonparm_nonlinearity(mdl, x)
    % calculate nonparm_nonlinearity parameters
    
        
    % find out parameters:
    pred=[]; resp=[];
    for sf = x.training_set,
        sf=sf{1};
        % It's an error if there is more than one channel:
        if size(x.dat.(sf).(mdl.input_stim)(:), 3) > 1 || ...
           size(x.dat.(sf).(mdl.input_resp)(:), 3) > 1,
            error('NPNL is a 1D method!');
        end
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

function x = do_np_nonlinearity(mdl, x, stack, xxx)    
    [phi, ~] = init_nonparm_nonlinearity(mdl, x);
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
        sf = fns{ii};
        [T, S, C] = size(x.dat.(sf).(mdl.input_stim));
        y = raw_nl(phi, x.dat.(sf).(mdl.input_stim)(:));                      
        y(isnan(y)) = 0; % TODO: Find better solution than zeroing to avoid singularities        
        x.dat.(sf).(mdl.output) = reshape(y,[T,S,C]);
    end
end

function help_plot_npnl(sel, mdls, xins, xouts)
    for ii = 1:length(mdls)
        [phi,outbinserr] = init_nonparm_nonlinearity(mdls{ii}, xins{ii}{end});
        xouts{ii}.dat.(sel.stimfile).npnlstim = phi{1};
        xouts{ii}.dat.(sel.stimfile).npnlpred = phi{2};
        xouts{ii}.dat.(sel.stimfile).temperr  = outbinserr;
    end
    
    hold on;    
    do_plot(xouts, 'npnlstim', 'npnlpred', ...
            sel, 'NPNL Input [-]', 'RespAvg Prediction [Hz]');   
	hold off;
end

function do_plot_scatter_and_nonlinearity(sel, stack, xxx)    
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1));    
    do_plot_scatter(sel, xins, mdls{1}.input_stim, mdls{1}.input_resp);  
    help_plot_npnl(sel, mdls, xins, xouts);
end

function do_plot_smooth_scatter_and_nonlinearity(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    do_plot_scatter(sel, xins, mdls{1}.input_stim, mdls{1}.input_resp, 100); 
    help_plot_npnl(sel, mdls, xins, xouts);
end

end