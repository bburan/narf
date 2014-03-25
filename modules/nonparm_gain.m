function m = nonparm_gain(args)
% Applies a 2D nonparametric static nonlinear function to the inputs

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @nonparm_gain;
m.name = 'nonparm_gain';
m.fn = @do_nonparm_gain;
m.pretty_name = 'Nonparametric Gain Scaling';
m.editable_fields = {'input_stim1', 'input_stim2', 'input_resp', 'time',...
                     'output', 'bincount', 'smoothwindow'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input_stim1 = 'stim';   % base input
m.input_stim2 = 'stim2';   % input to gain scale term
m.input_resp = 'respavg';
m.time = 'stim_time';
m.output = 'stim';
m.bincount = 20;
m.smoothwindow=1;

% Optional fields
m.is_splittable = true;
m.plot_fns = {};
m.auto_plot = @do_plot_smooth_scatter_and_nonlinearity;
m.plot_fns{1}.fn = @do_plot_scatter_and_nonlinearity;
m.plot_fns{1}.pretty_name = 'Stim/Resp Scatter';
m.plot_fns{2}.fn = @do_plot_smooth_scatter_and_nonlinearity;
m.plot_fns{2}.pretty_name = 'Stim/Resp Smooth Scatter';
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

function [phi,outbinserr] = init_nonparm_gain(mdl, x)
    
    % find out parameters:
    pred1=[]; pred2=[]; resp=[];
    for sf = x.training_set,
        sf=sf{1};
        % It's an error if there is more than one channel:
        if size(x.dat.(sf).(mdl.input_stim1), 3) > 1 || ...
           size(x.dat.(sf).(mdl.input_stim2), 3) > 1 || ...
           size(x.dat.(sf).(mdl.input_resp), 3) > 1,
            error('NPNL gain requires two single channels.');
        end
        pred1=cat(1,pred1,x.dat.(sf).(mdl.input_stim1)(:));
        pred2=cat(1,pred2,x.dat.(sf).(mdl.input_stim2)(:));
        resp=cat(1,resp,x.dat.(sf).(mdl.input_resp)(:));
    end
    
    keepidx=find(~isnan(resp));
    pred1=pred1(keepidx);  % this is the output of the STRF
    pred2=pred2(keepidx);  % this is the term that will feed the gain
    resp=resp(keepidx);
    
    bincount=mdl.bincount;
    smoothwindow=mdl.smoothwindow;
    pp=zeros(bincount,1);
    rr=zeros(bincount,1);
    rre=zeros(bincount,1);
    
    if nansum(resp),
        [ss2,si2]=sort(pred2);
        tb=bincount;
        b2=[];
        
        while length(b2)<bincount && tb<length(si2),
            tb=tb+1;
            edges2=round(linspace(1,length(si2)+1,tb));
            if std(ss2(edges2(1:(end-1))))>0
                [b2,ui]=unique(ss2(edges2(1:(end-1)))');
            else
                b2=ss1(edges2(1:(end-1)))';
                ui=1:length(b2);
            end
        end
        edgeend=edges2(end);
        edges2=edges2(ui);
        
        if length(b2)<bincount,
            b2=[b2 repmat(b2(end),[1 bincount-length(b2)])];
            edges2=[edges2 repmat(edges2(end),[1 bincount-length(edges2)])];
        end
        edges2=[edges2 edgeend];
        b2(bincount+1)=ss2(end).*1.1;
        
        for bb=1:bincount,
            pp(bb,1)=mean(ss2(edges2(bb):(edges2(bb+1)-1)));
        end
        
        for bb=1:bincount,
           ff=find(pred2>=b2(bb) & pred2<b2(bb+1));
           rr(bb)=mean(resp(ff))./mean(pred1(ff));
           rre(bb)=length(ff);
        end
        pp(isnan(pp))=0;
        rr(isnan(rr))=0;
        if smoothwindow
           rr=gsmooth(rr,[smoothwindow]);
        end
        
    end
    
    phi={pp,rr};
    outbinserr=rre;
end

function x = do_nonparm_gain(mdl, x, stack, xxx)    
    
    [phi, ~] = init_nonparm_gain(mdl, x);
    
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
        sf = fns{ii};
        [T, S, C] = size(x.dat.(sf).(mdl.input_stim1));
        s1=x.dat.(sf).(mdl.input_stim1)(:);
        s2=x.dat.(sf).(mdl.input_stim2)(:);
        g = raw_nl(phi, s2);                  
        g(isnan(g)) = 0;
        x.dat.(sf).(mdl.output) = reshape(g.*s1,[T,S,C]);
    end
    
    % I realized it's easier to cache the nonlinearity for later plotting
    % FIXME: storing the coefs here won't work with existing splitters/unifiers
    x.mycoefmatrix = phi; 
    
end

function help_plot_npnl(sel, mdls, xins, xouts)
    for ii = 1:length(mdls)
        [phi,outbinserr] = init_nonparm_gain(mdls{ii}, xins{ii}{end});
        xouts{ii}.dat.(sel.stimfile).npnlstim = phi{1};
        xouts{ii}.dat.(sel.stimfile).npnlpred = phi{2};
        xouts{ii}.dat.(sel.stimfile).temperr  = outbinserr;
    end
    
    hold on;    
    tsel=sel;
    tsel.stim_idx=1;
    do_plot(xouts, 'npnlstim', 'npnlpred', ...
            tsel, 'NPNL Input [-]', 'RespAvg Prediction [Hz]');   
	hold off;
end

function do_plot_scatter_and_nonlinearity(sel, stack, xxx)    
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1));    
    do_plot_scatter(sel, xins, mdls{1}.input_stim2, mdls{1}.input_resp);  
    help_plot_npnl(sel, mdls, xins, xouts);
end

function do_plot_smooth_scatter_and_nonlinearity(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    %do_plot_scatter(sel, xins, mdls{1}.input_stim2, mdls{1}.input_resp, 100); 
    help_plot_npnl(sel, mdls, xins, xouts);
end


end