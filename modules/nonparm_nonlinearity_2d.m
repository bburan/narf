function m = nonparm_nonlinearity_2d(args)
% Applies a 2D nonparametric static nonlinear function to the inputs

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @nonparm_nonlinearity_2d;
m.name = 'nonparm_nonlinearity_2d';
m.fn = @do_np_nonlinearity_2d;
m.pretty_name = 'Nonparametric Nonlinearity 2D';
m.editable_fields = {'input_stim1', 'input_stim2', 'input_resp', 'time',...
                     'output', 'bincount'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input_stim1 = 'stim1';
m.input_stim2 = 'stim2';
m.input_resp = 'respavg';
m.time = 'stim_time';
m.output = 'stim';
m.bincount = 20;

% Optional fields
m.plot_fns = {};
m.auto_plot = @do_plot_2d_nonlinearity;
m.plot_fns{1}.fn = @do_plot_2d_nonlinearity; 
m.plot_fns{1}.pretty_name = 'Stim/Resp Smooth Scatter 2D';
m.plot_fns{2}.fn = @do_plot_single_default_output;
m.plot_fns{2}.pretty_name = 'Output Channel (Single)';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% ------------------------------------------------------------------------
% Methods

function [phi,outbinserr] = init_nonparm_nonlinearity_2d(mdl, x, stack, xxx)
    
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

function x = do_np_nonlinearity_2d(mdl, x, stack, xxx)    
    
    [phi2d, ~] = init_nonparm_nonlinearity_2d(mdl, x, stack, xxx);
    
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
        sf = fns{ii};
        [T, S, C] = size(x.dat.(sf).(mdl.input_stim));
        y = raw_nl2d(phi2d, x.dat.(sf).(mdl.input_stim)(:));                      
        y(isnan(y)) = 0; 
        x.dat.(sf).(mdl.output) = reshape(y,[T,S,C]);
    end
    
    % I realized it's easier to cache the nonlinearity for later plotting
    % FIXME: storing the coefs here won't work with existing splitters/unifiers
    x.mycoefmatrix = phi2d; 
    
end

function do_plot_2d_nonlinearity(sel, stack, xxx)    
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    
    % For now, let's NOT make the plot function work for multiple paramsets
    % We'll just assume that we have only a single paramset for the NPNL2
    % FIXME: Later on we can correct this as needed
    mdl = mdls{1};
    xin = xins{1};
    xout = xouts{1};
    
	hold on;
    
    % FIXME: Plot a heat map here of the 2D nonlinearity
    imagesc(xout.mycoefmatrix);
                
    % Scatter plot the data points
    do_plot_scatter(sel, mdl, xins{ii}{end}, ...
                    mdl.input_stim1, mdl.input_stim2);
           
    do_xlabel('Input1 [-]');
    do_ylabel('Input2 [-]');
    legend off;
    axis tight;
    hold off;
end


end