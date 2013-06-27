function m = nonparm_nonlinearity_2d(args)
% Applies a 2D nonparametric static nonlinear function to the inputs

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @nonparm_nonlinearity_2d;
m.name = 'nonparm_nonlinearity_2d';
m.fn = @do_np_nonlinearity_2d;
m.pretty_name = 'Nonparametric Nonlinearity 2D';
m.editable_fields = {'input_stim1', 'input_stim2', 'input_resp', 'time',...
                     'output', 'bincount', 'smoothwindow'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input_stim1 = 'stim1';
m.input_stim2 = 'stim2';
m.input_resp = 'respavg';
m.time = 'stim_time';
m.output = 'stim';
m.bincount = 20;
m.smoothwindow=1;

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
    pred1=[]; pred2=[]; resp=[];
    for sf = x.training_set,
        sf=sf{1};
        % It's an error if there is more than one channel:
        if size(x.dat.(sf).(mdl.input_stim1), 3) > 1 || ...
           size(x.dat.(sf).(mdl.input_stim2), 3) > 1 || ...
           size(x.dat.(sf).(mdl.input_resp), 3) > 1,
            error('NPNL is a 1D method!');
        end
        pred1=cat(1,pred1,x.dat.(sf).(mdl.input_stim1)(:));
        pred2=cat(1,pred2,x.dat.(sf).(mdl.input_stim2)(:));
        resp=cat(1,resp,x.dat.(sf).(mdl.input_resp)(:));
    end
    
    keepidx=find(~isnan(resp));
    pred1=pred1(keepidx);
    pred2=pred2(keepidx);
    resp=resp(keepidx);
    
    bincount=mdl.bincount;
    smoothwindow=mdl.smoothwindow;
    pp=zeros(bincount,2);
    rr=zeros(bincount,bincount);
    rre=zeros(bincount,bincount);
    
    if nansum(resp)>0,
        [ss1,si1]=sort(pred1);
        tb=bincount;
        b1=[];
        
        while length(b1)<bincount && tb<length(si1),
            tb=tb+1;
            edges1=round(linspace(1,length(si1)+1,tb));
            if std(ss1(edges1(1:(end-1))))>0
                [b1,ui]=unique(ss1(edges1(1:(end-1)))');
            else
                b1=ss1(edges1(1:(end-1)))';
                ui=1:length(b1);
            end
               
        end
        edgeend=edges1(end);
        edges1=edges1(ui);
        
        if length(b1)<bincount,
            b1=[b1 repmat(b1(end),[1 bincount-length(b1)])];
            edges1=[edges1 repmat(edges1(end),[1 bincount-length(edges1)])];
        end
        edges1=[edges1 edgeend];
        b1(bincount+1)=ss1(end)+eps;
        
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
        b2(bincount+1)=ss2(end)+eps;
        
        for bb=1:bincount,
            pp(bb,1)=mean(ss1(edges1(bb):(edges1(bb+1)-1)));
            pp(bb,2)=mean(ss2(edges2(bb):(edges2(bb+1)-1)));
        end
        
        for bb1=1:bincount,
            for bb2=1:bincount,
                ff=find(pred1>=b1(bb1) & pred1<b1(bb1+1) &...
                        pred2>=b2(bb2) & pred2<b2(bb2+1));
                rr(bb1,bb2)=mean(resp(ff));
                nn=length(ff);
                if nn,
                    %rre(bb1,bb2)=std(resp(ff))./sqrt((nn+(nn==0)));
                    rre(bb1,bb2)=nn;
                end
            end
        end
        pp(isnan(pp))=0;
        rr(isnan(rr))=0;
        
        rr=gsmooth(rr,[smoothwindow smoothwindow]);
        
    end
    
    phi={pp,rr};
    outbinserr=rre;
    %keyboard
 end

function x = do_np_nonlinearity_2d(mdl, x, stack, xxx)    
    
    [phi2d, ~] = init_nonparm_nonlinearity_2d(mdl, x, stack, xxx);
    
    fns = fieldnames(x.dat);
    for ii = 1:length(fns)
        sf = fns{ii};
        [T, S, C] = size(x.dat.(sf).(mdl.input_stim1));
        y = raw_nl(phi2d, [x.dat.(sf).(mdl.input_stim1)(:) ...
                     x.dat.(sf).(mdl.input_stim2)(:)]);                      
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
    imagesc(xout.mycoefmatrix{1}(:,1),...
            xout.mycoefmatrix{1}(:,2),...
            xout.mycoefmatrix{2});
    hold on;
    % Scatter plot the data points
    %do_plot_scatter(sel, xins, ...
    %                 mdl.input_stim1, mdl.input_stim2);
   
    colorbar
    
    do_xlabel('Input1 [-]');
    do_ylabel('Input2 [-]');
    legend off;
    axis tight;
    hold off;
end


end