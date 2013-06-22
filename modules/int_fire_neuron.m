function m = int_fire_neuron(args)
% simulates a very simple integrate-and-fire neuron with single
% excitatory and inhibitory inputs

% Module fields that must ALWAYS be defined
m = [];
m.mdl = @int_fire_neuron;
m.name = 'int_fire_neuron';
m.fn = @do_int_fire_neuron;
m.pretty_name = 'Integrate-and-fire neuron';
m.editable_fields = {'input_stim1', 'input_stim2', 'input_resp', 'time',...
                     'output', 'Vrest', 'V0', 'gL','rectify_inputs','subsample'};
m.isready_pred = @isready_always;

% Module fields that are specific to THIS MODULE
m.input_stim1 = 'stim1';
m.input_stim2 = 'stim2';
m.input_resp = 'respavg';
m.time = 'stim_time';
m.output = 'stim';
m.Vrest=0;
m.V0=[40 -20];
m.gL=25;
m.rectify_inputs=0;
m.subsample=10;

% Optional fields
m.plot_fns = {};
m.auto_plot = @do_plot_irn_output;
m.plot_fns{1}.fn = @do_plot_irn_output; 
m.plot_fns{1}.pretty_name = 'IRN Input/Output';
m.plot_fns{2}.fn = @do_plot_irn_input; 
m.plot_fns{2}.pretty_name = 'IRN Inputs';

% Overwrite the default module fields with arguments 
if nargin > 0
    m = merge_structs(m, args);
end

% ------------------------------------------------------------------------
% Methods


function x = do_int_fire_neuron(mdl, x, stack, xxx)    
    
    fns = fieldnames(x.dat);
    
    Vrest=mdl.Vrest;
    V0=mdl.V0;
    gL=mdl.gL;
    fs=stack{1}{1}.raw_stim_fs;
    
    for ii = 1:length(fns)
         sf=fns{ii};
         % Compute the size of the INPUT
         gE=x.dat.(sf).(mdl.input_stim1);
         gI=x.dat.(sf).(mdl.input_stim2);
         if mdl.rectify_inputs,
             gE(gE<0)=0;
             gI(gI<0)=0;
         end
         
         % Compute the size of the INPUT
         [T, S] = size(gE);
         
         V = zeros(T, S);        
         
         % Filter!
         V(1,:)=Vrest;
         for tt = 2:T
             tV=V(tt-1,:);
             for ii=1:mdl.subsample,
                 
                 dV= (-gE(tt,:).*(tV-V0(1)) - gI(tt,:).*(tV-V0(2)) ...
                      - gL*(tV-Vrest))./(fs./2)./mdl.subsample;
                 tV=tV+dV;
             end
             V(tt,:)=tV;
             %V(tt,V(tt,:)>max(V0))=max(V0);
             %V(tt,V(tt,:)<min(V0))=min(V0);
             
             %if (V(tt,1)>10000),
             %if tt>50,
             %    sfigure(2)
             %plot([V(1:tt,1) gE(1:tt,1) gI(1:tt,1)]);
             %fprintf('%10.3f %10.3f %10.3f %10.3f \n',[dV(1) V(tt,1) gE(tt,1) gI(tt,1)]);
             %keyboard
             %end
         end
         
         % The output is the sum of the filtered channels
         x.dat.(sf).(mdl.output) = V; 
    end
    %keyboard
    
end

function do_plot_irn_output(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    do_plot(xouts, mdls{1}.time, ...
            {mdls{1}.input_stim1 mdls{1}.input_stim2 mdls{1}.output}, ...
            sel, 'Time [s]', 'IRN Output [-]');
end

function do_plot_irn_input(sel, stack, xxx)
    [mdls, xins, xouts] = calc_paramsets(stack, xxx(1:end-1)); 
    sel.chan_idx = []; % when chan_idx is empty, do_plot plots all channels
    
    do_plot(xouts, mdls{1}.time, {mdls{1}.input_stim1 mdls{1}.input_stim2} , ...
            sel, 'Time [s]', 'IRN Inputs [-]');
end

end