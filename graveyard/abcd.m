% Test of an ABCD System

damping = 1;  % Must be greater than 1
osc_frq = 5;    % In Hz
input_frq = 5;  % in Hz
input_gain = 1;
output_gain = 1;

x_0 = [0 1]';

A_inertia = zeros(2,2);
A_damping = [0 0; 0 1]*(-damping);
A_spin    = skewdec(2,0) * osc_frq*2*pi;
A = A_inertia + A_damping + A_spin;
B = [0 0]';
C = [1 0; 0 1];
D = [0 0]';

%sys = ss(A,B,C,D);
%DelayT(1) = struct('delay',0.5,'a',0,'b',2,'c',1,'d',0);
%DelayT(2) = struct('delay',1.2,'a',-1,'b',0,'c',0,'d',0);

% Basically, specify a new ABCD at a particular delay point. 
delayterms = struct('delay',0.5,'a',[0 0; 0 0],'b',[0 1]','c', [0 0; 0 0],'d', [0 0]'); 
delayterms(2) = struct('delay', 3,'a',[0 0; 0 0],'b',[0 1]','c', [0 0; 0 0],'d', [0 0]'); 

sys=delayss(A,B,C,D,delayterms); % With delays

t = 0:0.001:10.0;

u = cos(input_frq*2*pi*t);

figure; impulse(sys, t);

%figure; step(sys, t);

%figure; lsim(sys, u, t, x_0);

%figure; bode(sys);     
