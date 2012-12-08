% function y=expdecay(beta,x)
% 
% beta=[tau,A,offset]
%
function y=expdecay(beta,x)

tau=beta(1);
A=beta(2);
offset=beta(3);

y=A.*exp(-x/tau)+offset;
