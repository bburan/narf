function y=expreal(beta,x);

expon=real(beta(1));
offset=beta(2);
gain=real(beta(3));

offset=0;

y=(x-offset).^expon .* gain;

