% function r=adapt_post(beta,pred);
function r=adapt_post(beta,pred);

p1=pred;
p2=conv2(pred(:),[0 0 0 0 0 .4 .3 .2 .1]','same');
%p2=conv2(pred(:),[0 0 0 1/2 1/2]','same');

theta=beta(1);
gamma=beta(2);
A=beta(3);

r=A.*(p1-theta).*(p1>theta) .* exp((p1-p2)*gamma);

