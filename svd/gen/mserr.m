% function e=mserr(fitdata,obsdata);
function e=mserr(fitdata,obsdata);

d=sum(sum(sum(abs(fitdata-obsdata).^2)));
s=sum(sum(sum(abs(fitdata).^2)));
e=d/s;
