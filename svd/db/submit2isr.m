% function submit2isr(queueidx)

function submit2isr(queueidx)

isrhome='/homes/svd/';

queuefile=sprintf('%squeue/mat.%d.sh',isrhome,queueidx);

cmd=['ssh svd@seil.umd.edu qsub ',queuefile];
unix(cmd);
