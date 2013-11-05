function f = fibo(n)
% Recursive fibonacci. Slow. 
%
% >> tic; fibo(25); toc
% Elapsed time is 1.634411 seconds.

if n <= 1  
    f = 1;
else
    f = fibo(n-1) + fibo(n-2);
end

end