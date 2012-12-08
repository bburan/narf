pid=getenv('MYPID')
[s,w] = unix(['ps h -p ',num2str(pid)]) 
