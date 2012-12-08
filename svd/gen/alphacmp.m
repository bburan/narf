% function r=alphacmp(s1,s2);
%
% r=-1  if s1<s2
%  = 0  if s1==s2
%  = 1  if s1>s2
%
function r=alphacmp(s1,s2);

l1=length(s1);
l2=length(s2);

if l1>l2,
   ts1=s1(1:l2);
   ts2=s2;
else
   ts1=s1;
   ts2=s2(1:l1);
end

d=double(ts1)-double(ts2);

if sum(abs(d))==0,
   if l1>l2,
      r=1;
   elseif l1<l2,
      r=-1;
   else
      r=0;
   end
elseif d(min(find(d~=0)))<0,
   r=-1;
else
   r=1;
end
