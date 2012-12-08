%function isrnewcells
baphy_set_path;


mainbatch=179;
%condbatch=[133 135 136 137 138 156 157 120 154 155];
condbatch=[131:138 150 151 154:196];

fprintf('processing batch %d...\n',mainbatch);
r=dbbatchfill(mainbatch,[],1);
israddruns(r);

for bb=condbatch,
   fprintf('processing batch %d...\n',bb);
   r=dbbatchfill(bb,mainbatch,1);
   israddruns(r);
end


if 0,
   r=dbbatchfill(130,[],1);
   % don't need to measure strfs for 130 (AFM)
   r=dbbatchfill(131,[],1);
   israddruns(r);
   r=dbbatchfill(129,131,1);
   israddruns(r);
   
   r=dbbatchfill(113,[],1);
   israddruns(r);
   % speech includes buggy filtered speech (snr>100)
   r=dbbatchfill(115,[],1);
   israddruns(r);
   r=dbbatchfill(116,[],1);
   israddruns(r);
   r=dbbatchfill(118,[],1);
   israddruns(r);
   r=dbbatchfill(132,[],1);  % speech really clean (snr>101)
   israddruns(r);
   
   r=dbbatchfill(125,[],1);
   israddruns(r);
   r=dbbatchfill(126,[],1);
   israddruns(r);
   
   %r=dbbatchfill(117,115,1);
   r=dbbatchfill(117,126,1);
   israddruns(r);
   r=dbbatchfill(122,113,1);
   israddruns(r);
   r=dbbatchfill(124,126,1);
   israddruns(r);
end

