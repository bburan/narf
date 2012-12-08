function im2rriki(h,caption,filename);

if ~exist('h','var'),
   h=gcf;
end

dbopen
mysql('use rriki');


imdata=mysql('SELECT max(id) as id FROM images');

maxid=max([imdata.id 0]);
if ~exist('filename','var'),
   filename=sprintf('image%d.jpeg',maxid+1);
end
filename=['/auto/www/rriki/public/images/',basename(filename)];

if ~exist('caption','var'),
   caption='Matlab figure printed via im2rriki';
end

print(['-f',num2str(h)],'-djpeg','-r75',filename);

sqlinsert('images','caption',caption,...
          'file_name',basename(filename),...
          'content_type','image/jpeg');

