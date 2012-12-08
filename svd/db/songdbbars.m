

barblack=zeros(10,10,3);

barred=barblack;
barred(:,:,1)=255;

bargreen=barblack;
bargreen(:,:,2)=255;

barblue=barblack;
barblue(:,:,3)=255;

barorange=barblack;
barorange(:,:,1)=255;
barorange(:,:,2)=163;

barpurple=barblack;
barpurple(:,:,1)=196;
barpurple(:,:,3)=196;

barlblue=barblack;
barlblue(:,:,2)=164;
barlblue(:,:,3)=255;

cd /auto/k1/david/html/celldb
imwrite(barblack,'black.jpg','jpeg');
imwrite(barred,'red.jpg','jpeg');
imwrite(bargreen,'green.jpg','jpeg');
imwrite(barblue,'blue.jpg','jpeg');
imwrite(barorange,'orange.jpg','jpeg');
imwrite(barpurple,'purple.jpg','jpeg');
imwrite(barlblue,'lblue.jpg','jpeg');
