% Script to test out the recursomatic thingy
x = [];
x{1}.depth = 1;
x{1}.name = 'A';
x{1}.vals = {1 2 3};
x{2}.depth = 2;
x{2}.name = 'B';
x{2}.vals = {'a' 'b'};
x{3}.depth = 3;
x{3}.name = 'C';
x{3}.vals = {'$' '%'};

fn = @disp;

recursomatic(x, fn, []);

disp('done');