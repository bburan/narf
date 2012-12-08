% SHUFFLE Randomly shuffle the elements in an array.
%    S = SHUFFLE(A) rearranges the elements in A randomly.
%
%    See also RAND, RANDN, SORT.

function shuffled_array = shuffle (input_array)

     [y,ii]=   sort(rand(size(input_array)));

     shuffled_array = input_array(ii);
