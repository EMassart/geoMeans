function [ma] = idealMap(D)
%IDEALMAP(D) returns a cycle used in the Iterative mean algorithm [1].
%The cycle is built with a greedy algorithm is order to try to maximize its length.
%D is a matrix containing the lengths of the edges (which are usually the 
%distances between the matrices corresponding to the nodes).
%ma is a permutation of the nodes corresponding to the cycle returned

%References : 
%[1] Miklós Pálfia. The riemann barycenter computation and means of several matrices.
%Int. J. Comput. Math. Sci, 3(3):128–133, 2009.

n = size(D,1);
d = n*(n-1)/2;
i = 1;
j = 2;
r = zeros(d,3);                     %r will contain for each pair of matrices the distance separating them (first column) 
                                    %as well as the indices of the two matrices (second and third columns)
                                    
ma = zeros(n,1);                    %sequence of nodes to return (cycle obtained)

%build the matrix r starting from D
for k = 1:d
    r(k,1) = D(i,j);
    r(k,2) = i;
    r(k,3) = j;
    if j == n
        j = n-i+2;
        i = 1;
    else
        i = i+1;
        j = j+1;
    end
end

%the nodes are then chosen in order to favour the most distant matrices
%(greedy algorithm)
[~,I] = sort(r(:,1),'descend');
r = r(I,:);

%the nodes of the first pair appearing in r will be the two first nodes
%used to build the cycle
ma(1) = r(1,2);
j = r(1,3);                             %the variable j will contain at each moment the last node inserted in the cycle

%all the pairs in which the first node appears are removed from the
%matrix r in order to avoid prematurly closing the cycle
I2 = (r(:,2)== ma(1));
r(I2,:) = [];
I3 = (r(:,3)== ma(1));
r(I3,:) = [];

%try to find the next node to insert in the cycle : the one which is the
%most distant to the last node entered
for k = 2:n-1
    indx = find((r(:,2)==j | r(:,3)==j),1);             %indices of the pairs in which the last node insterted appears
    
    %among the remaining pairs, the first one is chosen (the pairs are ordered by decreasing distance between the nodes)
    if r(indx,2) == j
        j = r(indx,3);                          
        ma(k) = r(indx,2);                              
    else 
        j = r(indx,2);
        ma(k) = r(indx,3);
    end
    
    %all the pairs in which the previous node appears are removed from the
    %matrix r in order to avoid prematurly closing the cycle
    I2 = (r(:,2)== ma(k));
    r(I2,:) = [];
    I3 = (r(:,3)== ma(k));
    r(I3,:) = [];
end
ma(n) = j;


end