function [ A ] = meanTree_sub( A )
%MEANTREE_SUB(A) estimates the Karcher mean using the Tree merging algorithm.
%This algorithm consists in merging successively the matrices, until there
%is only one matrix left. The pair of matrices merged are chosen according
%to an upside down tree pattern.

%A is a 3D array containing the matrices whose mean has to be estimated
%along the third dimension
%A is the estimation of the Karcher mean

[m,~,n] = size(A);
k = floor(log2(n));
nStar = 2^k;                        %number of nodes to reach after the preprocessing step
count = 0;

%preprocessing step : nodes are merged until the number of remaining nodes
%becomes a power of two
if(nStar ~= n)
    l = n-nStar;
    indx = randperm(n);
    w = zeros(1,nStar);             %weights to give to the nodes, after the preprocessing step
    B = zeros(m,m,nStar);           %matrices remaining after the preprocessing step
    
    %two nodes are chosen and merged
    for i = 1:l                     
        i1 = indx(2*i-1);
        i2 = indx(2*i);
        w(i) = 2;
        B(:,:,i) = MGeom(A(:,:,i1),A(:,:,i2),0.5);
        count = count + 1;
    end
    
    %when enough nodes have been merged, the remainging ones are just
    %copied in the 3D array B.
    indx = indx(2*l+1:end);
    for i = l+1:nStar               
        i1 = indx(i-l);
        B(:,:,i) = A(:,:,i1);
        w(i) = 1;
    end
    
    A = B;
else 
    w = ones(1,n);
end
%end of the preprocessing step


i = k;

while i>0
    
    nCur = 2^(i);                       %number of nodes at the beginning of the current iteration
    nNew = 2^(i-1);                     %number of nodes at the end of the current iteration
    B = zeros(m,m,nNew);                %will contain the matrices remaining at the end of the current iteration
    indx = randperm(nCur);              %order in which the matrices will be merged
    wnew = zeros(1,nNew);               %will contain the weights to give during next iteration to the matrices contained in B 
    
   for j = 1:nNew
       i1 = indx(2*j-1);
       i2 = indx(2*j);
       s = w(i1)+w(i2);
       wcur = w(i2)/s;
       wnew(j) = s;
       B(:,:,j) = MGeom(A(:,:,i1), A(:,:,i2),wcur);
       count = count+1;
   end
   
   A = B;
   w = wnew;
   i = i-1;
   
end

end