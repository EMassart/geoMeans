function [ M, tTot, varargout] = meanStar( A )
%MEANSTAR(A) estimates the Karcher mean of a set of positive definite matrices
%with a Star algorithm. The Star algorithm is an iterative algorithm in which
%the matrices are represented as the nodes of a graph. At each iteration, a
%Star graph passing through these nodes is built, and the two nodes incident to each
%edge are replaced by their mean. The replacement of the nodes is immediately 
%performed, meaning that the replacement of the nodes at the beginning of
%the iteration will have an impact on the remaining nodes during the same
%iteration.
%A is a 3D array containing the matrices whose mean has to be estimated
%along the third dimension
%M is the estimation of the Karcher mean and tTot is the CPU time used by the function


tStart = cputime;
[~,~,n] = size(A);
D = zeros(n,n);
nIterMax = 1000;
iter = 0;
count = 0;
ok = 0;                                         %these two variables will check that each node has been at least once the center of the star
                                                %before stopping the algorithm
isTaken = zeros(1,n);
%tRecord = zeros(1,nIterMax);

%the computation of the distance separating the matrices is not formally
%required by the algorithm, therefore the CPU time required will be removed
%from the total 

tCumulStart = cputime;
for i = 1:n-1
    for j = i+1:n
        D(i,j) = norm(A(:,:,i)-A(:,:,j),'fro');
    end
end
tCumulStop = cputime - tCumulStart;


while ((max(max(D))>= 1e-8 || ~ok )&& iter<nIterMax)
    iter = iter+1;
%     t = cputime;
    c = randi(n);                                %center of the star
    isTaken(c) = isTaken(c)+1;
    extNodes = [1:c-1 c+1:n];
    for i=extNodes
        A(:,:,c) = MGeom(A(:,:,c),A(:,:,i),0.5);
        A(:,:,i) = A(:,:,c);
    end
    
    count = count + n-1;
    if(isempty(find(isTaken==0,1))), ok = 1; end
%     tRecord(iter) = (cputime-t)/(n-1);
    
    tCumulStart = cputime;
    D = zeros(n,n);
    for i = 1:n-1
        for j = i+1:n
            D(i,j) = norm(A(:,:,i)-A(:,:,j),'fro');
        end
    end
    tCumulStop = tCumulStop + cputime - tCumulStart;
end
    
M = A(:,:,1);
varargout{1} = count;
% tRecord = tRecord(1:iter);
% varargout{2} = tRecord;

if(iter == nIterMax),  fprintf('Maximum number of iterations reached \n'); end

tTot = cputime - tStart - tCumulStop;

end