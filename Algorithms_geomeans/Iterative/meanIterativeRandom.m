function [ M, tTot, varargout ] = meanIterativeRandom( A )
%MEANITERATIVERANDOM(A) estimates the Karcher mean using the Iterative algorithm proposed in [1].
%The cycle is chosen randomly at each iteration.
%A is a 3D array containing the matrices whose mean has to be estimated
%along the third dimension
%M is the estimation of the Karcher mean and tTot is the CPU time used by
%the function

%References : 
%[1] Miklós Pálfia. A multivariable extension of two-variable matrix means. SIAM Journal
%on Matrix Analysis and Applications, 32(2):385–393, 2011.

tStart = cputime;
[~,~,n] = size(A);
D = zeros(n,n);
nIterMax = 1000;
iter = 0;
count = 0;
% tRecord = zeros(1,nIterMax);

%the computation of the distance separating the matrices is not formally
%required by the algorithm, therefore the CPU time required will be removed
%from the total 
tCumulStart = cputime; 
for i = 1:n-1
    for j = i+1:n
        D(i,j) = norm(A(:,:,i)-A(:,:,j),'fro');
    end
end
D = D+D';
tCumulStop = cputime - tCumulStart;


while( max(max(D))>=1e-8 && iter<nIterMax)
    
    iter = iter+1;
%     t = cputime;
    G = randperm(n);
    B = A;
    
    for i = 1:n-1
        A(:,:,i) = MGeom(B(:,:,G(i)), B(:,:,G(i+1)), 0.5);
    end
    
    A(:,:,n) = MGeom(B(:,:,G(n)), B(:,:,G(1)), 0.5);
    count = count + n;
    
%     tRecord(iter) = (cputime-t)/n;
    
    tCumulStart = cputime; 
    D = zeros(n,n);
    for i = 1:n-1
        for j = i+1:n
            D(i,j) = norm(A(:,:,i)-A(:,:,j),'fro');
        end
    end   
    D = D+D';
    tCumulStop = tCumulStop + (cputime - tCumulStart);
end

M = A(:,:,1);
varargout{1} = iter;
% tRecord = tRecord(1:iter);
if(iter == nIterMax),  fprintf('Maximum number of iterations reached \n'); end

tTot = cputime-tStart-tCumulStop;
end





