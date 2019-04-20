function [ M, tTot, varargout ] = CrudeIterative( A )
%CRUDEITERATIVE(A) : estimation of the Karcher mean by using a combination
%of the Crude and Iterative means. The cycle chosen at each iteration is 
%chosen randomly.
%A is a 3D array containing the matrices whose mean has to be estimated
%along the third dimension
%M is the estimation of the Karcher mean and tTot is the CPU time used by
%the function

%References : 
%Ben Jeuris and Raf Vandebril. Geometric mean algorithms based on harmonic and
%arithmetic iterations. In Geometric Science of Information, pages 785–793. Springer,
%2013.
% Miklós Pálfia. A multivariable extension of two-variable matrix means. SIAM Journal
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
tCumulStop = cputime - tCumulStart;

%the algorithm works with two sequences of nodes : the first one correponds
%to the computation of arithmetic means, at a given iteration, while the
%second one correponds to the computation of harmonic means
B = A;
C = A;

while( max(max(D))>=1e-8 && iter<nIterMax)
    
    iter = iter+1;
%     t = cputime;
    G = randperm(n);
    B_prev = B;
    C_prev = C;
    
    for i = 1:n-1
        IB = inv(B_prev(:,:,G(i)));
        IC = inv(C_prev(:,:,G(i+1)));
        B(:,:,i) = inv((IB+IC)./2);
        C(:,:,i) = (B_prev(:,:,G(i))+C_prev(:,:,G(i+1)))./2;   
    end
    
    IB = inv(B_prev(:,:,G(n)));
    IC = inv(C_prev(:,:,G(1)));
    B(:,:,n) = inv((IB+IC)./2);
    C(:,:,n) = (B_prev(:,:,G(n))+C_prev(:,:,G(1)))./2; 
    
    count = count + n;
    
    tCumulStart = cputime;
    D = zeros(2*n,2*n);
    A = cat(3,B,C);
    for i = 1:2*n-1
        for j = i+1:2*n
            D(i,j) = norm(A(:,:,i)-A(:,:,j),'fro');
        end
    end    
    tCumulStop = tCumulStop + cputime - tCumulStart;

end

M = A(:,:,1);
varargout{1} = iter;
% tRecord = tRecord(1:iter);
if(iter == nIterMax),  fprintf('Maximum number of iterations reached \n'); end

tTot = cputime-tStart - tCumulStop ;
end


