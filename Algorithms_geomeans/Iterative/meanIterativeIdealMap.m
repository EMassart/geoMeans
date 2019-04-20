function [ M, tTot, varargout] = meanIterativeIdealMap( A )
%MEANITERATIVEIDEALMAP(A) estimates the Karcher mean using the Iterative algorithm proposed in [1].
%The cycle is chosen at each iteration according to the IdealMap procedure.
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
count = 0;                                      %records the number of means of two matrices evaluated
% tRecord = zeros(1,nIterMax);


for i = 1:n-1
    for j = i+1:n
        D(i,j) = norm(A(:,:,i)-A(:,:,j),'fro');
    end
end
D = D+D';


while( max(max(D))>=1e-8 && iter<nIterMax)
    
    iter = iter+1;
%     t = cputime;
    G = idealMap(D);                            %computes the new cycle from the distance matrix D
    B = A;
    
    for i = 1:n-1
        A(:,:,i) = MGeom(B(:,:,G(i)), B(:,:,G(i+1)), 0.5);
    end
    
    A(:,:,n) = MGeom(B(:,:,G(n)), B(:,:,G(1)), 0.5);
    count = count + n;
    
%     tRecord(iter) = (cputime-t)/n;
    
    D = zeros(n,n);
    for i = 1:n-1
        for j = i+1:n
            D(i,j) = norm(A(:,:,i)-A(:,:,j),'fro');
        end
    end    
    D = D+D';

end

M = A(:,:,1);
varargout{1} = iter;
% tRecord = tRecord(1:iter);
% varargout{2} = mean(tRecord(1:iter));

if(iter == nIterMax),  fprintf('Maximum number of iterations reached \n'); end

tTot = cputime-tStart;

end





