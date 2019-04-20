function [X,tTot,varargout]=cheap(A_input)

% [X,iter]=CHEAP(A1,A2,...,AP) computes the mean defined by Bini and Iannazzo
%  in [1]
% [X,iter]=CHEAP(A{1:p}) for a cell-array input
%
% varargin: positive definite matrix arguments A1,...,AP
% X: the Cheap mean of A1,...,AP
% iter: the number of iterations needed by the outer iteration
% 
% References
% [1] D.A. Bini and B. Iannazzo, "A note on computing matrix geometric 
% means", Adv. Comput. Math., 35-2/4 (2011), pp. 175-192.

tStart = cputime;
[n,~,p]=size(A_input);            %size and number of input matrices
tol=1e-9; maxiter=1000;
d = zeros(p,p);

for h=1:p
  A{h}=A_input(:,:,h);
end

for k=1:maxiter
  for h=1:p
    R{h}=chol(A{h});
    RI{h}=inv(R{h});
  end
  for h=1:p
    S=zeros(n);
    for ell=1:p
      if (ell~=h)
        % computes S=S+logm(RI{h}'*A{ell}*RI{h})
        Z=R{ell}*RI{h};
        [U V]=schur(Z'*Z);
        T=U*diag(log(diag(V)))*U';
        S=S+(T+T')/2;
      end
    end
    [U V]=schur(1/p*S);
    T=diag(exp(diag(V))).^(1/2)*U'*R{h};
    A1{h}=T'*T;
  end
  
  for h = 1:p-1
      for h_loc = h+1:p
          d(h,h_loc) = dist(A1{h},A1{h_loc});
      end
  end
  
  ni=norm(A1{1}-A{1})/norm(A{1});
  if (max(max(d))<1e-8)
    iter=k; break;
  end
  if (k==maxiter)
    disp('Max number of iterations reached');
    iter=k; break;
  end
  for h=1:p
    A{h}=A1{h};
  end
    
end

X = A1{1};
varargout{1} = iter;
tTot = cputime - tStart;