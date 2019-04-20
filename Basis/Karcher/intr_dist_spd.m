function y=intr_dist_spd(A,B)
%Authors : B. Jeuris, R. Vandebril and B. Vandereycken
%Implementation coming from : http://people.cs.kuleuven.be/~raf.vandebril/

% intr_dist_spd: Intrinsic distance on the manifold of positive definite
% matrices
%   intr_dist_spd(A,B) returns the intrinsic distance between PD matrices A 
%   and B when the manifold is endowed with inner product inpro_spd.

y=norm(logm(sqrtm(A)\B/sqrtm(A)),'fro');