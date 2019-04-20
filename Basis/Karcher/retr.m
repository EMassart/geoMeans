function Y=retr(X,E,retract)
%Authors : B. Jeuris, R. Vandebril and B. Vandereycken
%Implementation coming from : http://people.cs.kuleuven.be/~raf.vandebril/
% retr: retraction
%   retr(X,E,retract) returns the retraction from the point X in the 
%   direction of tangent vector E.
%   The parameter retract can take one of the following options (strings):
%       'symm'    -> first order retraction (natural for inpro_symm)
%       'spd'     -> natural retraction for inpro_spd
%       'approx2' -> second order retraction 
%                    (second order approximation to spd)
%   When the given option for retract is not recognized, 'symm' is taken as
%   default.

switch lower(retract)
    case 'symm'
        Y=X+E;
    case 'spd'
        Y=retr_spd(X,E);
    case 'approx2'
        Y=retr_approx2(X,E);
    otherwise
        display('not a valid retraction: switching to symm')
        Y=X+E;
end

Y=(Y+Y')/2;

end


function Y=retr_spd(X,E)

Y=X*expm(X\E);

end


function Y=retr_approx2(X,E)

Y=X+E+(E/X*E)/2;

end