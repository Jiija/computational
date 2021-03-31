function [V,P1,P2,P3,P4]=Abs22D(A,B,c,d)
%
% The function  Abs22D  visualizes the solution set 
% for the system of inequalities  |Ax-c| <= B|x|+d, where
%   A,B  are real matrices with  m  rows and  2 columns,
%   c,d  are real vectors of the length  m.
%
% Output arguments:
% V is a matrix of orientation points.
%   In the package, the orientation points are understood as vertices 
%   of intersections of the solution set with orthants.
% Pk, k=1,2,3,4.
%   If intersection of the solution set with the k-th orthant is bounded,
%   Pk is a closed clockwise path around the piece of the solution set 
%   in the k-th orthant. 
%   Otherwise, Pk is the same path but constructed around intersection 
%   of the cut box B with the piece of the solution set in the k-th orthant. 
%   The program chooses the cut box  B  automatically.
%
% Written by Irene A. Sharaya (Institute of Computational Technologies, Novosibirsk)
% Date:  December 12, 2013

    % turn to interval inclusion [uC,oC] x \subseteq [ud,od] 
    % in Kaucher arithmetic
    uC=A+B;
    oC=A-B;
    ud=c-d;
    od=c+d;

    % visualize the solution set 
    % to the inclusion [uC,oC] x \subseteq [ud,od]
    % and, if necessary, display V and Pk (k=1,2,3,4)
    if nargout < 1
        Cxind2D(uC,oC,ud,od);
    end
    if nargout == 1
        [V]=Cxind2D(uC,oC,ud,od);
    end
    if nargout > 1
        [V,P1,P2,P3,P4]=Cxind2D(uC,oC,ud,od);
    end
end
