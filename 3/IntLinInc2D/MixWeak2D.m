function [V,P1,P2,P3,P4]=MixWeak2D(infA,supA,infb,supb,relations)
%
% The function  MixWeak2D  visualizes  
% the united solution set  (= weak solution set)
% for the system of interval linear relations  
%                A x relations b,
% where 
%   A  is an interval matrix with  m  rows and  2 columns;
%   b  is an interval vector of the length  m;
%   rows with relation 'equal' have ÀÅ-order of quantifiers.
%
% Input arguments:
% infA, supA  are  left and right endpoints of the interval matrix  A  
%    respectively;
% Aq is a quantifier matrix with  m  rows and  2 columns, 
%    its elements are letters 'A' and 'E', which correspond to universality 
%    and existence quantifiers in the description of the quantifier solution.
% infb, supb  are  left and right endpoints of the interval vector  b  
%    respectively;
% bq is a quantifier vector of the length  m, whose elements are letters 
%    'A' and 'E' (similar to  Aq);
% relations  is vector of the length m, whose elements are signs '=', '>', '<',
%    that designate relations "equal", "greater than or equal" and 
%    "less than or equal" respectively.
%
% Output arguments:
% V is a matrix of orientation points.
%    In the package, the orientation points are understood as vertices 
%    of intersections of the solution set with orthants.
% Pk, k=1,2,3,4.
%    If intersection of the solution set with the k-th orthant is bounded,
%    Pk is a closed clockwise path around the piece of the solution set 
%    in the k-th orthant. 
%    Otherwise, Pk is the same path but constructed around intersection 
%    of the cut box B with the piece of the solution set in the k-th orthant. 
%    The program chooses the cut box  B  automatically.

    % turn to interval inclusion [uC,oC] x \subseteq [ud,od] 
    % in Kaucher arithmetic
    uC=supA;
    oC=infA;
    ud=infb;
    od=supb;
    m=size(infb,1);
    relLE=[relations==ones(m,1)*'<'];
    ud(relLE)=-Inf;
    relGE=[relations==ones(m,1)*'>'];
    od(relGE)=Inf;

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
