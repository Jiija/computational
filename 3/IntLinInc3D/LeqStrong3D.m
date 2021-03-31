function [V]=LeqStrong3D(infA,supA,infb,supb,OrientPoints,transparency,varargin)
%
% The function  LeqStrong3D  visualizes the strong solution set
% for the system of interval linear inequalities  A x <= b, where
%   A  is an interval matrix with  m  rows and  3 columns,
%   b  is an interval vector of the length  m.
%
% Input arguments:
%   infA, supA  are  left and right endpoints of the interval matrix  A  
%     respectively;
%   infb, supb  are  left and right endpoints of the interval vector  b  
%     respectively;
%   OrientPoints is a parameter that controls drawing of the orientation points; 
%     if it has the value 0, the program does not plot these points; 
%     if it has the value 1, the orientation points are plotted; 
%   transparency is a parameter responsible for the transparency 
%     of the real faces, its values may be 0 or 1;   
%     0 means the absence of the transparency, 
%     1 means that the real faces are transparent.
%
% Optional input arguments:
%   varargin is a prescribed cut box, 
%     it is inputted as 6 numbers, xb, xe, yb, ye, zb, ze, where 
%     xb,yb,zb are coordinates of the lower endpoint of the box, 
%     xe,ye,ze are coordinates of the upper endpoint of the box.
%
% Output arguments:
%   V is a matrix of the orientation points.
%     In the package, the orientation points are understood as vertices 
%     of intersections of the solution set with orthants.

   % turn to interval inclusion [uC,oC] x \subseteq [ud,od] 
   % in Kaucher arithmetic
   uC=infA;
   oC=supA;
   m=size(infb,1);
   ud=-ones(m,1)*Inf;
   od=infb;

   % visualize the set of formal solutions for
   % the interval inclusion [uC,oC] x \subseteq [ud,od] 
   var=(length(varargin)==6);
   if var 
      [V]=Cxind3D(uC,oC,ud,od,OrientPoints,transparency,varargin{1},varargin{2},varargin{3},varargin{4},varargin{5},varargin{6});
   else
      [V]=Cxind3D(uC,oC,ud,od,OrientPoints,transparency);
   end

end
