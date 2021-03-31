function [V] = OrientationPoints(S)
%
% The function OrientationPoints constructs the matrix  V  
% of orientation points from the matrix  S  of boundary intervals.

   V=S(:,1:2); 

   xbf=isfinite(V(:,1));
   ybf=isfinite(V(:,2));
   bf=~min(xbf,ybf); 

   V(bf,:)=[]; 
   V=NonRepeatRows(V);

end