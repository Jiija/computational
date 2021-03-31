function [W]=DrawingBox(V)
%
% The function  DrawingBox  chooses the drawing box  W, so that 
% it contains all the points whose coordinates are rows of the matrix  V.

   % Interval hull of the points from V.

   ihV=[min(V(:,1)) min(V(:,2)); max(V(:,1)) max(V(:,2))];
   center=[(ihV(1,1)+ihV(2,1))/2 (ihV(1,2)+ihV(2,2))/2];
   radius=[(ihV(2,1)-ihV(1,1))/2 (ihV(2,2)-ihV(1,2))/2];
   rm=max(radius);

   % The drawing box  W  must be wider then the interval hull
   % and depends on its degenerate components.

   if isequal(ihV(1,:),ihV(2,:))   % all components are degenerate
      W=[center(1)-1 center(2)-1; center(1)+1 center(2)+1];
   elseif all(ihV(1,:)<ihV(2,:))   % all components are nondegenerate
      W=ihV+.2*[-radius; radius];
   else                            % one component is degenerate
      W=[center(1)-1.2*rm  center(2)-1.2*rm ; center(1)+1.2*rm  center(2)+1.2*rm];
   end

end