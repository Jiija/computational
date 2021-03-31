function [W]=ChooseDrawingBox(V,cfinite)
% 
% The function ChooseDrawingBox chooses the cut-drawing box W 
% so that W contains all points - rows of the matrix V.
% The input argument  cfinite  may take values 
% 0 (the visualized set is unbounded) and 1 (the visualized set is bounded).

%% In all the matrices, the first column is abscissa, 
%% the second column is ordinate and the third column is applicate.
%% In the matrices ihV and W, the first row is the lower endpoint and 
%% the second one is the upper endpoint of corresponding boxes.

   % computing interval hull of the points from V
   ihV=[min(V(:,1))  min(V(:,2)) min(V(:,3)); max(V(:,1)) max(V(:,2)) max(V(:,3))];
   center=(ihV(1,:)+ihV(2,:))/2;
   radius=(ihV(2,:)-ihV(1,:))/2;
   rm=max(radius);

   % a slight widening of the cut-drawing box as compared with the interval hull
   if isequal(ihV(1,:),ihV(2,:)) % all components of the hull are degenerate
      W=[center-1; center+1];
   else       % the interval hull has at least one nondegenerate component
      W=[center-1.2*rm; center+1.2*rm];
   end

   % if the visualized set is unbounded, then we enlarge the cut-drawing box W
   if ~cfinite
      radius=(W(2,:)-W(1,:))/2;
      W=W+2*[-radius; radius];
   end

end