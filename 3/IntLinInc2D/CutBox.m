function [W,B]=CutBox(W)
%
% The function CutBox extends the drawing box W 
% (so that the infinity is perceivable)
% and chooses the cut box B.
 
   radius=[(W(2,1)-W(1,1))/2 (W(2,2)-W(1,2))/2];
   W=W+.4*[-radius; radius];

   radius=[(W(2,1)-W(1,1))/2 (W(2,2)-W(1,2))/2];
   B=W+2*[-radius; radius];

end