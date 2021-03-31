function [V,P1,P2,P3,P4]=Cxind2D(uC,oC,ud,od)
%
% The function Cxind2D visualizes the set of formal solutions 
% for the interval inclusion [uC,oC] x \subseteq [ud,od] 
% in Kaucher arithmetic, 
% where:
%   uC and oC  are real matrices with  m  rows and  2 columns
%     (left and right endpoints of the interval matrix [uC,oC] respectively),
%   ud and od  are vectors of the length m, 
%     whose components may be real numbers or -Inf and +Inf 
%     (left and right endpoints of the interval vector [ud,od] respectively),
%   x is the vector of unknowns, having the length 2.
%
% Output arguments:
% If intersection of the solution set with the k-th orthant is bounded,
%   Pk is a closed clockwise path around the piece of the solution set 
%   within the k-th orthant. 
%   Otherwise, Pk is the same path but around intersection of the cut box B 
%   with the piece of the solution set within the k-th orthant 
%   (the program chooses the cut box  B  automatically).
% V is the matrix of orientation points 
%   (Throughout the package, orientation points are vertices 
%    of intersections of the solution set with orthants).
%
% Written by Irene A. Sharaya (Institute of Computational Technologies, Novosibirsk)
% Date:  March 7, 2012

   P1=zeros(0,2);
   P2=zeros(0,2);
   P3=zeros(0,2);
   P4=zeros(0,2);

   % The solution set in the orthant  k  is described by a system of
   % inequalities  Ak*x >= bk  with the following matrix  Ak  and vector  bk 
   A1(:,:)=[uC; -oC; 1 0; 0 1]; 
   A2(:,:)=[oC(:,1) uC(:,2); -uC(:,1) -oC(:,2); -1 0; 0 1];
   A3(:,:)=[oC; -uC; -1 0; 0 -1];
   A4(:,:)=[uC(:,1) oC(:,2); -oC(:,1) -uC(:,2); 1 0; 0 -1];
   b1=[ud; -od; 0; 0];
   b2=b1;
   b3=b1;
   b4=b1;


   % For each orthant  k:
   % 1) remove, from the above systems, rows for which the solution 
   %    set is obviously empty, and rows for which the solution set 
   %    is the whole of the plane R^2;
   % 2) then compute the matrix Sk of the boundary intervals for 
   %    the solution set in the orthant k.

   [A1,b1,cnmty]=ClearRows(A1,b1);
   if cnmty
      [S1]=BoundaryIntervals(A1,b1);
   else
      [S1]=zeros(0,5);
   end

   [A2,b2,cnmty]=ClearRows(A2,b2);
   if cnmty
      [S2]=BoundaryIntervals(A2,b2);
   else
      [S2]=zeros(0,5);
   end

   [A3,b3,cnmty]=ClearRows(A3,b3);
   if cnmty
      [S3]=BoundaryIntervals(A3,b3);
   else
      [S3]=zeros(0,5);
   end

   [A4,b4,cnmty]=ClearRows(A4,b4);
   if cnmty
      [S4]=BoundaryIntervals(A4,b4);
   else
      [S4]=zeros(0,5);
   end

   % put together all the matrices Sk into the matrix S 
   S=[S1; S2; S3; S4];

   if size(S,1)==0
      disp('Solution set is empty (it does not have boundary intervals)') 
      V=[];
      return;
   end

   % compile a matrix  V  of orientation points 
   V=OrientationPoints(S);

   % choose a drawing box W   
   [W]=DrawingBox(V);

   % If the solution set is unbounded
   % (= the matrix S has infinite elements),
   % it is necessary 
   % 1) to extend the drawing box W (so that the infinity is perceivable),
   % 2) to choose the cut box B and 
   % 3) to recalculate the matrices Sk having infinity elements.

   bounded=1;
   if ~all(all(isfinite(S)))
   bounded=0;

      % to extend the drawing box W and to choose the cut box B
      [W,B]=CutBox(W);

      % add, to the system  Ak*x >= b, two inequalities from the cut box  B
      % and recalculate the matrices  Sk  having infinity elements

      if ~all(all(isfinite(S1)))
          A1=[A1;    -1 0;    0 -1];
          b1= [b1;  -B(2,1); -B(2,2)];
          [S1]=BoundaryIntervals(A1,b1);
      end

      if ~all(all(isfinite(S2)))
          A2=[A2;    1 0;    0 -1];
          b2= [b2;  B(1,1); -B(2,2)];
          [S2]=BoundaryIntervals(A2,b2);
      end

      if ~all(all(isfinite(S3)))
          A3=[A3;    1 0;    0 1];
          b3= [b3;  B(1,1); B(1,2)];
          [S3]=BoundaryIntervals(A3,b3);
      end

      if ~all(all(isfinite(S4)))
          A4=[A4;   -1 0;    0 1];
          b4= [b4; -B(2,1); B(1,2)];
          [S4]=BoundaryIntervals(A4,b4);
      end
      
    end
      
    % construct a closed path Pk around the piece of the solution set 
    % in the k-th orthant from the matrices Sk (k=1,2,3,4).
    % The rows of Pk are vertices of polytope according to clockwise movement.
    [P1] = Intervals2Path(S1);
    [P2] = Intervals2Path(S2);
    [P3] = Intervals2Path(S3);
    [P4] = Intervals2Path(S4);

    % draw the solution set
    SSinW

    % prepare the matrix V for displaying and 
    % display the number of orientation points
    V=V';        
    fprintf('Number of orientation points = %d\n\n',size(V,2)) 

end                
