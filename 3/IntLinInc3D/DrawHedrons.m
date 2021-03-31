function []=DrawHedrons(A,b,m,transparency,var)
%
% The function  DrawHedrons  draws the intersection 
% îf the solution set to the system  A x >= b  with 
% the plane  A(i,:) x = b(i)  for each required i  
% (i.e., it draws required i-suppots).

   mn=size(b,1);

   for i=4:mn  % do not draw i-supports 
               % that correspond to coordinate planes

      % switch from the initial coordinates system (0,x1,x2,x3)
      % to the intrinsic coordinates system (pr,y1,y2) of plane i
      % (more exactly, of the plane A(i,:)x=b(i))
      [pr,At,bt,y1,y2]=ChangeVariable(A,b,i);

      % delete zero rows from  At
      [At,bt,cnmtyt]=ClearZeroRows(At,bt);
    
      if ~cnmtyt
         continue; % The plane i does not intersect the solution set.
      end

      [St]=BoundaryIntervals(At,bt);

      if size(St,1)>0
         [Pt] = Intervals2Path(St);

         % rewrite matrix Pt in the initial coordinates (x)
         P=zeros(0,3);
         for k=1:size(Pt,1)
             P=[P ; pr+y1*Pt(k,1)+y2*Pt(k,2) ];
         end

         % assign transparency and color to faces
         %   (green for real faces,
         %   yellow for prescribed cut faces
         %   and red for self-acting cut faces)
         transp=1;
         c='g';
         if transparency
            transp= 0.4; 
         end
         if i>m
            transp=0.4; 
            if var
               c='y';
            else
               c='r';
            end
         end

         % draw i-support
         fill3(P(:,1),P(:,2),P(:,3),c,'FaceAlpha',transp);

      end

   end    

end