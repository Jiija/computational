function [V,cfinite]=AddV(V,cfinite,A,b,m,cnmty)
%
% The function AddV adds, to the list V, vertices of the solution set H
% of the linear inequalities system A x >= b. 
% It is known beforehand that H lies within a single orthant.
%
% Written by Irene A. Sharaya (Institute of Computational Technologies, Novosibirsk)
% Date:  September 15, 2012

   if ~cnmty % there are no points of the set H within the orthant 
      return;
   end

   for i=1:m

      % switch from the initial coordinates system (0,x1,x2,x3)
      % to the intrinsic coordinates system (pr,y1,y2) of plane i
      % (more exactly, of the plane A(i,:)x=b(i))
      [pr,At,bt,y1,y2]=ChangeVariable(A,b,i);

      % delete zero rows from  At
      [At,bt,cnmtyt]=ClearZeroRows(At,bt);

      if ~cnmtyt
         continue; % the plane i does not intersect the set H
      end

      [St]=BoundaryIntervals(At,bt);

      if size(St,1)==0 % 
         continue; % the plane i does not intersect the set H
      end

      % gather the vertices from the plane i to form the matrix Vt

        Vt=St(:,1:2);  

        xbf=isfinite(Vt(:,1));
        ybf=isfinite(Vt(:,2));
        binf=~min(xbf,ybf); % the inequalities St(binf,5) give 
                            % boundary intervals starting in the infinity 

        if max(binf)
           cfinite=0;
        end

        Vt(binf,:)=[]; 

      % Vtt is the matrix Vt in the original coordinates (x)
      Vtt=zeros(0,3); 
      for k=1:size(Vt,1)
          Vtt=[Vtt ; pr+y1*Vt(k,1)+y2*Vt(k,2) ];
      end


      V=[V;Vtt];

   end

end