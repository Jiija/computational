function [S] = BoundaryIntervals(A,b)
%
% The function BoundaryIntervals computes the matrix S of boundary 
% intervals for the system of linear inequalities A*x >= b, 
% where A is a real matrix with 2 columns and m nonzero rows (m >= 1),  
% b is a real m-vector, x is a real 2-vector of unknowns.
 
   m = size(A,1);
   S = zeros(0,5);

   for i=1:m
      si=true;     
      q=[-Inf Inf];  
      dotx= A(i,:)*b(i)/(A(i,:)*A(i,:)');  
      p=[-A(i,2) A(i,1)]; 

      for k=1:m        
         Akx=A(k,:)*dotx';  
         c=A(k,:)*p';       
         switch sign(c)
            case -1
               q(2)=min(q(2),(b(k)-Akx)/c);
            case 1
               q(1)=max(q(1),(b(k)-Akx)/c);
            otherwise
               if single(Akx)<single(b(k))
                  if single(A(k,:)*A(i,:)')>0
                     si=false; 
                     break 
                  else 
                     return 
                  end
               end
         end
      end

      if single(q(1))>single(q(2))
         si=false;
      end 

      if si
      S=[ S ; dotx+p*q(1) , dotx+p*q(2) , i];
      end

   end
               
end