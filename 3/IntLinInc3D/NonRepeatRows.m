function [A] = NonRepeatRows(A)
%
% The function  NonRepeatRows  deletes repeated rows from the input matrix. 

   [m] = size(A,1);
   i=1;
   mi=m;

   while i < mi,
      p = [];
      
      for j=(i+1):mi
%         if isequal([A(i,:),A(j,:)]
         if max(abs(A(i,:)-A(j,:)))<1.e-12
            p = [p j];
         end
      end
   
      A(p,:) = [];
      mi=size(A,1);
      i=i+1;
   end       
      
end                
