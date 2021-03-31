function [A,b,cnmty] = ClearRows(A,b)
%
% The function ClearRows revises rows in the system of linear 
% inequalities A*x >= b, where 
%   A is a rectangular real matrix, 
%   b is a vector from the extended real axis (R with +Inf and -Inf).
% Output arguments are
%   an equivalent system (a new matrix A and a new vector b)
%   and a parameter  cnmty  (Control of NoneMpTiness). 
%   The parameter cnmty may take values 
%     'false' (if the solution set is obviously empty) and
%     'true' (if we do not know whether the solution set is empty or not).
%   If the output argument cnmty equals 'true', then 
%   the output matrix A does not have zero rows and 
%   the output vector b does not have components +Inf and -Inf.

   [m,n] = size(A);
   cnmty = true;
   p = [];

   for i=1:m

      if b(i)==Inf
         cnmty = false;
         return;
      end
      if b(i)==-Inf
         p = [p i];
         continue;
      end

%      if single(A(i,:)) == zeros(1,n)
      if max(abs(A(i,:)))<1.e-8
         p = [p i];
%         if 0 < b(i)
         if b(i)>1e-12
            cnmty = false;
            return;
         end
      end
   end
   
   A(p,:) = [];
   b(p,:) = [];
               
end                
