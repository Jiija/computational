function [A,b,cnmty] = ClearZeroRows(A,b)
%
% The function ClearZeroRows 
% revises rows in the system of linear inequalities A*x >= b, where 
%   A is a rectangular real matrix, 
%   b is a real vector, 
%   x is a real vector of unknowns,
% and deletes rows that correspond to zero rows of the matrix A.
% 
% Output arguments are
%   a new system of inequalities (a new matrix A and a new vector b)
%   and a parameter  cnmty  (Control of NoneMpTiness). 
%   The parameter cnmty may take values
%     'false' (if the solution set to the initial system is obviously empty) and
%     'true' (if we do not know whether this solution set is empty or not).
%   If  cnmty is 'true', then 
%     1) the new system is equivalent to the initial one,
%     2) the absence of rows in the new system means that
%        the solution set of the initial system coincides with the whole space
%        and therefore its boundary is empty.

   [m,n] = size(A);
   cnmty = true;
   p = [];

   for i=1:m
%      if single(A(i,:)) == zeros(1,n)
      if max(abs(A(i,:)))<1.e-8
         p = [p i];
%         if 0 < b(i)
         if b(i)>1e-12
            cnmty = false;
         end
      end
   end
   
   A(p,:) = [];
   b(p,:) = [];
               
end                
