function [pr,At,bt,y1,y2]=ChangeVariable(A,b,i)
%
% The function ChangeVariable transforms the system of linear relations
% A x relations b, given in the initial coordinates (0,x1,x2,x3),
% into the system  At y relations* bt, written in intrinsic coordinates 
% (pr,y1,y2) of the plane A(i,:)x=b(i).
% (The vector 'relations' is a vector with components '=', '<', '>', '<=' or '>='.
% The vector 'relations*' is the vector 'relations' without the i-th component.)
%
% Written by Irene A. Sharaya (Institute of Computational Technologies, Novosibirsk)
% Date:  September 15, 2012

   % as the origin  pr  of the intrinsic coordinates system, 
   % we take the projection of the origin of the initial coordinates system
   % onto the plane A(i,:)x=b(i) 
   pr= A(i,:)*b(i)/(A(i,:)*A(i,:)');

   % initialization
   At=zeros(0,2);
   bt=zeros(0,1);
   y1=zeros(1,3);
   y2=zeros(1,3);

   % determining the orts y1 and y2

     % y1 
     % we take the longest projection of the vector A(i,:) onto the coordinate planes
     % and rotate it by 90 degrees within the projection plane 
     [c,l]=min(abs(A(i,:))); 
     l1=mod(l,3)+1;
     l2=mod(l1,3)+1;
     y1(l1)=A(i,l2);
     y1(l2)=-A(i,l1);

     % y2 
     % computing cross product of the vectors A(i,:) and y1
     y2=cross(A(i,:),y1);

   % calculation of At and bt
   A(i,:)=[];
   b(i,:)=[];
   At=[A*y1' A*y2'];
   bt=b-A*pr';

end
 