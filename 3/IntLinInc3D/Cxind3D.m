function [V]=Cxind3D(uC,oC,ud,od,OrientPoints,transparency,varargin) 
%
% The function Cxind3D visualizes the set of formal solutions 
% to the interval inclusion [uC,oC] x \subseteq [ud,od] 
% in Kaucher arithmetic.
%
% Input arguments:
%   uC, oC  are real matrices with  m  rows and  3 columns
%     (the left and right endpoints of the interval matrix C=[uC,oC] respectively);
%   ud and od  are vectors of the length m, 
%     whose components may be real numbers or -Inf or +Inf 
%     (the left and right endpoints of the interval vector d=[ud,od] respectively),
%   OrientPoints is a parameter that controls drawing of the orientation points; 
%     if it has the value 0, the program does not plot these points; 
%     if it has the value 1, the orientation points are plotted; 
%   transparency is a parameter responsible for the transparency 
%     of the real faces, its values may be 0 or 1;   
%     0 means the absence of the transparency, 
%     1 means that the real faces are transparent.
%
% Optional input arguments:
%   varargin is a prescribed cut box, 
%     it is inputted as 6 numbers, xb, xe, yb, ye, zb, ze, where 
%     xb,yb,zb are coordinates of the lower endpoint of the box, 
%     xe,ye,ze are coordinates of the upper endpoint of the box.
%
% Output arguments:
%   V is a matrix of the orientation points.
%     In the package, the orientation points are understood as vertices 
%     of intersections of the solution set with orthants.
%
% Written by Irene A. Sharaya (Institute of Computational Technologies, Novosibirsk)
% Date:  September 15, 2012

   % In orthant  k,  the inclusion [uC,oC] x \subseteq [ud,od] 
   % is equivalent to the system of linear inequalities Ak(:,:) x >= bk 
   % where Ak è bk are as follows: 
   Appp(:,:)=[eye(3); uC; -oC ];
   Appm(:,:)=[  1 0 0; 0  1 0; 0 0 -1; uC(:,1) uC(:,2) oC(:,3); -oC(:,1) -oC(:,2) -uC(:,3)];
   Apmp(:,:)=[  1 0 0; 0 -1 0; 0 0  1; uC(:,1) oC(:,2) uC(:,3); -oC(:,1) -uC(:,2) -oC(:,3)];
   Ampp(:,:)=[ -1 0 0; 0  1 0; 0 0  1; oC(:,1) uC(:,2) uC(:,3); -uC(:,1) -oC(:,2) -oC(:,3)];
   Apmm(:,:)=[  1 0 0; 0 -1 0; 0 0 -1; uC(:,1) oC(:,2) oC(:,3); -oC(:,1) -uC(:,2) -uC(:,3)];
   Ampm(:,:)=[ -1 0 0; 0  1 0; 0 0 -1; oC(:,1) uC(:,2) oC(:,3); -uC(:,1) -oC(:,2) -uC(:,3)];
   Ammp(:,:)=[ -1 0 0; 0 -1 0; 0 0  1; oC(:,1) oC(:,2) uC(:,3); -uC(:,1) -uC(:,2) -oC(:,3)];
   Ammm(:,:)=[-eye(3); oC; -uC ];
   bppp=[ 0; 0; 0; ud; -od ];
   bppm=bppp;
   bpmp=bppp;
   bmpp=bppp;
   bpmm=bppp;
   bmpm=bppp;
   bmmp=bppp;
   bmmm=bppp;

   % For each orthant, we revise rows in the respective system of linear inequalities.
   [Appp,bppp,cnmtyppp]=ClearRows(Appp,bppp);
   mppp=size(bppp,1);
   [Appm,bppm,cnmtyppm]=ClearRows(Appm,bppm);
   mppm=size(bppp,1);
   [Apmp,bpmp,cnmtypmp]=ClearRows(Apmp,bpmp);
   mpmp=size(bpmp,1);
   [Ampp,bmpp,cnmtympp]=ClearRows(Ampp,bmpp);
   mmpp=size(bmpp,1);
   [Apmm,bpmm,cnmtypmm]=ClearRows(Apmm,bpmm);
   mpmm=size(bpmm,1);
   [Ampm,bmpm,cnmtympm]=ClearRows(Ampm,bmpm);
   mmpm=size(bmpm,1);
   [Ammp,bmmp,cnmtymmp]=ClearRows(Ammp,bmmp);
   mmmp=size(bmmp,1);
   [Ammm,bmmm,cnmtymmm]=ClearRows(Ammm,bmmm);
   mmmm=size(bmmm,1);

   % initialization of  V  and  cfinite 
   V=zeros(0,3);
   cfinite=1; % we suppose that the solution set is bounded 

 
   % If the intersection of the solution set with an orthant is not empty,
   % we add vertices of the intersection into the matrix V  of the orientation points.
   % If the intersection is unbounded, we set  cfinite=0.
   % We do that for each orthant.
   [V,cfinite]=AddV(V,cfinite,Appp,bppp,mppp,cnmtyppp);
   [V,cfinite]=AddV(V,cfinite,Appm,bppm,mppm,cnmtyppm);
   [V,cfinite]=AddV(V,cfinite,Apmp,bpmp,mpmp,cnmtypmp);
   [V,cfinite]=AddV(V,cfinite,Ampp,bmpp,mmpp,cnmtympp);
   [V,cfinite]=AddV(V,cfinite,Apmm,bpmm,mpmm,cnmtypmm);
   [V,cfinite]=AddV(V,cfinite,Ampm,bmpm,mmpm,cnmtympm);
   [V,cfinite]=AddV(V,cfinite,Ammp,bmmp,mmmp,cnmtymmp);
   [V,cfinite]=AddV(V,cfinite,Ammm,bmmm,mmmm,cnmtymmm);

   V = NonRepeatRows(V);
   numV=size(V,1);
   fprintf('Number of orientation points = %d\n\n',numV) 
   if numV==0
      disp('The solution set is empty');
      return;
   end

   % We calculate the cut-drawing box W.
   var=(length(varargin)==6);
   if var % If a user inputs a prescribed cut box 
          % as the optional input argument  varargin,
          % then W is set equal to the prescribed cut box.
      W=[varargin{1} varargin{3} varargin{5}; varargin{2} varargin{4} varargin{6}];
      disp('The cut box is prescribed');
%      if ~( W(1,:)<=W(2,:) )
%         disp('The cut box is prescribed incorrectly');
%         return;
%      end
   else % Otherwise, we choose the cut-drawing box W 
        % so that W contains all the orientation points.
      [W]=ChooseDrawingBox(V,cfinite);
   end

   % If the solution set is unbounded or the cut box is prescribed, we inset 
   % six cut inequalities into the system Ak(:,:)x>=bk for each orthant k.
   if ~cfinite | var
      Appp=[Appp; eye(3); -eye(3) ];
      bppp=[bppp; W(1,:)';-W(2,:)'];
      Appm=[Appm; eye(3); -eye(3) ];
      bppm=[bppm; W(1,:)';-W(2,:)'];
      Apmp=[Apmp; eye(3); -eye(3) ];
      bpmp=[bpmp; W(1,:)';-W(2,:)'];
      Ampp=[Ampp; eye(3); -eye(3) ];
      bmpp=[bmpp; W(1,:)';-W(2,:)'];
      Apmm=[Apmm; eye(3); -eye(3) ];
      bpmm=[bpmm; W(1,:)';-W(2,:)'];
      Ampm=[Ampm; eye(3); -eye(3) ];
      bmpm=[bmpm; W(1,:)';-W(2,:)'];
      Ammp=[Ammp; eye(3); -eye(3) ];
      bmmp=[bmmp; W(1,:)';-W(2,:)'];
      Ammm=[Ammm; eye(3); -eye(3) ];
      bmmm=[bmmm; W(1,:)';-W(2,:)'];
   end


   % Drawing the solution set in all orthants, one after another.

%     figure    % each start of the package generates a new figure for drawing
     figure(1) % the package draws in figure 1 for all starts 
     clf(1);  % clearing figure 1 
              % so as the new picture does not overlap the previous one
     hold on
     axis vis3d; % auto, on, normal, fill, image, tight, equal
     axis([W(1,1),W(2,1),W(1,2),W(2,2),W(1,3),W(2,3)]);
     view(3);
     grid on;
     xlabel('x1');
     ylabel('x2');
     zlabel('x3');

     % If OrientPoints is equal to 1, we plot the orientation points.
     if OrientPoints 
        if var
           nVinW = [V(:,1)>=W(1,1) & V(:,1)<=W(2,1)... 
                  & V(:,2)>=W(1,2) & V(:,2)<=W(2,2)...
                  & V(:,3)>=W(1,3) & V(:,3)<=W(2,3)];
           VinW = V(nVinW,:);
           scatter3(VinW(:,1),VinW(:,2),VinW(:,3),15,'k','filled');
        else
           scatter3(V(:,1),V(:,2),V(:,3),15,'k','filled');
        end
     end

     DrawHedrons(Appp,bppp,mppp,transparency,var);
     DrawHedrons(Appm,bppm,mppm,transparency,var);
     DrawHedrons(Apmp,bpmp,mpmp,transparency,var);
     DrawHedrons(Ampp,bmpp,mmpp,transparency,var);
     DrawHedrons(Apmm,bpmm,mpmm,transparency,var);
     DrawHedrons(Ampm,bmpm,mmpm,transparency,var);
     DrawHedrons(Ammp,bmmp,mmmp,transparency,var);
     DrawHedrons(Ammm,bmmm,mmmm,transparency,var);

     hold off

end                
