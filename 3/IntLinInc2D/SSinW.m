% The procedure  SSinW  draws the solution set.

   figure(1) 
   clf(1);   
%% or
%  figure   % if you want to save all pictures

   hold on % start drawing

      axis([W(1,1),W(2,1),W(1,2),W(2,2)])  
      if bounded
         axis equal
      end

      % grid on

      % dotted lines for axis
      %% abscissa axis
      if (0>=W(1,2)) & (0<=W(2,2))
      plot(W(:,1),[0,0],':b','LineWidth',1)
      end
      %% ordinate axis
      if (0>=W(1,1)) & (0<=W(2,1))
      plot([0,0],W(:,2),':b','LineWidth',1)
      end

      % draw the solution set by orthants
      % (interior is green, boundary is black, vertices are small circles)
      fill(P1(:,1),P1(:,2),'g') 
      plot(P1(:,1),P1(:,2),'o','MarkerFaceColor','k','MarkerSize',3)
      fill(P2(:,1),P2(:,2),'g') 
      plot(P2(:,1),P2(:,2),'o','MarkerFaceColor','k','MarkerSize',3)
      fill(P3(:,1),P3(:,2),'g') 
      plot(P3(:,1),P3(:,2),'o','MarkerFaceColor','k','MarkerSize',3)
      fill(P4(:,1),P4(:,2),'g') 
      plot(P4(:,1),P4(:,2),'o','MarkerFaceColor','k','MarkerSize',3)

   hold off % finish drawing
                