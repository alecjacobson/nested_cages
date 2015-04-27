function Pf = peng2004_2d_point_polygon(x0,y0,P,steps,h,a,signed,vis)
  % PENG2004_2D_POINT_POLYGON moves a list of points with initial coordinates
  % (x0,y0) towards inside a polygon by
  % flowing towards the gradient direction of the p-distance defined in
  % "Interactively Modeling of Topologically Complex Geometric Detail"
  % Peng et a. [2004].
  %
  % Pf = peng2004_2d_point_polygon(x0,y0,a0,b0,a1,b1,steps,h,a)
  %
  % Input:
  %   x0:  #moving_points by 1 x-ccordinates of the initial moving points
  %   y0:  #moving_points by 1 y-ccordinates of the initial moving points
  %   P:  #vertices by 2 list of x-y coordinates of the polygon that defines
  %   the distance field
  %   steps:   number of steps to be taken in the gradient flow
  %   h:    time step
  %   a:    exponent that defines the p-distance in Peng et al. [2004]
  %   signed:   true or false to determine when to use singed distance
  %             (meaningful only if p is odd, not clear what it means when
  %              p is even though. Complex numbers show up)
  %   vis:   can be either 'vertices' to visualize vertices
  %          moving individually or 'polygon' to see the polygon defined
  %          by the vertices moving
  % Output:
  %   Pf:   #moving_points by 2 vector with final positions
  
  X = [x0 y0];
  
  % visualoization settings
  if (strcmp(vis,'vertices'))
      X_handle = plot(X(:,1),X(:,2),'.','MarkerSize',7);
      hold on
      plot([P(:,1);P(1,1)],[P(:,2);P(1,2)],'r','Linewidth',2)
      input('')
  elseif (strcmp(vis,'polygon'))
      X_handle = plot([X(:,1);X(1,1)],[X(:,2);X(1,2)],'Linewidth',1.5);
      hold on
      plot([P(:,1);P(1,1)],[P(:,2);P(1,2)],'r','Linewidth',1.5)
      input('')
  end

  for k=1:steps
      
      if (a==inf)
          [grad,~,~] = signed_distance_direction(X,P,[[(1:size(P,1))'] [(2:size(P,1))';1]]);
          grad = -grad;
      else
          
          grad_int = zeros(size(x0,1),2);
          int = zeros(size(x0,1),1);
          
          % loop over all segments
          for j=1:size(P,1)
              
              if j==size(P,1)
                  a0 = P(size(P,1),1); b0 = P(size(P,1),2); a1 = P(1,1); b1 = P(1,2);
              else
                  a0 = P(j,1); b0 = P(j,2); a1 = P(j+1,1); b1 = P(j+1,2);
              end
              
              [~,int_seg,grad_seg,~] = peng2004_2d_point_segment_integral_gradient(X,a0,b0,a1,b1,a);
              
              % sum to current gradient of integral and current integral
              grad_int = grad_int + grad_seg;
              int = int+int_seg;
              
          end
          
          % calculate the gradient of the distance
          grad = -[int.^(-1/a-1) int.^(-1/a-1)].*grad_int;
          % normalize (important! to avoid blowing up near the surface)
          grad = normalizerow(grad);
          
          % flip gradients accorindg to inside/outside
          if (signed)
              signs = -(2*inpolygon(X(:,1),X(:,2),P(:,1),P(:,2))-1);
              grad(:,1) = signs.*grad(:,1);
              grad(:,2) = signs.*grad(:,2);
          end
          
      end
      
      % evolve surface
      X = X+h*(-grad);
      
      % draw
      if (strcmp(vis,'vertices'))
          delete(X_handle)
          X_handle = plot(X(:,1),X(:,2),'.','MarkerSize',7);
          input('')
          drawnow
      elseif (strcmp(vis,'polygon'))
          delete(X_handle)
          X_handle = plot([X(:,1);X(1,1)],[X(:,2);X(1,2)],'Linewidth',1.5);
          input('')
          drawnow;
      end
      
  end
  
  Pf = X;
  
end