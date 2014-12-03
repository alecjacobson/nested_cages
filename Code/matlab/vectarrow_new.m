function h = vectarrow_new(P,V)
% VECTARROW_NEW plot a list of points and vectors
%
% vectarrow(P,V)
% 
% Inputs
%   P  list of points #x3 or #x2
%   V  list of vectors #x3 or #x2
% Output 
%   h  3D plot handle

hold on
% initialize handle with something
if size(P,2)==3
    h = plot3(1000,1000,1000);
else
    h = plot(1000,1000);
end

xi = [];
yi = [];
if size(P,2)==3
    zi = [];
end

for index = 1:size(P,1)
    
    p0 = P(index,:);
    p1 = P(index,:) + V(index,:);
    
    if max(size(p0))==3
        if max(size(p1))==3
            x0 = p0(1);
            y0 = p0(2);
            z0 = p0(3);
            x1 = p1(1);
            y1 = p1(2);
            z1 = p1(3);
            
            p = p1-p0;
            alpha = 0.1;  % Size of arrow head relative to the length of the vector
            beta = 0.1;  % Width of the base of the arrow head relative to the length
            
            % Nan was added to avoid connection between arrowa when
            % plotting (cool trick!)
            hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps)); nan];
            hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps)); nan];
            hw = [z1-alpha*p(3);z1;z1-alpha*p(3); nan];
            
            xi = [xi;x0;x1;hu(:)];
            yi = [yi;y0;y1;hv(:)];
            zi = [zi;z0;z1;hw(:)];
            
        else
            error('p0 and p1 must have the same dimension')
        end
    elseif max(size(p0))==2
        if max(size(p1))==2
            
            x0 = p0(1);
            y0 = p0(2);
            x1 = p1(1);
            y1 = p1(2);
            
            p = p1-p0;
            alpha = 0.1;  % Size of arrow head relative to the length of the vector
            beta = 0.1;  % Width of the base of the arrow head relative to the length
            
            % Nan was added to avoid connection between arrowa when
            % plotting (cool trick!)
            hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps)); nan];
            hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps)); nan];
            
            xi = [xi;x0;x1;hu(:)];
            yi = [yi;y0;y1;hv(:)];
            
        else
            error('p0 and p1 must have the same dimension')
        end
    else
        error('this function only accepts 2D or 3D vector')
    end
    hold off;
end
if size(P,2)==3
    set(h,'Xdata',xi,'Ydata',yi,'Zdata',zi);
else
    set(h,'Xdata',xi,'Ydata',yi);
end