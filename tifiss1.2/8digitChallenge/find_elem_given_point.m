function [el] = find_elem_given_point(xp,yp,xy,evt)
%FIND_ELEM_GIVEN_POINT  identifies an element that contains given point
%   [el] = find_elem_given_point(xp,yp,xy,evt);
%   input
%          xp         x-coordinate of the point
%          yp         y-coordinate of the point
%          xy         vertex coordinate vector  
%          evt        element mapping matrix
%    TIFISS function: AB; 15 August 2019.
% Copyright (c) 2018 Alex Bespalov, Leonardo Rocchi

  nvtx = size(xy,1);    % Number of vertices

  xxp = repmat(xp,nvtx,1);
  yyp = repmat(yp,nvtx,1);
  distances = sqrt( ( xxp - xy(:,1) ).^2 + ( yyp - xy(:,2) ).^2 ); 
  
% find the mesh-point closer to the input point (xp,yp)  
  [~,vtx] = min( distances );
  
% Note that if ~isequal(nvtx,nnz(distances)), 
% then the input point (xp,yp) is a mesh-node, i.e. it is 
% exactly the clvtx-th vertex of the mesh
  
 
% get the elements sharing the vtx-th mesh-point   
  elems = find( sum( ismember(evt, vtx) , 2 ) );


% Check which element of elems the input node (xp,yp) belongs to 
% Using barycentric coordinates for each element in elems:
%
%   xp = a*x1 + b*x2 + c*x3
%   yp = a*y1 + b*y2 + c*y3
%   1  = a + b + c  
%
% where (x1,y1), (x2,y2), and (x3,y3) are the coordinates of the (local) 
% 1st, 2nd, and 3rd vertex of a triangle in elems;
%
% Solution:
%   det := ((y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3));
%   a = ((y2 - y3)*(x - x3) + (x3 - x2)*(y - y3)) / det;
%   b = ((y3 - y1)*(x - x3) + (x1 - x3)*(y - y3)) / det;
%   c = 1 - a - b;
%
% * If 0<=a<=1, 0<=b<=1, and 0<=c<=1 then (xp,yp) belongs to the element.
% Note that this holds also if (xp,yp) belongs to an edge of the element as
% well as if it is a vertex of the element

% We take account of this last possibility:
  answers = zeros(length(elems),1); 

  for elk = 1:length(elems)
      %
      coordel = xy( evt( elems(elk) ,:), : );
      xx = coordel(:,1);
      yy = coordel(:,2);
      %
      % solving the associated linear system
      det = (yy(2) - yy(3))*(xx(1) - xx(3)) + (xx(3) - xx(2))*(yy(1) - yy(3));
      a = ( (yy(2) - yy(3))*(xp - xx(3)) + (xx(3) - xx(2))*(yp - yy(3)) ) / det;
      b = ( (yy(3) - yy(1))*(xp - xx(3)) + (xx(1) - xx(3))*(yp - yy(3)) ) / det;
      c = 1 - a - b;
      % 
      if ( 0.0 <= a && a <= 1.0 ) && ( 0.0 <= b && b <= 1.0 ) && ( 0.0 <= c && c <= 1.0 )
          % the point (xp,yp) belong to the elems(elk)-th triangle 
          answers(elk) = 1;
      end    
  end

% find the element which the input node (xp,yp) belongs to
  el = elems( answers==1 );
  
% if length(el) == 1, the input node belong to the el-th element  
  if length(el) > 1
      % = 2, if the input node belong to the edge shared by the two elements in el
      % = 4, if input node is a mesh-point, i.e. it belong to all elements which share it
      el = el(1);
  end

end % end function

