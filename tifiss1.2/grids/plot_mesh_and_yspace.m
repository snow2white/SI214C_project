function plot_mesh_and_yspace(evt,xy,evtY,xyY,subdivPar,iplot_elem,iplot_vtx)
%PLOT_MESH_AND_YSPACE plots mesh and its uniform (red/bisec3) refinement 
%
% input:
%              evt      element mapping matrix
%               xy      vertex coordinate vector  
%             evtY      element mapping matrix for Y-space
%              xyY      vertex coordinate vector for Y-space 
%        subdivPar      red or bisec3 refinement switch
%
%       iplot_elem      (optional) element's number plotting switch
%        iplot_vtx      (optional) vertices's number plotting switch
%
% Plots the mesh given by the coordinates and connections of the nodes 
% and the uniform refinement associated with the Y-space (either red
% or bisec3 subdivision)
%
% Example:
%    xy   = [-0.5 -0.5; 0.5 -0.5; 0.5 0.5; -0.5 0.5; 0.0 0.0];
%    evt  = [1 5 4; 2 5 1; 3 5 2; 4 5 3];
%    xyY  = [-0.25 -0.25; 0.25 -0.25; -0.25 0.25; 0.25 0.25; ...
%            -0.50  0.00; 0.00 -0.50;  0.50 0.00; 0.00 0.50];
%    evtY = [3 5 1; 1 6 2; 2 7 4; 4 8 3];
%    plot_mesh_and_yspace(evt,xy,evtY,xyY,2,1,1);
%
% See also PLOT_MESH, PLOTCOLOR_MARKEDELEM
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2017 A. Bespalov, L. Rocchi
 
  if nargin < 7
      iplot_vtx = 0;        % no vertices' numbers
      if nargin < 6
          iplot_elem = 0;   % no elements' numbers
          if nargin < 5
             error('Insufficient input number!');
          end  
      end 
  end
            
  nel   = size(evt,1);      % number of elements
  nvtx  = size(xy,1);       % total number of vertices (X space) 
  nvtxY = size(xyY,1);      % total number of midpoints (Y space)
  
% List of available colours
  redcol     = [255   0   0]./255; % red
  bisec3col  = [255 140   0]./255; % orange   
  darkBlue   = [0   141 226]./255;
  darkGreen  = [0.0 100 0.0]./255;
  lightGrenn = [0.0 153 0.0]./255;
  brownPeru  = [205 133  63]./255;
  chocolate  = [210 105  30]./255;
    
% -------------------------------------------------------------  
% Adjacency matrix
% -------------------------------------------------------------  
% This is a sparse adjacency matrix for the mesh given by evt matrix
  adj = sparse(nvtx,nvtx);
  evtt = [evt,evt(:,1)];
  for j = 1:3
      ncol1 = evtt(:,j);
      ncol2 = evtt(:,j+1);
      adj = adj + sparse(ncol1,ncol2,1,nvtx,nvtx);
  end
  
% -------------------------------------------------------------  
% Plot mesh
% -------------------------------------------------------------  
  figure;
  hold on;
  gplot(adj,xy,'b');

% Display element's numbers  
  if iplot_elem
      xl_v = zeros(nel,3); 
      yl_v = zeros(nel,3); 
      for ivtx = 1:3
          xl_v(:,ivtx) = xy(evt(:,ivtx),1); % x-coordinates of all nodes
          yl_v(:,ivtx) = xy(evt(:,ivtx),2); % y-coordinates of all nodes
      end
      % Element's centroid coordinates
      xyc(:,1) = sum(xl_v,2) / 3;
      xyc(:,2) = sum(yl_v,2) / 3;
      %
      % Write element's number
      elem = (1:nel)';
      elenum = int2str(elem);      
      text(xyc(:,1),xyc(:,2),elenum,'Color','blue','Fontsize',12);
  end

% Display vertices' numbers
  if iplot_vtx      
      vertices = (1:nvtx)';
      vtxnum = int2str(vertices);
      text(xy(:,1),xy(:,2),vtxnum,'Color','black','Fontsize',14);
  end   

% -------------------------------------------------------------------
% Plot Y space w.r.t. the uniform subdivision
% -------------------------------------------------------------------   
  vtxY = (1:nvtxY)';
  vtxnum = int2str(vtxY);
  if subdivPar == 1
      %
      % Red subdivision
      %
      % Write vertices' numbers
      %plot(xyY(:,1),xyY(:,2),'o','Color',redcol,'MarkerSize',6);
      text(xyY(:,1),xyY(:,2),vtxnum,'Color',redcol,'Fontsize',12);
      for elem = 1:nel
          Ycoord = xyY(evtY(elem,:),:);
          plot([Ycoord(2,1) Ycoord(1,1)], [Ycoord(2,2) Ycoord(1,2)],'-.','Color',redcol);    % Line Y_2 - Y_1
          plot([Ycoord(2,1) Ycoord(3,1)], [Ycoord(2,2) Ycoord(3,2)],'-.','Color',redcol);    % Line Y_2 - Y_3    
          plot([Ycoord(3,1) Ycoord(1,1)], [Ycoord(3,2) Ycoord(1,2)],'-.','Color',redcol);    % Line Y_3 - Y_1
      end
      
  else
      %
      % Bisec3 subdivision
      %
      % Write vertices' numbers
      %plot(xyY(:,1),xyY(:,2),'o','Color',bisec3col,'MarkerSize',6);
      text(xyY(:,1),xyY(:,2),vtxnum,'Color',bisec3col,'Fontsize',12);
      for elem = 1:nel
          Xcoord = xy(evt(elem,:),:);
          Ycoord = xyY(evtY(elem,:),:);
          plot([Ycoord(2,1) Ycoord(1,1)], [Ycoord(2,2) Ycoord(1,2)],'-.','Color',bisec3col); % Line Y_2 - Y_1
          plot([Ycoord(2,1) Xcoord(2,1)], [Ycoord(2,2) Xcoord(2,2)],'-.','Color',bisec3col); % Line Y_2 - X_2
          plot([Ycoord(2,1) Ycoord(3,1)], [Ycoord(2,2) Ycoord(3,2)],'-.','Color',bisec3col); % Line Y_2 - Y_3    
      end
  end

  axis square;
  set(gca,'XTick',[],'YTick',[]); 
  hold off;

end  % end function