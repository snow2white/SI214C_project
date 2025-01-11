function plotcolor_markedelem(M,evt,xy,iplot_elem,plot_title)
%PLOTCOLOR_MARKEDELEM plots marked elements of a mesh 
%   
% plotcolor_markedelem(M,evt,xy,iplot_elem,plot_title)
%
% input:
%               M      vector of marked elements
%             evt      element mapping matrix
%              xy      vertex coordinate vector 
%      iplot_elem      (optional) element's number plotting switch
%      plot_title      (optional) figure's title
%
% The function plots the mesh (given by evt and xy) with the filled coloured 
% marked elements (given by the vector M). 
% It is also possible to display the individual elements numbers as well as a 
% title.
%
% Default colour is rgb orange [1.0 0.8 0.4].
%
% ----------------------------------------------------------
% Example 1. Calling without any optional input:
%    plotcolor_markedelem(M,evt,xy)  
%
% Example 2. Calling with only elements' numbers:
%    plotcolor_markedelem(M,evt,xy,1);
%
% Example 3:
%    xy = [-0.5 -0.5; 0.5 -0.5; 0.5 0.5; -0.5 0.5; 0.0 0.0];
%    evt = [1 5 4; 2 5 1; 3 5 2; 4 5 3];
%    M = [1; 3];
%    plotcolor_markedelem(M,evt,xy,1,'Mesh example');
% ----------------------------------------------------------
%
% See also PLOT_MESH
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi

  if nargin < 5
      plot_title = '';     % no title          
      if nargin < 4
          iplot_elem = 0;  % no elements' numbers
          if nargin < 3
              error('Insufficient input number!');
          end
      end
  end
                
  nel  = size(evt,1);   % Number of elements
  nvtx = size(xy,1);    % Number of vertices
  
% List of available colours
  darkBlue   = [0   141 226]./255;
  darkGreen  = [0.0 100 0.0]./255;
  lightGrenn = [0.0 153 0.0]./255;
  brownPeru  = [205 133  63]./255;
  chocolate  = [210 105  30]./255;
  orange     = [255 140   0]./255;
  
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
  gplot(adj,xy,'b');
  hold on       

% Plot coloured marked elements
  LM = length(M);
  XY = zeros(3,2,LM);
  for i = 1:LM
      elem = M(i);
      vtxs = evt(elem,:);
      coordelem = xy(vtxs,:);
      XY(:,:,i) = coordelem;
      fill(XY(:,1,i),XY(:,2,i),orange);
  end
  
% Also display element's numbers  
  if iplot_elem
      xl_v = zeros(nel,3); 
      yl_v = zeros(nel,3); 
      for ivtx = 1:3
          xl_v(:,ivtx) = xy(evt(:,ivtx),1); % x-coordinates of all vertices
          yl_v(:,ivtx) = xy(evt(:,ivtx),2); % y-coordinates of all vertices
      end
      % Element's centroid coordinates
      xyc(:,1) = (1/3) * sum(xl_v,2);
      xyc(:,2) = (1/3) * sum(yl_v,2);
      %
      % Write element's number
      elem = (1:nel)';
      elenum = int2str(elem);      
      text(xyc(:,1),xyc(:,2),elenum,'Color','blue','Fontsize',12);
  end
  
  axis('square');
  title(plot_title,'Fontsize',16);
  hold off;
   
end % end function