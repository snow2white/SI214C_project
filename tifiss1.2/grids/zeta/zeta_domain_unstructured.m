%ZETA_DOMAIN_UNSTRUCTURED Z-shaped unstructured grid generator
%
% The function takes the following inputs from the user:
%
%                     h     the mesh size
%    mesh_density_ratio     relative mesh density close the corners
%                 slope     diagonal edge's slope
%
%    Calls DISTMESH function: distmesh2d
%
% Grid defining data is saved to the file: zeta_grid_unstruct.mat
%
%   TIFISS scriptfile: AB; 07 June 2017.
% Copyright (c) 2017 A. Bespalov, L. Rocchi
%
%
% NOTE 
% -------------------------------------------------------------------------
% 1. A smaller 'mesh_density_ratio' leads to finer grids close corners
% 2. 'slope' represents the slope of the diagonal line (y=(slope)x) of the
%    Z-domain:
%             slope >> 1 (~60) creates a L-shaped domain
%             slope == 1       creates a pure Zeta domain
%             slope << 1       creates a crack domain
%             -----------------------------------------
%            (slope < 0        creates a square domain)
%
% To reproduce a crack domain insert slope~0.02 (less or equal): in this
% case no mesh generator is called. The crack domain mesh is already stored
% for slope=0.02 and saved in crack_mesh.mat file in /grids folder


  global slope                % Used by the function trapezoid_fd
  global mesh_density_ratio   % Used by the function zeta_fh
  
% Instruction for the user
% -------------------------------------------------------------------------
  fprintf('You are now asked for mesh size, mesh density near corner (0,0), and\n');
  fprintf('the diagonal edge''s slope defining the domain:\n');
  fprintf('   slope >> 1     nearly a L-shaped domain (~60)\n');
  fprintf('   slope == 1     pure Z-shaped domain\n');
  fprintf('   slope << 1     nearly a crack domain (~0.02)\n');
  
% Set distmesh2d parameters
% -------------------------------------------------------------------------
  slope = default('Diagonal edge''s slope (default 0.7)',0.7);
  if slope <= 0.02
      fprintf('\nGenerating triangulation for crack-domain...');
  else
      h = default('Target mesh size (default 0.2)',0.2);
      mesh_density_ratio = default('Mesh density near corner (default 0.4)',0.4);
      fprintf('\nGenerating triangulation...');
  end

% Fix some boundary points for better precision of the domain
% -------------------------------------------------------------------------
  if slope > 1
      % 'Internal' Zeta-domain
      pfix = [1.0, -1.0; 1.0, 0.0; 1.0, 1.0; 0.0, 1.0; -1.0, 1.0; ...
                -1.0, 0.0; -0.5, 0.0; 0.0, 0.0; -1/slope, -1.0; 0.0, -1.0]; 
      
  elseif slope == 1  
      % Pure Zeta-domain
      pfix = [1.0, -1.0; 1.0, 0.0; 1.0, 1.0; 0.0, 1.0; -1.0, 1.0; ...
                    -1.0, 0.0; -0.5, 0.0; 0.0, 0.0; -1.0, -1.0; 0.0, -1.0];
                
  elseif slope < 1
      if slope <= 0.02
          % Load unstructured crack domain
          load crack_mesh_unstruct;
      else
          % Zeta domain
          pfix = [1.0, -1.0; 1.0, 0.0; 1.0, 1.0; 0.0, 1.0; -1.0, 1.0; -1.0, 0.0; ...
                  -0.5, 0.0; 0.0, 0.0; -1.0, -slope; -1.0, -1.0; 0.0, -1.0];
      end
  end

% Calling mesh generator
% -------------------------------------------------------------------------
  if slope > 0.02
      % Call distmesh2d
      [xy,evt] = distmesh2d(@zeta_fd, @zeta_fh, h, [-1,-1;1,1],pfix);
      fprintf('done\n');
  end

% -------------------------------------------------------------------------  
  
% Get boundary nodes and edges
  [eboundt,bound] = find_boundary_elements(evt);
  
% Check domain
% -------------------------------------------------------------------------
% If the corner node (0.0) is not a boundary node the mesh has not been
% generated correctly. The reason is a too large mesh size or few mesh 
% density near corners
  xybnd = xy(bound,:);
  ii = find(xybnd(:,1)==0.0);
  if isempty(find( xybnd(ii,2)==0.0, 1))
      fprintf('WARNING: Z-domain not correctly generated!');
      fprintf('\nDecrease either mesh size or mesh density near corners...\n');
  end
  
% Check the mesh regularity
  [xy,evt,eboundt,bound] = mesh_regularity_check(xy,evt,eboundt,bound);
  
% reset boundary nodes: original lines from ell_domain function
% -------------------------------------------------------------------------
% NOTE: this does not work for ZETA-domain, it works for L-shaped domain
% -------------------------------------------------------------------------
% neboundt=length(eboundt(:,1));
% for i = 1:neboundt
%     iel = eboundt(i,1);
%     nodes = evt(iel,:);
%     if eboundt(i,2)==1, nodeb=evt(iel,[2,3]); xybd=xy(nodeb,:); end
%     if eboundt(i,2)==2, nodeb=evt(iel,[3,1]); xybd=xy(nodeb,:); end
%     if eboundt(i,2)==3, nodeb=evt(iel,[1,2]); xybd=xy(nodeb,:); end
%     h=norm(xybd(1,:)-xybd(2,:),2);
%     % vertical bound or horizontal bound
%     if abs(xybd(1,2)-xybd(2,2))>abs(xybd(1,1)-xybd(2,1))
%        if xybd(1,1)<-0.5, xy(nodeb,1)=-1;
%        elseif xybd(1,1)<0.5, xy(nodeb,1)=0;
%        else xy(nodeb,1)=1;
%        end
%     else
%        if xybd(1,2)<-0.5, xy(nodeb,2)=-1;
%        elseif xybd(1,2)<0.5, xy(nodeb,2)=0;
%        else  xy(nodeb,2)=1;
%        end
%     end        
% end

% Reset boundary nodes
% -------------------------------------------------------------------------
% Some x or y coordinates on the boundaries need to be rewritten with the 
% correct value (\pm 1, or 0) as distmesh2d could move a little bit the 
% points. Usually the boundary points along the square are computed well. 
% The movements could occour along the horizontal negative bound with y~=0
  xyhor = xy(xy(:,1)<=0.0,:);
  xy(abs(xyhor(:,2)-0.0)<1e-05,2) = 0.0;
  
% Save data
  gohome
  cd datafiles
  save zeta_grid_unstruct.mat xy evt bound eboundt 
