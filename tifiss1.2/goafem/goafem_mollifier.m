function goafem_mollifier(dom_type,x0,y0,r,plotmoll)
%GOAFEM_MOLLIFIER computes the normalization constant of [PO99] mollifier
%
% input: 
%     dom_type     domain type
%      (x0,y0)     internal fixed point
%            r     radius of the ball B((x0,y0),r)
%     plotmoll     (optional) if =1 plot the mollifier over the domain
%
% The function implements the mollifier defined in [P099] for three different 
% domains: 
% - 1 unit square domain (0,1)^2
% - 2 L-shaped domain    (-1,1)^2 \ (-1,0]^2
% - 3 large crack domain (-1,1)^2 \ (-1,0)x{0} (crack on the left)
% In particular, the normalization constant is computed and saved in: 
% datafiles/constmollifier_po99.dat
% 
% To be used as: 
% - goafem_mollifier(dom_type,x0,y0,r)    % does not show the mollifier
% - goafem_mollifier(dom_type,x0,y0,r,1)  % shows the mollifier
%
% Reference:
% [P099] Prudhomme, Oden, On goal-oriented error estimation for elliptic 
% problems: application to the control of pointwise errors, Comput. Methods 
% Appl. Engrg., 176(1-4)313-331, 1999.
%
% Function(s) called: triquad
%
% See also GOAFEM_PO09_L2GOAL
%
%   TIFISS function: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

  if nargin < 5
      plotmoll = 0; % do not plot the mollifier
  end
    
% -----------------------------------------------------------------------------    
% (Not normalized) Mollifier
% -----------------------------------------------------------------------------    
                
% Support of the ball B((x0,y0);r)                
  supp = @(X,Y) sqrt( (x0 - X).^2  + (y0 - Y).^2 ) < r;  

% Exponential function 
  expfunc = @(X,Y) exp( -r^2 ./ ( r^2 - ((x0 - X).^2  + (y0 - Y).^2) ) );
            
% Not normalized mollifier (without C): expfunc inside the ball and zero outside
  moll = @(X,Y) expfunc(X,Y) .* supp(X,Y);
      
% -----------------------------------------------------------------------------    
% Split the domain in (four) triangles to handle integrals 
% -----------------------------------------------------------------------------    
  if dom_type == 1
      % Unit square domain (0,1)^2
      triangle1 = [1.0 0.0; 0.5 0.5; 0.0 0.0];
      triangle2 = [1.0 1.0; 0.5 0.5; 1.0 0.0];
      triangle3 = [0.0 1.0; 0.5 0.5; 1.0 1.0];
      triangle4 = [0.0 0.0; 0.5 0.5; 0.0 1.0];
  elseif dom_type == 2
      % L-shaped domain (-1,1)^2 \ (-1,0]^2
      triangle1 = [-1.0 0.0; 0.0 0.0; -1.0  1.0];
      triangle2 = [-1.0 1.0; 0.0 0.0;  1.0  1.0];
      triangle3 = [ 1.0 1.0; 0.0 0.0;  1.0 -1.0];
      triangle4 = [-1.0 1.0; 0.0 0.0;  0.0 -1.0];
  elseif dom_type == 3
      % Large crack domain (-1,1)^2 \ (-1,0)x{0} 
      load('crack_grid.mat','xP1','yP1','xP2','yP2','xf','yf'); 
      triangle1 = [-1.0 -1.0;  1.0 -1.0;  1.0  1.0];
      triangle2 = [-1.0  1.0;  0.0  0.0;  1.0  1.0];
      triangle3 = [  xf   yf; -1.0  1.0;  xP2  yP2];
      triangle4 = [  xf   yf;  xP1  yP1; -1.0 -1.0];
  end
 
% -----------------------------------------------------------------------------    
% Integral over the domain (the four triangles). The normalization constant 
% C satisfies \int_D C moll(x)dx = 1. Hence: 
% C = (\int_D moll(x)dx)^{-1} = (\int_B((x0,y0);r) moll(x)dx)^{-1} 
% -----------------------------------------------------------------------------    

% Gaussian nodes in each triangle  
  gnodes = 150;

% First triangle
  [X1,Y1,Wx1,Wy1] = triquad(gnodes,triangle1);
  mollev = moll(X1,Y1); mollev(isnan(mollev)) = 0.0;
  int1 = Wx1' * mollev * Wy1;
  
% Second triangle
  [X2,Y2,Wx2,Wy2] = triquad(gnodes,triangle2);
  mollev = moll(X2,Y2); mollev(isnan(mollev)) = 0.0;
  int2 = Wx2' * mollev * Wy2;
     
% Third triangle
  [X3,Y3,Wx3,Wy3] = triquad(gnodes,triangle3);
  mollev = moll(X3,Y3); mollev(isnan(mollev)) = 0.0;
  int3 = Wx3' * mollev * Wy3;

% Fourth triangle
  [X4,Y4,Wx4,Wy4] = triquad(gnodes,triangle4);
  mollev = moll(X4,Y4); mollev(isnan(mollev)) = 0.0;
  int4 = Wx4' * mollev * Wy4;
  
% Sum of the four contributing integrals 
  sumint = int1 + int2 + int3 + int4;
  
% Normalization constant
  C = (sumint)^(-1);

% Show the mollifier?
  if plotmoll%==1     
      plotmollifier(dom_type,C,moll);
  end 
  
% Save data
  gohome; cd datafiles;
  save constmollifier_po99.mat C x0 y0 r;
  
end % end function


% -----------------------------------------------------------------------------
% Child function
% -----------------------------------------------------------------------------
function plotmollifier(dom_type,C,moll)

% Load data
  if dom_type == 1
      % taken from gohome/datafiles
      load('square_grid.mat','xy','mv','bound','mbound'); 
      [evt,eboundt] = p1grid(xy,mv,bound,mbound,0);
  elseif dom_type == 2
      % taken from gohome/datafiles
      load('ell_grid.mat','xy','mv','bound','mbound'); 
      [evt,eboundt] = p1grid(xy,mv,bound,mbound,0);
  elseif dom_type == 3 
      % taken from gohome/datafiles
      load('crack_grid.mat','xy','evt','bound','eboundt');
  end
% Uniform refinements for 
  for i = 1:3
      [xyc,evtc,~,~,~] = uniform_refinement(xy,evt,bound,eboundt,2);
  end
  % (Normalized) mollifier
  normoll = @(X,Y) C * moll(X,Y);
  figure;
  trimesh(evtc,xyc(:,1),xyc(:,2),normoll(xyc(:,1),xyc(:,2)));
  axis square;

end % end child function