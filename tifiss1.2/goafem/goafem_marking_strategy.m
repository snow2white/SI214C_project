function [M] = goafem_marking_strategy(errsprimal,errsdual,strategy,threshold)
%GOAFEM_MARKING_STRATEGY Doerfler marking of elements/edges in goal-oriented framework 
%
% [M] = goafem_marking_strategy(errsprimal,errsdual,strategy,threshold)
%
% input: 
%     errsprimal    vector of element/edge indicators (primal solution)
%       errsdual    vector of element/edge indicators (dual solution)
%       strategy    marking strategy type (based on Doerfler)
%      threshold    threshold marking parameter 
%
% output: 
%              M    set of marked elements/edges
%
% Marking (of either elements or edges) uses Doerfler marking;
% The input 'strategy' allows to choose among six different strategies:
% - Feischl, Praetorius, Van Der Zee [FPZ16]
% - Mommer, Stevenson                [MS09]
% - Holst, Pollock                   [HP16]
% - Becker, Estecahandy, Trujillo    [BET11]
% - Only primal problem 
% - Only dual problem
% The first four strategies join the two sets of marked elements/edges of the 
% primal and the dual problem;
% The fifth strategy marks elements/edges associated with the primal solution;
% The sixth strategy marks elements/edges associated with the dual solution;
%
% References: 
% [FPZ16] Feischl, Praetorius, van der Zee, An abstract analysis of optimal 
% goal-oriented adaptivity, SIAM J. Numer. Anal., 54(3)1423-1448, 2016;
%
% [MS09] Mommer, Stevenson, A goal-oriented finite element method with
% convergence rates, SIAM J. Numer. Anal., 47(2)861-866, 2009;
%
% [HP16] Holst, Pollock, Convergence of goal-oriented adaptive finite element 
% methods for nonsymmetric problems, Numer. Meth. Part. D. E., 32(2)479-509, 2016;
%
% [BET11] Becker, Estecahandy, Trujillo, Weighted marking for goal-oriented 
% adaptive finite element methods, SIAM J. Numer. Anal., 49(6)2451-2469, 2016.
%
% Function(s) called: marking_strategy_fa 
%
%   TIFISS function: MR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri  

  if nargin < 3
      % Default marking parameter
      threshold = 0.5;
      if nargin < 2
          error('Insufficient input arguments!');
      end
  end

% Strategies
  if strategy == 1 
      % [FPZ16]
      % M_primal and M_dual contain the indices of the marked elements for
      % the primal and the dual problem, respectively; the indices are sorted 
      % in decreasing order w.r.t the value of corresponding indicators
      M_primal = marking_strategy_fa(errsprimal,2,threshold);
      M_dual   = marking_strategy_fa(errsdual,  2,threshold);
      n_marked = min( size(M_primal,1), size(M_dual,1) );
      M = sort( unique([M_primal(1:n_marked); M_dual(1:n_marked)]) );
    
  elseif strategy == 2 
      %[MS09]
      M_primal = marking_strategy_fa(errsprimal,2,threshold);
      M_dual   = marking_strategy_fa(errsdual,  2,threshold);
      if size(M_primal,1) > size(M_dual,1)
          M = sort(M_dual);
      else
          M = sort(M_primal);
      end
      
  elseif strategy == 3 
      % [HP16]
      M_primal = marking_strategy_fa(errsprimal,2,threshold);
      M_dual   = marking_strategy_fa(errsdual,  2,threshold);
      M        = sort( unique([M_primal;M_dual]) );

  elseif strategy == 4 
      % [BET11]
      eta2_primal = sum(errsprimal.^2);
      eta2_dual   = sum(errsdual.^2);
      rho         = sqrt(eta2_dual*errsprimal.^2 + eta2_primal*errsdual.^2);
      M           = marking_strategy_fa(rho,2,threshold);
      M           = sort(M);

  elseif strategy == 5 
      % Only primal
      M = marking_strategy_fa(errsprimal,2,threshold);
      M = sort(M);

  elseif strategy == 6 
      % Only dual
      M = marking_strategy_fa(errsdual,2,threshold);
      M = sort(M);
  end 

end % end function