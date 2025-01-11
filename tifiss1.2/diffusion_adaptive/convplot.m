%CONVPLOT plots the computed energy error estimates versus the overall number of dof
%   TIFISS scriptfile: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, L. Rocchi

% Colors
  darkBlue   = [0   141 226] ./ 255;
  darkGreen  = [0.0 100 0.0] ./ 255;
  lightGrenn = [0.0 153 0.0] ./ 255;
  brownPeru  = [205 133  63] ./ 255;
  chocolate  = [210 105  30] ./ 255;

% Set axis limit according to the model problem
  xLimInf = (5)*10^(0);  
  xLimSup = (2)*10^(5);
  if sn == 1
      % Square domain (0,1)^2
      if pmethod == 1, yLimInf = (8)*10^(-4.0); yLimSup = (2)*10^(-1.0); C = 0.7; 
      else             yLimInf = (2)*10^(-5.0); yLimSup = (2)*10^(-2.0); C = 1.9; 
      end 
  elseif sn == 2 
      % Square domain (0,1)^2
      if pmethod == 1, yLimInf = (1)*10^(-3.0); yLimSup = (1)*10^(-1);   C = 0.3;
      else             yLimInf = (4)*10^(-5.0); yLimSup = (2)*10^(-1.0); C = 1.1;
      end
  elseif sn == 3
      % Square domain (0,1)^2
      if pmethod == 1, yLimInf = (2)*10^(-4.0); yLimSup = (3)*10^(-2);   C = 0.1;
      %P2 for variable coefficients not supported yet
      %else            yLimInf = (4)*10^(-5.0); yLimSup = (3)*10^(-2.0); C = 0.6;
      end    
  elseif sn == 4
      % L-shaped domain
      if pmethod == 1, yLimInf = (4)*10^(-3.0); yLimSup = (3)*10^(-1.0); C = 1.8;
      else             yLimInf = (3)*10^(-4.0); yLimSup = (1)*10^(-1.0); C = 1.2;
      end
  elseif sn == 5
      % L-shaped domain
      if pmethod == 1, yLimInf = (4)*10^(-3.0); yLimSup = (3)*10^(-1.0); C = 1.8;
      else             yLimInf = (2)*10^(-4.0); yLimSup = (3)*10^(-1.0); C = 1.2;
      end   
  else
      % Crack domain
      if pmethod == 1, yLimInf = (6)*10^(-3.0); yLimSup = (3)*10^(-1.0); C = 2.1;
      else             yLimInf = (4)*10^(-4.0); yLimSup = (1)*10^(-1.0); C = 1.2;
      end
  end  
    
% Plot the refinement path with dof
  figure;
  loglog(intdofs + nnzerobdofs,comperror,'b-o','MarkerSize',13); 
  grid on; hold on; 
  title(['P',num2str(pmethod),' | Refinement path'],'Fontsize',16);
  xlabel('degrees of freedom N','Fontsize',16);  %,'interpreter','latex');
  ylabel('energy error estimate','Fontsize',16); %,'interpreter','latex');
  
% Reference optimal rates dof^(-1/2) (for P1) and dof^(-1) (for P2)
  if pmethod == 1
      plot([xLimInf xLimSup],C*[xLimInf xLimSup].^(-1/2),'-','Color',[210 105 30]./255);
      hl = legend('energy error estimate','N^{-1/2}');
  else
      plot([xLimInf xLimSup],C*[xLimInf xLimSup].^(-1),'-','Color',[210 105 30]./255);
      hl = legend('energy error estimate','N^{-1}');
  end
% Axis  
  axis([xLimInf xLimSup yLimInf yLimSup]);
% Setting
  set(hl,'FontSize',17); %,'interpreter','latex');
  set(gca,'XTick',[10^1 10^2 10^3 10^4 10^5 10^6],...
          'XTickMode','manual',...
          'XMinorTick','on', 'YMinorTick','on',...
          'XMinorGrid','off','YMinorGrid','off',...
          'GridLineStyle','--','FontSize',16);
  set(gcf,'units','normalized','Position',[0.25 0.05 0.55 0.8]);
  hold off;

% end scriptfile 