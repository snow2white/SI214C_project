%GOAFEM_CONVPLOT plots the computed error estimates versus the overall number of dofs
%   TIFISS scriptfile: LR; 22 June 2018
% Copyright (c) 2018 A. Bespalov, D. Praetorius, L. Rocchi, M. Ruggeri

% Set sizes for legend, axis, markers, and linewidth
  fontSizeLegend = 14;
  fontSizeAxis   = 16;
  markerSizePlot = 12;
  lindeWidthPlot = 1.5;
  
% Set of colors
  darkGreen  = [0.0 100 0.0] ./ 255;
  lightGrenn = [0.0 153 0.0] ./ 255;
  brownPeru  = [205 133  63] ./ 255;
  chocolate  = [210 105  30] ./ 255;  
  
% Set limit for degrees of freedom
  xLimInf = (1)*10^(1.0);
  xLimSup = (1)*10^(5.0);
  if sn == 1
      yLimInf = 10^(-6.0);
      yLimSup = 10^(-1);     
  elseif ismember(sn,[2,3])
      yLimInf = 10^(-5.0);
      yLimSup = 10^(-0.5);
  else%sn==4
      yLimInf = 10^(-4.0);
      yLimSup = 10^(-0.5);
  end  

% ----------------------------------------------------------------------------- 
% Plot
% -----------------------------------------------------------------------------
  figure;
  
% Global error estimates (primal) - circles, blue
  loglog(dof,errors_primal,'bo-','Linewidth',lindeWidthPlot,'MarkerSize',markerSizePlot); 
  hold on;

% Global error estimates (dual) - squares, red
  loglog(dof,errors_dual,'rs-','Linewidth',lindeWidthPlot,'MarkerSize',markerSizePlot); 

% Global product estimates - triangles, green
  loglog(dof,errors,'^-','Color',lightGrenn,'Linewidth',lindeWidthPlot,'MarkerSize',markerSizePlot);
  
  grid on;
  axis([xLimInf xLimSup yLimInf yLimSup]);
  xlabel('degrees of freedom $N$','FontSize',fontSizeAxis,'interpreter','latex');
  ylabel('error estimates','FontSize',fontSizeAxis,'interpreter','latex');  

% Reference rates
  if sn == 1
      loglog([xLimInf xLimSup],0.5*[xLimInf xLimSup].^(-1/2),'k-');
      loglog([xLimInf xLimSup],0.03*[xLimInf xLimSup].^(-1.0),'-.','Color',darkGreen);
  elseif ismember(sn,[2,3])
      loglog([xLimInf xLimSup],2.5*[xLimInf xLimSup].^(-1/2),'k-');
      loglog([xLimInf xLimSup],0.1*[xLimInf xLimSup].^(-1.0),'-.','Color',darkGreen);
  else
      loglog([xLimInf xLimSup],2.5*[xLimInf xLimSup].^(-1/2),'k-');
      loglog([xLimInf xLimSup],0.3*[xLimInf xLimSup].^(-1.0),'-.','Color',darkGreen);
  end
  hl = legend('$\mu_\ell$ (primal)','$\zeta_\ell$ (dual)','$\mu_\ell\zeta_\ell$','$N^{ -1/2}$','$N^{ -1}$'); 
  set(hl,'FontSize',fontSizeLegend,'interpreter','latex');
   
% ----------------------------------------------------------------------------- 
% Setup axis, position, and size of the plot windows
% -----------------------------------------------------------------------------
  set(gca,'XTick',[10^1 10^2 10^3 10^4 10^5 10^6],...
          'XTickMode','manual',...
          'XMinorTick','on', 'YMinorTick','on',...
          'XMinorGrid','off','YMinorGrid','off',...
          'GridLineStyle','--','FontSize',fontSizeAxis);
  set(gcf,'units','normalized','Position',[0.25 0.05 0.55 0.8]);
  hold off;

% end scriptfile