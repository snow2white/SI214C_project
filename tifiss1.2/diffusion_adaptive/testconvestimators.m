%TESTCONVESTIMATORS plots the convergence of energy error estimators (example 3 and 5)
%    TIFISS scriptfile: LR; 09 July 2018
% Copyright (c) 2017 D.J. Silvester, A. Bespalov

% Legend Estimators: 
%     errbub: hierarchical eY with quadratic functions
%    errhat1: hierarchical eY with p/w linear functions (elementwise error residuals)
%    errhat2: hierarchical eY with p/w linear functions (fully linear system)
%    errhat3: 2-level estimator with p/w linear functions

  close all;

% Colors
  darkBlue   = [0   141 226] ./ 255;
  darkGreen  = [0.0 100 0.0] ./ 255;
  lightGrenn = [0.0 153 0.0] ./ 255;
  brownPeru  = [205 133  63] ./ 255;
  chocolate  = [210 105  30] ./ 255;
  
%% Test problem 3: square domain, regular solution
%  coeff: a(x,y) = 1 + C*fx*fy
%    rhs: f(x,y) = (2-x^2-y^2)/8 
% using:
%   tolerance = 1e-03
%    strategy = 2 
%   threshold = 0.5
%    grid par = 3

% Axis
  xLimInf = (1)*10^(1);  
  xLimSup = (6)*10^(4);
  yLimInf = (3)*10^(-4);
  yLimSup = (3)*10^(-2);
%
% Hierarchical P2 bubbles     
  intdofbub  = [49    62    87   115   155   212   290   391   515   696   941  1247  1658  2189  2900  3882  5081  6688  8647 11330]; 
  errbub     = [0.0097135844 0.0069064726 0.0057730114 0.0051378408 0.0043521274 0.0039071478 0.0032323908 0.0029313474 0.0024152565 0.0020710884 0.0019125408 0.0016288844 0.0014009302 0.0012032556 0.0010398477 0.0009526496 0.0008148404 0.0007135747 0.0006136657 0.0005335536];
%
% Hierarchical P1 bubbles - elementwise (only marking edges)
  intdofhat1 = [49   60   80  103  139  186  255  347  449  596  803 1098 1477 1937 2570 3428 4664 6175];
  errhat1    = [0.0072705592 0.0055261045 0.0045025305 0.0038618751 0.0033320275 0.0029781014 0.0025185477 0.0022498791 0.0019111880 0.0016268051 0.0014609670 0.0012555688 0.0010879943 0.0009279781 0.0008022758 0.0007188365 0.0006167633 0.0005346929];
%
% Hierarchical P1 bubbles - fully assembled linear system (marking edges)
  intdofhat2 = [49   59    77   101   138   187   264   357   467   623   857  1185  1566  2062  2753  3793  5168];
  errhat2    = [0.0062515836 0.0050375844 0.0042600515 0.0034966174 0.0029953377 0.0026300307 0.0022983804 0.0019702222 0.0016789471 0.0014480681 0.0012662813 0.0011080116 0.0009517848 0.0008119406 0.0007020317 0.0006161818 0.0005341619];
% 
% 2-level estimator (marking edges)
  intdofhat3 = [49    58    82   120   175   271   414   619   918  1402  2054  3042  4580];
  errhat3    = [0.0057323133 0.0046144968 0.0037044407 0.0029593499 0.0024038642 0.0019901418 0.0016628127 0.0013533983 0.0011038953 0.0009110530 0.0007535248 0.0006139137 0.0005037473];

% Plot 
  figure;
  loglog(intdofbub,errbub,'r-^','MarkerSize',11); hold on; grid on;
  loglog(intdofhat1,errhat1,'bo-','MarkerSize',11); 
  loglog(intdofhat2,errhat2,'s-','MarkerSize',11,'Color',darkBlue);
  loglog(intdofhat3,errhat3,'d-','MarkerSize',11,'Color',darkGreen);
  plot([xLimInf xLimSup],0.1*[xLimInf xLimSup].^(-1/2),'-','Color',[210 105 30]./255); 
  xlabel('degrees of freedom','Fontsize',16);
  ylabel('estimated energy error','Fontsize',16);
  hl = legend('Hierarchical eY (elementwise) P2 bubbles',        ...
              'Hierarchical eY (elementwise) P1 linear bubbles', ...
              'Hierarchical eY (full system) P1 linear bubbles', ...
              '2-Level estimator P1 linear bubbles',             ...
              'N^{-1/2}');   
  title('adaptive test problem 3','Fontsize',18);
  axis([xLimInf xLimSup yLimInf yLimSup]);
% Setting
  set(hl,'FontSize',17);
  set(gca,'XTick',[10^1 10^2 10^3 10^4 10^5],...
          'XTickMode','manual','XMinorTick','on', 'YMinorTick','on',...
          'XMinorGrid','off','YMinorGrid','off','GridLineStyle','--','FontSize',16);
  set(gcf,'units','normalized','Position',[0.25 0.05 0.55 0.8]);
  hold off;   
 
%% Test problem 5: L-shaped domain, singular solution; 
%% see ADAPT_DIFF_TESTPROBLEM, example 5: --> problem=2; <--
%  singular solution: u(r,theta) = r^(2/3)sin((pi + 2theta)/3) with
%  coeff: a(x,y) = 1
%    rhs: f(x,y) = 0
% -> Lap( u(r,theta) ) = 0
% using:
%   tolerance = 3e-02
%    strategy = 2 
%   threshold = 0.5
%    grid par = 3 

% Axis
  xLimInf = (1)*10^(1);  
  xLimSup = (5)*10^(4);
  yLimInf = (3)*10^(-4);
  yLimSup = (1)*10^(0);
 
% Hierarchical P2 bubbles (only marking elements)
  intdofbub  = [33    35    38    49    62    81   106   148   208   281   382   526   734  1013  1351  1827  2434  3296  4375  5765  7601  9906 12957];      
  errbub     = [0.2379601726 0.2098318888 0.1817686512 0.1579158063 0.1360038942 0.1210786032 0.1018579213 0.0874061946 0.0748844325 0.0656881472 0.0553565339 0.0478396672 0.0409017736 0.0348717397 0.0299325506 0.0261877132 0.0224838439 0.0197105849 0.0169009904 0.0147880843 0.0128554421 0.0113839955 0.0098785874];
  energybub  = [1.3699499153 1.3670088505 1.3641436246 1.3619863690 1.3606112035 1.3595665017 1.3587307022 1.3580738103 1.3575934000 1.3572619159 1.3570097624 1.3564976171 1.3558664055 1.3556421589 1.3555769510 1.3555007104 1.3553925968 1.3552441787 1.3552104842 1.3551965663 1.3551767722 1.3551514046 1.3551164927];
%
% Hierarchical P1 bubbles - elementwise (only marking edges)
  intdofhat1 = [33    36    40    49    64    85   123   161   219   297   398   545   732  1010  1364  1848  2486  3331  4460  5916  7791];
  errhat1    = [0.1669427393 0.1472998523 0.1293089154 0.1137109511 0.0998255949 0.0857953257 0.0722568340 0.0653338933 0.0545596459 0.0474867419 0.0407773116 0.0356887288 0.0299877970 0.0263478732 0.0225254813 0.0194019833 0.0167201394 0.0144271933 0.0125862722 0.0108724826 0.0095495264];  
  energyhat1 = [1.3699499153 1.3663383364 1.3639765631 1.3619863690 1.3605524965 1.3594032995 1.3585271585 1.3580050306 1.3575372273 1.3572181746 1.3569869789 1.3562145256 1.3558521238 1.3556456116 1.3555761984 1.3555261808 1.3553110423 1.3552283021 1.3552094927 1.3551950825 1.3551684175];
%
% Hierarchical P1 bubbles - fully assembled linear system (marking edges)
  intdofhat2 = [33    35    38    43    52    73    99   139   193   262   355   483   658   893  1215  1619  2163  2905  3884  5187];
  errhat2    = [0.1431758253 0.1286609064 0.1128924167 0.1018848310 0.0898893160 0.0784365506 0.0679673061 0.0578586697 0.0498478385 0.0427944017 0.0368882891 0.0320476479 0.0274400424 0.0236868994 0.0202912795 0.0176630159 0.0153317237 0.0132223067 0.0114769665 0.0099256191];
  energyhat2 = [1.3699499153 1.3670088505 1.3641436246 1.3624568275 1.3609672289 1.3597847155 1.3588790929 1.3581694117 1.3576719148 1.3571353452 1.3563267132 1.3559312529 1.3557771681 1.3556765274 1.3555544426 1.3553719981 1.3552736032 1.3552406054 1.3552187506 1.3551830457];
%
% Hierarchical P1 bubbles - fully assembled linear system (marking elements)
  intdofhat2_elem = [33    35    38    43    52    73    99   139   193   262   355   483   658   893  1215  1619  2163  2905  3884  5187];
  errhat2_elem    = [0.1431758253 0.1286609064 0.1128924167 0.1018848310 0.0898893160 0.0784365506 0.0679673061 0.0578586697 0.0498478385 0.0427944017 0.0368882891 0.0320476479 0.0274400424 0.0236868994 0.0202912795 0.0176630159 0.0153317237 0.0132223067 0.0114769665 0.0099256191];
  energyhat2_elem = [1.3699499153 1.3670088505 1.3641436246 1.3624568275 1.3609672289 1.3597847155 1.3588790929 1.3581694117 1.3576719148 1.3571353452 1.3563267132 1.3559312529 1.3557771681 1.3556765274 1.3555544426 1.3553719981 1.3552736032 1.3552406054 1.3552187506 1.3551830457];
%
% 2-level estimator (marking edges)
  intdofhat3 = [33    40    52    74   117   189   311   516   858  1418  2365  3922  6484];
  errhat3    = [0.1322551777 0.1120226030 0.0940571647 0.0771212676 0.0627339618 0.0494093881 0.0389866128 0.0306243357 0.0237283658 0.0184860220 0.0143980114 0.0111606134 0.0086799064];
  energyhat3 = [1.3699499153 1.3650789357 1.3620386222 1.3598341389 1.3584784605 1.3576198173 1.3571068921 1.3567882884 1.3559339096 1.3555561350 1.3554499514 1.3552272006 1.3551895211];
%
% 2-level estimator (marking elements)
  intdofhat3_elem = [33    36    42    50    64    84   119   161   223   310   422   576   785  1063  1447  1945  2619  3518  4719  6318];
  errhat3_elem    = [0.1322551777 0.1180910649 0.1060897547 0.0951448969 0.0842528908 0.0741869401 0.0632434977 0.0549857412 0.0469996256 0.0405662359 0.0349555271 0.0300800112 0.0257347095 0.0220121589 0.0188745157 0.0163720852 0.0141149973 0.0121862919 0.0104906917 0.0090794653];
  energyhat3_elem = [1.3699499153 1.3663383364 1.3637817175 1.3619528731 1.3605739583 1.3595059300 1.3585955456 1.3580093853 1.3575437353 1.3572113547 1.3569728278 1.3563284932 1.3558471460 1.3556378007 1.3555698903 1.3555132125 1.3553473517 1.3552408658 1.3552078878 1.3551937971];
  
% EXTRA P2 APPROXIMATIONS
  intdofP2   = [161   165   177   181   185   199   203   212   225   239   259   297   337   399   493   601   755   935  1195  1497  1835  2271  2809  3489  4319  5363  6663];
  errP2      = [0.0743761520 0.0684597584 0.0563806431 0.0490173829 0.0455795617 0.0376907189 0.0333792308 0.0298149856 0.0245936356 0.0218049683 0.0188752135 0.0157272619 0.0130806002 0.0107641532 0.0086790411 0.0071408516 0.0057908042 0.0046896782 0.0037827597 0.0030566074 0.0024777156 0.0020044305 0.0016144029 0.0013073722 0.0010539133 0.0008472771 0.0006817326];
  energyP2   = [1.3577368248 1.3572077745 1.3565241552 1.3562056428 1.3560006556 1.3556919989 1.3555657173 1.3554436172 1.3553119307 1.3552468293 1.3552019216 1.3551618610 1.3551340855 1.3551142265 1.3551008567 1.3550920630 1.3550860320 1.3550821387 1.3550794267 1.3550777713 1.3550768228 1.3550762604 1.3550756549 1.3550751299 1.3550747974 1.3550746631 1.3550745776];

% Plot 
  figure;
  loglog(intdofbub,errbub,'r-^','MarkerSize',11); hold on; grid on;
  loglog(intdofhat1,errhat1,'bo-','MarkerSize',11); 
  loglog(intdofhat2,errhat2,'s-','MarkerSize',11,'Color',darkBlue);
  loglog(intdofhat3,errhat3,'d-','MarkerSize',11,'Color',darkGreen);  
  loglog(intdofP2,errP2,'*-','MarkerSize',11,'Color',brownPeru);
  plot([xLimInf xLimSup],2*[xLimInf xLimSup].^(-1/2),'-','Color',[210 105 30]./255); 
  xlabel('degrees of freedom','Fontsize',16);
  ylabel('estimated energy error','Fontsize',16);
  hl = legend('Hierarchical eY - P2 bubbles (elementwise)', ...
              'Hierarchical eY - P1 bubbles (elementwise)', ...
              'Hierarchical eY - P1 bubbles (full system)', ...
              '2-Level estimator - P1 bubbles',             ...
              'P2 approximations',                          ...
              'N^{-1/2}'); 
  title('adaptive test problem 5','Fontsize',18);
  axis([xLimInf xLimSup yLimInf yLimSup]);
% Setting
  set(hl,'FontSize',15);
  set(gca,'XTick',[10^1 10^2 10^3 10^4 10^5],...
          'XTickMode','manual','XMinorTick','on', 'YMinorTick','on',...
          'XMinorGrid','off','YMinorGrid','off','GridLineStyle','--','FontSize',16);
  set(gcf,'units','normalized','Position',[0.25 0.05 0.55 0.8]);
  hold off;  
   
% Effectivity indices for example 5 - Solution is known:
% Energy norm of the solution (see above) on the L-shaped domain:
  refenergy = 1.355074407666495783786331230658106505870819091796875e+00;
  
  effbub  = errbub  ./ sqrt( abs(refenergy^2 - energybub.^2)  );
  effhat1 = errhat1 ./ sqrt( abs(refenergy^2 - energyhat1.^2) );
  effhat2 = errhat2 ./ sqrt( abs(refenergy^2 - energyhat2.^2) );  % (for marking edges)
  effhat3 = errhat3 ./ sqrt( abs(refenergy^2 - energyhat3.^2) );  % (for marking edges)
  effP2   = errP2   ./ sqrt( abs(refenergy^2 - energyP2.^2)   );
  
% Plot the effectivity indices  
  figure;
  loglog(intdofbub,effbub,  'r-^','Linewidth',1.5,'MarkerSize',10); hold on; grid on;
  loglog(intdofhat1,effhat1,'bo-','Linewidth',1.5,'MarkerSize',10); 
  loglog(intdofhat2,effhat2,'s-','Linewidth',1.5,'MarkerSize',10,'Color',darkBlue);
  loglog(intdofhat3,effhat3,'d-','Linewidth',1.5,'MarkerSize',10,'Color',darkGreen);  
  loglog(intdofP2,effP2,    '*-','Linewidth',1.5,'MarkerSize',10,'Color',brownPeru);
%
  axis([xLimInf xLimSup 0.4 1.35]);
  title('adaptive test problem 5 - effectivity indices','Fontsize',18);
  xlabel('degrees of freedom','Fontsize',16);
  ylabel('effectivity index','FontSize',16);
  hl = legend('Hierarchical eY - P2 bubbles (elementwise)', ...
              'Hierarchical eY - P1 bubbles (elementwise)', ...
              'Hierarchical eY - P1 bubbles (full system)', ...
              '2-Level estimator - P1 bubbles',             ...
              'P2 approximations');
% Setting
  set(hl,'FontSize',14);         
  set(gca,'XTickMode','manual','XMinorTick','on', 'YMinorTick','on','XMinorGrid','off','YMinorGrid','off','GridLineStyle','--','FontSize',14);      
  set(gcf,'units','normalized','Position',[0.25 0.05 0.55 0.8]);
  
% end scriptfile