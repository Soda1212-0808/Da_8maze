function plot3D(seq, xspec, varargin,nPlotMax)
%
% plot3D(seq, xspec, ...)
%
% Plot neural trajectories in a three-dimensional space.
%
% INPUTS:
%
% seq        - data structure containing extracted trajectories
% xspec      - field name of trajectories in 'seq' to be plotted 
%              (e.g., 'xorth' or 'xsm')
%
% OPTIONAL ARGUMENTS:
%
% dimsToPlot - selects three dimensions in seq.(xspec) to plot 
%              (default: 1:3)
% nPlotMax   - maximum number of trials to plot (default: 20)
% redTrials  - vector of trialIds whose trajectories are plotted in red
%              (default: [])
%
% @ 2009 Byron Yu -- byronyu@stanford.edu

  dimsToPlot = 1:3;
%   nPlotMax   = 71;
  redTrials  = [1 ];
%    assignopts(who, varargin);

 
  
  if size(seq(1).(xspec), 1) < 3
    fprintf('ERROR: Trajectories have less than 3 dimensions.\n');
    return
  end

  f = figure;
  pos = get(gcf, 'position');
  set(f, 'position', [pos(1) pos(2)-300 1.3*pos(3) 1.3*pos(4)]);
  mm=max(length(seq), nPlotMax);

  
  for n = 1:max(length(seq), nPlotMax)
      datt = seq(n).(xspec)(dimsToPlot,:);
      T   = seq(n).T;
      
      if ismember(seq(n).trialId, redTrials)
          if strcmp(seq(n).event,'n1')
              col =  [0.5 0.5 1] ;
              lw=0.5;
          elseif strcmp(seq(n).event,'n2')
              col =  [0 0 0.5] ;
              lw=2;
          elseif strcmp(seq(n).event,'s1')
              col =  [1 0.9  0.9] ;
              lw=0.5;
          elseif strcmp(seq(n).event,'s2')
              col =  [0.5 0 0] ;
              lw=2;
          end
          
      elseif strcmp(seq(n).event,'n1')
          col =  [0.5 0.5 1] ;
          lw=0.5;
      elseif strcmp(seq(n).event,'n2')
          col =  [0 0 1] ;
          lw=0.5;
      elseif strcmp(seq(n).event,'s1')
          col =  [1 0.5  0.5] ;
          lw=0.5;
      elseif strcmp(seq(n).event,'s2')
          col =  [1 0 0] ;
          lw=0.5;
          
      end
      plot3(datt(1,:), datt(2,:), datt(3,:), '-', 'linewidth', lw, 'color', col);
      
        hold on;
     scatter3(datt(1,21), datt(2,21), datt(3,21),'o','filled','MarkerFaceColor',col)
%       scatter3(datt(1,1), datt(2,1), datt(3,1),'diamond','MarkerFaceColor',col)
  
      middle_dat(n,:)=datt(1:3,1);
    
      
  end
% plot3(middle_dat(:,1),middle_dat(:,2),middle_dat(:,3), '-', 'linewidth',3)
  
  axis equal;
  if isequal(xspec, 'xorth')
    str1 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(1));
    str2 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(2));
    str3 = sprintf('$$\\tilde{\\mathbf x}_{%d,:}$$', dimsToPlot(3));
  else
    str1 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(1));
    str2 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(2));
    str3 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(3));
  end
set(gca,'FontSize',16,'LineWidth',1.5,'FontName','arial')
  
  xlabel(str1, 'interpreter', 'latex', 'fontsize', 24);
  ylabel(str2, 'interpreter', 'latex', 'fontsize', 24);
  zlabel(str3, 'interpreter', 'latex', 'fontsize', 24);
