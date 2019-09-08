function [R, p,R_spearman, p_spearman] = Plot_Fit_ConfidenceInterval(xdata,ydata,degree,alpha)
%H = figure
% xdata = GA;
% ydata = Parenchyma;
% 
% degree = 2;		% Degree of the fit
% alpha = 0.05;	% Significance level


% Compute the fit and return the structure used by 
% POLYCONF.
[p,S] = polyfit(xdata,ydata,degree);

% Compute the real roots and determine the extent of the 
% data.
r = roots(p)'; 							% Roots as a row vector.
real_r = r(imag(r) == 0);               % Real roots.

% Assure that the data are row vectors.
xdata = reshape(xdata,1,length(xdata));
ydata = reshape(ydata,1,length(ydata));

% Extent of the data.
mx = min([real_r,xdata]);
Mx = max([real_r,xdata]);
my = min([ydata,0]);
My = max([ydata,0]);

% Scale factors for plotting.
sx = 0.05*(Mx-mx);
sy = 0.05*(My-my);

% Plot the data, the fit, and the roots.
% % hdata = plot(xdata,ydata,'kd','MarkerSize',5,...
% % 		'LineWidth',2);
hdata = plot(xdata,ydata,'b.','MarkerSize',30,...
    'LineWidth',4,'MarkerFaceColor','b','MarkerEdgeColor','b');
hold on
xfit = mx-sx:0.01:Mx+sx;
yfit = polyval(p,xfit);
hfit = plot(xfit,yfit,'k-','LineWidth',4);
%hroots = plot(real_r,zeros(size(real_r)),...
 %             'bo','MarkerSize',5,...
  %            'LineWidth',2,...
   %           'MarkerFaceColor','b');
grid on
%plot(xfit,zeros(size(xfit)),'k-','LineWidth',2)
axis([min(xdata) max(xdata)+sx my-sy My+sy])

% Add prediction intervals to the plot.
[Y,DELTA] = polyconf(p,xfit,S,'alpha',alpha);
hconf = plot(xfit,Y+DELTA,'k--','LineWidth',2);
plot(xfit,Y-DELTA,'k--','LineWidth',2)
set(gca,'LineWidth',2,'FontSize',24)

% Display the polynomial fit and the real roots.
%approx_p = round(100*p)/100; % Round for display.
%htitle = title(['{\bf Fit:   }',texlabel(polystr(approx_p))]);
%set(htitle,'Color','b')
%approx_real_r = round(100*real_r)/100; % Round for display.
%hxlabel = xlabel(['{\bf Real Roots:     }',num2str(approx_real_r)]);
%set(hxlabel,'Color','b')
[R p]=corrcoef(xdata,ydata);
[R_spearman p_spearman] = corr(xdata',ydata','Type','Spearman');
% Add a legend.
%legend([R,p],'R','p-value');

