function [fit_result, goodness] = population_fit ( time , population )
%%--------------------------------------------------------------------------------------------------------------------------------------------------
%
%population_fit - fits the E.coli growth model to experimental data.
%
%--------------------------------------------------------------------------------------------------------------------------------------------------
%
%   [fit_result, goodness] = population_fit ( time , population )
%
%--------------------------------------------------------------------------------------------------------------------------------------------------
%
%	Inputs:
%
%		time:		time points for measurements (unit = hour)
%					E.g. [0,1,2,3]
%
%		population:	OD600 or CFU readings, must be in the same length as 'time'
%					E.g. [0.1,0.4,0.7,0.9]
%
%--------------------------------------------------------------------------------------------------------------------------------------------------
%
%	Outputs:
%
%		fit_result:	fitting result given by non-linear least squares fitting method
%
%		goodness:	fitting goodness
%
%--------------------------------------------------------------------------------------------------------------------------------------------------
%
%	For more information, visit our GitHub: <a href="matlab: web('https://github.com/iGEM2017JLU')">iGEM2017JLU</a>.
%	Written by Shuai Wang
%	Jilin University iGEM 2017
%
%--------------------------------------------------------------------------------------------------------------------------------------------------

[xData, yData] = prepareCurveData( time , population );
ft = fittype( 'Nec_*exp(k_*x)*N0/(Nec_-N0)/(1+exp(k_*x)*N0/(Nec_-N0))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0]; opts.Upper = [2 5 2];
opts.StartPoint = [0.5 2 0.02];
[fit_result, goodness] = fit( xData, yData, ft, opts );

% Plot fitting results.
figure( 'Name', 'E.coli Growth Model Fit' );
h=plot( fit_result, xData, yData );
legend( h, 'population vs. time', 'model curve', 'Location', 'SouthEast' );
xlabel('time (hour)');
ylabel('population');
grid on


end