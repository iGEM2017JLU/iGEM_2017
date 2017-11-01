function [co,ls]=SAF_NFODE(varargin)
%%--------------------------------------------------------------------------------------------------------------------------------------------------
%
% SAF_NFODE - Simulated annealing fitting for coefficients in nonanalytic first-order ODE models.
%	This function fits coefficients in first-order ordinary differential equations with given datas using simulated annealing method.
%
%--------------------------------------------------------------------------------------------------------------------------------------------------
%
%	[ coef ] = SAF_NFODE ( eqns , var_list , var_ini , coef_list , coef_ini , der_list , region , NAME , RESULTS , ... , TAG , SET , ... )
%	[ coef , l_squ ] = SAF_NFODE ( eqns , var_list , var_ini , coef_list , coef_ini , der_list , region , NAME , RESULTS , ... , TAG , SET , ... )
%
%--------------------------------------------------------------------------------------------------------------------------------------------------
%
%	Inputs:
%
%		eqns: 		first-order ordinary differential equations.
%					E.g. { 'Dy1 == k1 * x + y1 * ( y2 + k2 )' , 'v2 == k2 * y1 + y2 - k1' } (in which Dy1 represents dy1/dx, v2 represents dy2/dx)
%
%		var_list: 	list of variables, has to be symbolic, helps to tell variables from coefficients.
%					E.g. { 'y1' , 'y2' , 'x' }
%                   (notice that 'RK4_STEP''RK4_TOL''SA_BGT''SA_EDT''SA_STEP''SA_TRYTIME''SA_CHANGE' are illegal variable names here)
%					the independent variable has to be the last term on this list
%
%		var_ini:		initial values for all variables, helps to set the initial point for numerical solving.
%					E.g. [ y1_val , y2_val , x_val ] (notice that values must be presented in the same order as var_list)
%
%		coef_list:	name list of coefficients, helps to tell coefficients from variables.
%					E.g. { 'k1' , 'k2' }
%
%		coef_ini:	initial guess for the value of coefficients, not essential, but will help to reach the best fitting value faster.
%					E.g. [ k1_val , k2_val ]
%
%		der_list:		list of derivatives, must be in the same order as var_list.
%					E.g. { 'Dy1' , 'v2' }
%
%		region:		region for the independent variable, helps the ode solver to identify its' solving interval.
%					E.g. [ x_min , x_max ] 
%
%		NAME:		name list of measured variables. At least two terms should be measured in the experiment.
%		RESULTS:	experimental results.
%					E.g. 'x' , [ x_1 , x_2 , x_3 ] , 'y1' , [ y1_1 , y1_2 , y1_3 ] (notice that RESULT must be given in the same order)
%
%		TAG:		additional options.
%		SET:		settings for options.
%					E.g. 'RK4_STEP' , 0.005		-- sets the Runge-Kutta original step to 0.005.
%					E.g. 'RK4_TOL' , 1E-14		-- sets the Runge-Kutta truncation error tolerance to 1E-14, step will be reset due to the tolerance
%					E.g. 'SA_BGT' , 700		-- sets the simulated annealing beginning temperature to 700.
%					E.g. 'SA_EDT' , 5		-- sets the simulated annealing ending temperature to 5.
%					E.g. 'SA_STEP' , 0.5		-- sets the simulated annealing temperature changing to 0.5 per step.
%					E.g. 'SA_TRYTIME' , 15		-- sets the simulated annealing try times to 15 before each step of cooling down.
%					E.g. 'SA_CHANGE', 0.01		-- sets the maximum range of coefficient changings to 0.01
%
%--------------------------------------------------------------------------------------------------------------------------------------------------
%
%	Outputs:
%
%		coef:		fitting results of the coefficients.
%		l_squ:		least squares for the distances from the experimental points to the fitting curves.
%
%--------------------------------------------------------------------------------------------------------------------------------------------------
%
%	For more information, visit our GitHub: <a href="matlab: web('https://github.com/iGEM2017JLU')">iGEM2017JLU</a>.
%	Written by Shuai Wang
%	Jilin University iGEM 2017
%
%--------------------------------------------------------------------------------------------------------------------------------------------------


%% Robustness (input check)
if nargin<6
	error('Input error: not enough arguments!');
else
	odes=varargin{1};varName=varargin{2};varInitial=varargin{3};coefName=varargin{4};
    if ~(iscell(odes)&&iscell(varName)&&isnumeric(varInitial)&&iscell(coefName))
        error('Input error: illegal data type!');
    end
	if length(varName)~=length(varInitial)
		error('Input error: unable to set initial point!');
	end
end
if isnumeric(varargin{5})
    coefInitial=varargin{5};
	if nargin<7
		error('Input error: not enough arguments!');
	else
		if ~iscell(varargin{6})
			error('Input error: illegal data type!');
		else
			derName=varargin{6};
		end
	end
	input_subscript=7;
else
    warning('Input warning: unable to get coeficients initial values.');
	coefInitial=rand(size(coefName));
	if ~iscell(varargin{5})
		error('Input error: illegal data type!')
	else
		derName=varargin{5};
	end
	input_subscript=6;
end
if isnumeric(varargin{input_subscript})&&(prod(size(varargin{input_subscript})==[1,2]))
	reg=varargin{input_subscript};
else
	error('Input error: illegal data type!');
end
rk4step=0.05;rk4tol=1E-6;sa_bgt=400;sa_edt=25;sa_step=5;sa_tt=15;sa_cr=0.002;
if mod(nargin-input_subscript,2)
	error('Input error: inputs don''t match!');
else
	input_subscript=input_subscript+1;
end
datacount=0;
for i=input_subscript:2:nargin
	switch varargin{i}
	case 'RK4_STEP'
		rk4step=varargin{i+1};
	case 'RK4_TOL'
		rk4tol=varargin{i+1};
	case 'SA_BGT'
		sa_bgt=varargin{i+1};
	case 'SA_EDT'
		sa_edt=varargin{i+1};
	case 'SA_STEP'
		sa_step=varargin{i+1};
	case 'SA_TRYTIME'
		sa_tt=varargin{i+1};
	case 'SA_CHANGE'
		sa_cr=varargin{i+1};
    otherwise
		for j=1:length(varName)
			if strcmp(varargin{i},varName{j})
				datacount=datacount+1;
				datasubscript(datacount)=j;
				try
					data(datacount,:)=varargin{i+1};
				catch
					error('Input error: unable to match up experimental results!');
				end
			end
		end
	end
end


%% Optimization
coef_o=coefInitial;
for T=sa_bgt:-sa_step:sa_edt
    for j=1:sa_tt
		coef_c=coef_o-sa_cr/2+rand(size(coef_o))*sa_cr;
		%% solve the equations with rk4 (non-stiff equations only), and calculate the distance from experimental datas to the model curve.
		var_value_o=varInitial;var_value_c=varInitial;x_NotAIndVarName=reg(1,1);len=length(var_value_o);steps=1;
		clear('solution_list_o'); solution_list_o(steps,:)=var_value_o;
		clear('solution_list_c'); solution_list_c(steps,:)=var_value_c;
		while x_NotAIndVarName<reg(1,2)
			%% run the first step of rk4
			[derO_NotADerName1,derC_NotADerName1]=der_cal(odes,coefName,varName,derName,coef_o,coef_c,var_value_o,var_value_c);
			var_value_1_o_m(1:len-1)=var_value_o(1:len-1)+rk4step*derO_NotADerName1/2; % m=middle stepsize
			var_value_1_o_d(1:len-1)=var_value_o(1:len-1)+rk4step*derO_NotADerName1; % d=double stepsize
			var_value_1_o_h(1:len-1)=var_value_o(1:len-1)+rk4step/4*derO_NotADerName1; % h=half stepsize
            var_value_1_o_m(len)=var_value_o(len)+rk4step/2;
            var_value_1_o_d(len)=var_value_o(len)+rk4step;
            var_value_1_o_h(len)=var_value_o(len)+rk4step/4;
            var_value_1_c_m(1:len-1)=var_value_c(1:len-1)+rk4step*derC_NotADerName1/2;
			var_value_1_c_d(1:len-1)=var_value_c(1:len-1)+rk4step*derC_NotADerName1;
			var_value_1_c_h(1:len-1)=var_value_c(1:len-1)+rk4step/4*derC_NotADerName1;
            var_value_1_c_m(len)=var_value_c(len)+rk4step/2;
            var_value_1_c_d(len)=var_value_c(len)+rk4step;
            var_value_1_c_h(len)=var_value_c(len)+rk4step/4;
            %% run the second step of rk4
            [derO_NotADerName2m,derC_NotADerName2m]=der_cal(odes,coefName,varName,derName,coef_o,coef_c,var_value_1_o_m,var_value_1_c_m);
			var_value_2_o_m(1:len-1)=var_value_o(1:len-1)+rk4step*derO_NotADerName2m/2;
			var_value_2_c_m(1:len-1)=var_value_c(1:len-1)+rk4step*derC_NotADerName2m/2;
			var_value_2_o_m(len)=var_value_o(len)+rk4step/2; var_value_2_c_m(len)=var_value_c(len)+rk4step/2;
			[derO_NotADerName2d,derC_NotADerName2d]=der_cal(odes,coefName,varName,derName,coef_o,coef_c,var_value_1_o_d,var_value_1_c_d);
			var_value_2_o_d(1:len-1)=var_value_o(1:len-1)+rk4step*derO_NotADerName2d;
			var_value_2_c_d(1:len-1)=var_value_c(1:len-1)+rk4step*derC_NotADerName2d;
			var_value_2_o_d(len)=var_value_o(len)+rk4step; var_value_2_c_d(len)=var_value_c(len)+rk4step;
			[derO_NotADerName2h,derC_NotADerName2h]=der_cal(odes,coefName,varName,derName,coef_o,coef_c,var_value_1_o_h,var_value_1_c_h);
			var_value_2_o_h(1:len-1)=var_value_o(1:len-1)+rk4step*derO_NotADerName2h/4;
			var_value_2_c_h(1:len-1)=var_value_c(1:len-1)+rk4step*derC_NotADerName2h/4;
			var_value_2_o_h(len)=var_value_o(len)+rk4step/4; var_value_2_c_h(len)=var_value_c(len)+rk4step/4;
			%% run the third step of rk4
			[derO_NotADerName3m,derC_NotADerName3m]=der_cal(odes,coefName,varName,derName,coef_o,coef_c,var_value_2_o_m,var_value_2_c_m);
			var_value_3_o_m(1:len-1)=var_value_o(1:len-1)+rk4step*derO_NotADerName3m/2;
			var_value_3_c_m(1:len-1)=var_value_c(1:len-1)+rk4step*derC_NotADerName3m/2;
			var_value_3_o_m(len)=var_value_o(len)+rk4step/2; var_value_3_c_m(len)=var_value_c(len)+rk4step/2;
			[derO_NotADerName3d,derC_NotADerName3d]=der_cal(odes,coefName,varName,derName,coef_o,coef_c,var_value_2_o_d,var_value_2_c_d);
			var_value_3_o_d(1:len-1)=var_value_o(1:len-1)+rk4step*derO_NotADerName3d;
			var_value_3_c_d(1:len-1)=var_value_c(1:len-1)+rk4step*derC_NotADerName3d;
			var_value_3_o_d(len)=var_value_o(len)+rk4step; var_value_3_c_d(len)=var_value_c(len)+rk4step;
			[derO_NotADerName3h,derC_NotADerName3h]=der_cal(odes,coefName,varName,derName,coef_o,coef_c,var_value_2_o_h,var_value_2_c_h);
			var_value_3_o_h(1:len-1)=var_value_o(1:len-1)+rk4step*derO_NotADerName3h/4;
			var_value_3_c_h(1:len-1)=var_value_c(1:len-1)+rk4step*derC_NotADerName3h/4;
			var_value_3_o_h(len)=var_value_o(len)+rk4step/4; var_value_3_c_h(len)=var_value_c(len)+rk4step/4;
			%% run the last step of rk4
			[derO_NotADerName4m,derC_NotADerName4m]=der_cal(odes,coefName,varName,derName,coef_o,coef_c,var_value_3_o_m,var_value_3_c_m);
			[derO_NotADerName4d,derC_NotADerName4d]=der_cal(odes,coefName,varName,derName,coef_o,coef_c,var_value_3_o_d,var_value_3_c_d);
			[derO_NotADerName4h,derC_NotADerName4h]=der_cal(odes,coefName,varName,derName,coef_o,coef_c,var_value_3_o_h,var_value_3_c_h);
			derO_NotADerNameWm=(derO_NotADerName1+2*derO_NotADerName2m+2*derO_NotADerName3m+derO_NotADerName4m);
			derO_NotADerNameWd=(derO_NotADerName1+2*derO_NotADerName2d+2*derO_NotADerName3d+derO_NotADerName4d);
			derO_NotADerNameWh=(derO_NotADerName1+2*derO_NotADerName2h+2*derO_NotADerName3h+derO_NotADerName4h);
			derC_NotADerNameWm=(derC_NotADerName1+2*derC_NotADerName2m+2*derC_NotADerName3m+derC_NotADerName4m);
			derC_NotADerNameWd=(derC_NotADerName1+2*derC_NotADerName2d+2*derC_NotADerName3d+derC_NotADerName4d);
			derC_NotADerNameWh=(derC_NotADerName1+2*derC_NotADerName2h+2*derC_NotADerName3h+derC_NotADerName4h);
			%% check stepsize
			flag_1=(max(abs(derO_NotADerNameWm-derO_NotADerNameWd))>rk4tol/15)&&(max(abs(derC_NotADerNameWm-derC_NotADerNameWd))>rk4tol);
			flag_2=(max(abs(derO_NotADerNameWm-derO_NotADerNameWh))<rk4tol/15)||(max(abs(derC_NotADerNameWm-derC_NotADerNameWh))<rk4tol);
			switch flag_1-flag_2
			case 1
				rk4step=rk4step/2;
			case -1
				rk4step=rk4step*2;
			otherwise
				steps=steps+1;
				var_value_o=var_value_o+[derO_NotADerNameWm,1]*rk4step; var_value_c=var_value_c+[derC_NotADerNameWm,1]*rk4step;
				solution_list_o(steps,:)=var_value_o; solution_list_c(steps,:)=var_value_c; x_NotAIndVarName=x_NotAIndVarName+rk4step;
			end
		end
		distance_o=0;distance_c=0;[~,b]=size(data);[a,~]=size(solution_list_o);
		for i=1:b
			distance_o=distance_o+min(sum((data(1:datacount,i)*ones(1,a)-solution_list_o').^2));
			distance_c=distance_c+min(sum((data(1:datacount,i)*ones(1,a)-solution_list_c').^2));
		end
		if exp((distance_o-distance_c)/T)>rand(1)
			coef_o=coef_c; distance_o=distance_c;
		end
    end
end
co=coef_o; ls=distance_o;
end



function [derO_NotADerName,derC_NotADerName]=der_cal(odes,coefName,varName,derName,coef_o,coef_c,var_value_o,var_value_c)
for i=1:length(varName)
	eval(['syms ',varName{i},';']);
end
for i=1:length(coefName)
	eval(['syms ',coefName{i},';']);
end
for i=1:length(derName)
	eval(['syms ',derName{i},';']);
end
cmdstr_o='der_o_NotADerName=solve('; cmdstr_c='der_c_NotADerName=solve(';
for i=1:length(odes)
	cmdstr_o=[cmdstr_o,odes{i},',']; cmdstr_c=[cmdstr_c,odes{i},','];
end
for i=1:length(coefName)
	cmdstr_o=[cmdstr_o,coefName{i},'==',num2str(coef_o(i)),','];
	cmdstr_c=[cmdstr_c,coefName{i},'==',num2str(coef_c(i)),','];
end
for i=1:length(varName)
	cmdstr_o=[cmdstr_o,varName{i},'==',num2str(double(var_value_o(i))),','];
	cmdstr_c=[cmdstr_c,varName{i},'==',num2str(double(var_value_c(i))),','];
end
for i=1:length(derName)
	cmdstr_o=[cmdstr_o,derName{i},',']; cmdstr_c=[cmdstr_c,derName{i},','];
end
for i=1:length(coefName)
	cmdstr_o=[cmdstr_o,coefName{i},',']; cmdstr_c=[cmdstr_c,coefName{i},','];
end
for i=1:length(varName)-1
	cmdstr_o=[cmdstr_o,varName{i},',']; cmdstr_c=[cmdstr_c,varName{i},','];
end
cmdstr_c=[cmdstr_c,varName{length(varName)},');'];
cmdstr_o=[cmdstr_o,varName{length(varName)},');'];
try
	eval(cmdstr_o);
	for i=1:length(derName)
		eval(['derO_NotADerName(',num2str(i),')=der_o_NotADerName.',derName{i},';']);
	end
catch
	error(['Calculation Ends: unable to solve derivatives']);
end
try
	eval(cmdstr_c);
	for i=1:length(derName)
		eval(['derC_NotADerName(',num2str(i),')=der_c_NotADerName.',derName{i},';']);
	end
end
end
