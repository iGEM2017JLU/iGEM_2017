function [ population, concentration, time ] = pollutionControl( doses, basics )
%%--------------------------------------------------------------------------------------------------------------------------------------------------
%
% pollutionControl - Giving predictions of pollutant concentration and engineered bacteria population changings from known kinetic laws and doses.
%	This function combines degradation kinetic equations with population growth model and give prediction with the Euler method.
%
%--------------------------------------------------------------------------------------------------------------------------------------------------
%
%	[ population, concentration, time ] = pollutionControl( doses, basics)
%
%--------------------------------------------------------------------------------------------------------------------------------------------------
%
%	Inputs:
%
%		doses:			pollutant and bacteria adding list
%						E.g. see file: settings.m
%
%		basics:			basic settings for numeric solving, including kinetic laws, initial settings, etc.
%						E.g. see file: settings.m
%
%--------------------------------------------------------------------------------------------------------------------------------------------------
%
%	Outputs:
%
%		population:		population predictions
%
%		concentration:	concentration predictions
%
%		time:			time list
%
%--------------------------------------------------------------------------------------------------------------------------------------------------
%
%	For more information, visit our GitHub: <a href="matlab: web('https://github.com/iGEM2017JLU')">iGEM2017JLU</a>.
%	Written by Shuai Wang & Shan Wang
%	Jilin University iGEM 2017
%
%--------------------------------------------------------------------------------------------------------------------------------------------------

a=length(doses);freeGrowth=basics.growthFunc0; inhibitingGrowth=basics.growthFunc1;
w1=basics.weightingFunc0; w2=basics.weightingFunc1;
targetDegradation=basics.degradationFunc;
dt=min([doses(2:a).time_point]-[doses(1:a-1).time_point])/20; time=min([doses(:).time_point]):dt:max([doses(:).time_point])+dt;
concentration=zeros(size(time)+100); concentration(1)=basics.initial_concentration;
population=zeros(1,(length(time)+100)); population(1)=basics.initial_population;
for i=1:length(time)
    concentration(i+1)=concentration(i)+targetDegradation(population(i),concentration(i),dt);
    population(i+1)=population(i)+w1(concentration(i))*freeGrowth(population(i),dt)+w2(concentration(i))*inhibitingGrowth(population(i),dt);
    time_t=time(i);
    check_dose=find(((([doses(:).time_point]-time_t)<=dt).*(([doses(:).time_point]-time_t)>0))==1);
    if ~isempty(check_dose)
        concentration(i+1)=concentration(i+1)+sum([doses(check_dose).concentration_dose]);
        population(i+1)=population(i+1)+sum([doses(check_dose).population_dose]);
    end
end
hold on
h_pop=figure('Name','Population - Time'); plot(time,population(1:length(time)),'b-'); xlabel('Time'); ylabel('E. coli Population');
hold off
hold on
h_con=figure('Name','Concentration - Time'); plot(time,concentration(1:length(time)),'b-'); xlabel('Time'); ylabel('Pollutants Concentration');
hold off
end