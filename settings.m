doses(1).time_point=20; doses(1).population_dose=0.0; doses(1).concentration_dose=0.2;
doses(2).time_point=40; doses(2).population_dose=0.0; doses(2).concentration_dose=0.2;
doses(3).time_point=260; doses(3).population_dose=0.0; doses(3).concentration_dose=0.2;
doses(4).time_point=370; doses(4).population_dose=0.0; doses(4).concentration_dose=0.2;
doses(5).time_point=420; doses(5).population_dose=0.0; doses(5).concentration_dose=0.2;
doses(6).time_point=550; doses(6).population_dose=0.0; doses(6).concentration_dose=0.2;
doses(7).time_point=900; doses(7).population_dose=0.0; doses(7).concentration_dose=0.2;

basics.growthFunc0=@freeGrowthSample; basics.weightingFunc0=@freeGrowthWeightSample;
basics.growthFunc1=@inhibitingGrowthSample; basics.weightingFunc1=@inhibitingGrowthWeightSample;
basics.degradationFunc=@chlorophenolDegradationSample; basics.initial_concentration=0.2; basics.initial_population=0.8;
pollutionControl( doses, basics );