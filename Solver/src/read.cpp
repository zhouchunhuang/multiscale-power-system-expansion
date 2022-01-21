#include "model.h"
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <map>

#include <string>
#include <sstream>
#include <random>

using namespace std;

void ReadWrite::readParameters()
{
	auto paramReg = ParamRegistry::instance();
	paramFile = inputDirectory + paramFile;
	ifstream inputFile(paramFile.c_str());
	if (!inputFile)
	{
		std::cout << "Couldn't find the parameter file. Will use the default values!" << std::endl;
	}
	else
	{
		while (!inputFile.eof())
		{
			string	fieldName;
			inputFile >> fieldName;
			if (fieldName == "END")
				break;
			else if (fieldName == "Periods")
			{
			   inputFile >> paramReg->numPeriods;
			   paramReg->numScenarios = int(pow(double(paramReg->numRealizations), paramReg->numPeriods - 1));
			   paramReg->prob_scenario = double(1.0 / double(paramReg->numScenarios));
			   continue;
			}
			else if (fieldName == "DmdScns")
			{
			   inputFile >> paramReg->numDmdScenarios;
			   paramReg->prob_dmd_scenario = double(1.0 / double(paramReg->numDmdScenarios));
			   paramReg->numSPsPerSeason = paramReg->numDmdScenarios;
			   continue;
			}
			else if (fieldName == "TimeLimit")
			{
			   inputFile >> paramReg->time_lmt;
				continue;
			}
			else if (fieldName == "SmpTimeLimit")
			{
			   inputFile >> paramReg->smp_timlmt;
				continue;
			}
			else if (fieldName == "SpTimeLimit")
			{
			   inputFile >> paramReg->sp_tmlmt;
				continue;
			}			
			else if (fieldName == "MpGapTol")
			{
			   inputFile >> paramReg->mp_gap;
				continue;
			}
			else if (fieldName == "SmpGapTol")
			{
			   inputFile >> paramReg->smp_gap;
				continue;
			}
			else if (fieldName == "MaxSmpGap")
			{
			   inputFile >> paramReg->max_smp_gap;
				continue;
			}
			else if (fieldName == "SpGapTol")
			{
			   inputFile >> paramReg->sp_gap;
				continue;
			}
			else if (fieldName == "KeepBendersCuts")
			{
			   inputFile >> paramReg->keep_benders_cuts;
				continue;
			}
			else if (fieldName == "RemoveCutsThreshold")
			{
			   inputFile >> paramReg->remove_cuts_threshold;
			   continue;
			}
			else {
				continue;
			}
		}
	}
}

void ReadWrite::readSystemData(Model &_model)
{
	busFile = inputDirectory + busFile;
	genFile = inputDirectory + genFile;
	lineFile = inputDirectory + lineFile;
	stgFile = inputDirectory + stgFile;
	weeklyLoadFile = inputDirectory + weeklyLoadFile;
	dailyLoadFile = inputDirectory + dailyLoadFile;
	hourlyLoadFile = inputDirectory + hourlyLoadFile;

	readBusData(_model);
	readGenerators(_model);
	readRenewable(_model);
	readLines(_model);
	readStorages(_model);
	readLoadData(_model);
}

void ReadWrite::readBusData(Model &_model)
{
	_model.buses.clear();
	ifstream f_bus(busFile.c_str());												//open the csv file
	if (!f_bus.good())
	{
		cerr << "CANNOT OPEN THE BUS FILE!" << endl;
		return;
	}
	int row = 0, col = 0;
	while (f_bus.good()) {
		row++;
		Bus newBus(_model.numNodes);
		string s;
		if (!getline(f_bus, s)) break;

		col = 0;
		istringstream ss(s);
		while (ss)
		{
			col++;
			string value;
			if (!getline(ss, value, ',')) break;
			if (row > 1)
			{
				switch (col)
				{
				case 1:
					newBus.id = atoi(value.c_str());
					break;
				case 2:
					newBus.set_peak_load(atof(value.c_str()));
					break;
				default:
					break;
				}
			}
		}
		if (row > 1)
		{			
			_model.buses.push_back(newBus);
		}
	}
	_model.numBus = static_cast<int>(_model.buses.size());
	f_bus.close();
}

void ReadWrite::readLoadData(Model &_model)
{
	if (_model.buses.empty()) {
		cerr << "Error: Load data have to be read after buses are initiated!";
		return;
	}
	auto &buses = _model.buses;

	const auto paramReg = ParamRegistry::instance();
	std::string line;
	int row = 0;
	std::vector<double> weekly_peak_dmd_percentage, daily_peak_dmd_percentage;
	std::vector<std::vector<double> > hrly_peak_dmd_percentage;

	std::ifstream wkly_load(weeklyLoadFile.c_str());
	if (!wkly_load.good()) {
		std::cerr << "CANNOT OPEN THE WEEKLY LOAD FILE!" << '\n';
		return;
	}
	line = "";
	row = 0;
	while (std::getline(wkly_load, line)) {
		row++;
		const std::vector<std::string> load_data = split(line);
		if (row >= 5) {
			weekly_peak_dmd_percentage.push_back(std::stof(load_data[1]) * 0.01);
		}
	}
	wkly_load.close();

	std::ifstream daily_load(dailyLoadFile.c_str());
	if (!daily_load.good()) {
		std::cerr << "CANNOT OPEN THE DAILY LOAD FILE!" << '\n';
		return;
	}
	line = "";
	row = 0;
	while (std::getline(daily_load, line)) {
		row++;
		const std::vector<std::string> load_data = split(line);
		if (row >= 5) {
			daily_peak_dmd_percentage.push_back(std::stof(load_data[1]) * 0.01);
		}
	}
	daily_load.close();

	std::ifstream hourly_load(hourlyLoadFile.c_str());
	if (!hourly_load.good()) {
		std::cerr << "CANNOT OPEN THE HOURLY LOAD FILE!" << '\n';
		return;
	}
	line = "";
	row = 0;
	while (std::getline(hourly_load, line)) {
		row++;
		const std::vector<std::string> load_data = split(line);
		if (row >= 7) {
			std::vector<double> percentages;
			for (int i = 1; i < 7; i++) {
				percentages.push_back(std::stof(load_data[i]) * 0.01);
			}
			hrly_peak_dmd_percentage.push_back(percentages);
		}
	}
	hourly_load.close();

	// generate hourly loads for a year per bus
	for (auto & bus : buses) 
	{
		std::vector<double> node_peak_load = std::vector<double>(_model.numNodes);
		bus.season_hourly_loads.clear();
		for (int n = 0; n < _model.numNodes; n++)
		{
			const auto &node = _model.scenarioNodes[n];
			const int period = node.period;
			// future load is uncertain
			const double factor = paramReg->demandChange[n%paramReg->numRealizations];
			if (n == 0) {
				node_peak_load[n] = bus.get_peak_load();
			}else{
			   node_peak_load[n] = node_peak_load[node.parent] * factor;
			}
			bus.peak_loads_in_scenario_nodes[n] = int(node_peak_load[n]);

			vector2int season_load(Season::numSeasons);
			// winter week
			for (int h = 0; h < _model.ucPeriods; h++) {
				//const int day = int(h / HOURSPERDAY);
				const int hour = int(h % HOURSPERDAY);
				// double daily_percentage = daily_peak_dmd_percentage[day];
				//double hrly_percentage = day < 5 ? hrly_peak_dmd_percentage[hour][0] : hrly_peak_dmd_percentage[hour][1];
				//hourly_load.push_back(int(node_peak_load[n] * daily_percentage * hrly_percentage));
				season_load[Season::WinterWkdy].push_back(int(node_peak_load[n] * hrly_peak_dmd_percentage[hour][0]));
				season_load[Season::WinterWknd].push_back(int(node_peak_load[n] * hrly_peak_dmd_percentage[hour][1] * WkndPeakLoadRatio));
			}
			//season_load[Season::Winter] = hourly_load;

			// summer week
			for (int h = 0; h < _model.ucPeriods; h++) {
				const int hour = int(h % HOURSPERDAY);
				season_load[Season::SummerWkdy].push_back(int(node_peak_load[n] * hrly_peak_dmd_percentage[hour][2]));
				season_load[Season::SummerWknd].push_back(int(node_peak_load[n] * hrly_peak_dmd_percentage[hour][3] * WkndPeakLoadRatio));
			}

			// fall and spring week
			for (int h = 0; h < _model.ucPeriods; h++) {
				const int hour = int(h % HOURSPERDAY);
				season_load[Season::SpringWkdy].push_back(int(node_peak_load[n] * hrly_peak_dmd_percentage[hour][4]));
				season_load[Season::SpringWknd].push_back(int(node_peak_load[n] * hrly_peak_dmd_percentage[hour][5] * WkndPeakLoadRatio));
				season_load[Season::FallWkdy].push_back(int(node_peak_load[n] * hrly_peak_dmd_percentage[hour][4]));
				season_load[Season::FallWknd].push_back(int(node_peak_load[n] * hrly_peak_dmd_percentage[hour][5] * WkndPeakLoadRatio));
			}
			bus.season_hourly_loads.push_back(season_load);			
		}
	}
}

void ReadWrite::readGenerators(Model &_model)
{
   const auto paramReg = ParamRegistry::instance();
   const auto numPeriods = paramReg->numPeriods;

	ifstream f_gen(genFile.c_str());												//open the csv file
	if (!f_gen.good())
	{
		std::cerr << "CANNOT OPEN THE GENERATOR FILE!" << endl;
		return;
	}
	int row = 0, col = 0;
	while (f_gen.good()) {
		row++;
		Generator newGen(_model.numNodes);
		string existing_or_potential;
		string s;
		if (!getline(f_gen, s)) break;

		col = 0;
		istringstream ss(s);
		while (ss)
		{
			col++;
			string value;
			if (!getline(ss, value, ',')) break;
			if (row == 1 && col == 2)
			{
				_model.numGen = std::atoi(value.c_str());
			}
			if (row > 2 && row <= 2 + _model.numGen)
			{
				switch (col)
				{
				case 1:
					newGen.id = atoi(value.c_str());
					break;
				case 2:
					newGen.name = value.c_str();
					break;
				case 3:
					newGen.loc = atoi(value.c_str());
					break;
				case 4:
					existing_or_potential = value.c_str();
					break;
				case 6:
					newGen.type = atoi(value.c_str());
					break;
				case 9:
					newGen.capacity = (atof(value.c_str()));
					break;
				case 12:
					newGen.minUpTime = atoi(value.c_str());
					break;
				case 13:
					newGen.minDownTime = atoi(value.c_str());
					break;
				case 14:
					newGen.rampUp = atof(value.c_str());
					break;
				case 15:
					newGen.rampDown = atof(value.c_str());
					break;
				case 21:
					newGen.startUpCost = atof(value.c_str());
					break;
				case 22:
					newGen.shutDownCost = atof(value.c_str());
					break;
				case 23:
					newGen.powerMin = atof(value.c_str());
					break;
				case 24:
					newGen.powerMax = atof(value.c_str());
					break;
				case 25:
					newGen.maxSpinRsv = atof(value.c_str());
					break;
				case 26:
					newGen.fuelCostA = atof(value.c_str());
					break;
				case 27:
					newGen.fuelCostB = atof(value.c_str());
					break;
				case 28:
					newGen.fuelCostC = atof(value.c_str());
					break;
				case 29:
					newGen.varCost = atof(value.c_str()) * 1000.0;
					break;
				default:
					break;
				}
			}

		}
		if (row > 2 && row <= 2 + _model.numGen)
		{
			newGen.currentCapacity = (existing_or_potential == "E") ? newGen.capacity : 0;
			_model.generators.push_back(newGen);
		}
	}
	f_gen.close();
}

void ReadWrite::readRenewable(Model &_model) 
{
	_model.numRnGen = 0;
	int row = 0, col = 0;
	// wind energy
	std::string rnGenFile;
	
	for (auto &gen : _model.generators) {
		if (gen.get_rn_type() == RenewableType::WIND) rnGenFile = inputDirectory + "/renewable/WIND.csv";
		else if(gen.get_rn_type() == RenewableType::SOLAR)	rnGenFile = inputDirectory + "/renewable/SOLAR.csv";
		else continue;
		
		std::ifstream f_rnGenFile(rnGenFile.c_str());
		gen.is_renewable = true;
		gen.minUpTime = 0;
		gen.minDownTime = 0;
		gen.rampUp = gen.capacity;
		gen.rampDown = gen.capacity;
		_model.rn_generators.push_back(&gen);
		_model.numRnGen++;
		gen.distn_params = vector2double(Season::numSeasons);
		int row = 0, col = 0;
		while (f_rnGenFile.good()) {
			row++;
			string s;
			if (!getline(f_rnGenFile, s)) break;
			col = 0;
			istringstream ss(s);
			while (ss)
			{
				col++;
				std::string value;
				if (!getline(ss, value, ',')) break;
				if (row < 2) continue;
				switch (col) {
				case 2:
					gen.distn_params[Season::SpringWkdy].push_back(std::atof(value.c_str()));
					gen.distn_params[Season::SpringWknd].push_back(std::atof(value.c_str()));
					break;
				case 3:
					gen.distn_params[Season::SummerWkdy].push_back(std::atof(value.c_str()));
					gen.distn_params[Season::SummerWknd].push_back(std::atof(value.c_str()));
					break;
				case 4:
					gen.distn_params[Season::FallWkdy].push_back(std::atof(value.c_str()));
					gen.distn_params[Season::FallWknd].push_back(std::atof(value.c_str()));
					break;
				case 5:
					gen.distn_params[Season::WinterWkdy].push_back(std::atof(value.c_str()));
					gen.distn_params[Season::WinterWknd].push_back(std::atof(value.c_str()));
					break;
				default:
					break;
				}
			}
		}
	}
}

void ReadWrite::readLines(Model &_model)
{
   const auto paramReg = ParamRegistry::instance();
	int col = 0;
	int row = 0;
	int _id = 0;
	ifstream f_lines(lineFile.c_str());//open the csv file
	if (!f_lines.good())
	{
		cerr << "CANNOT OPEN THE TRANSMISSION-LINE FILE!" << endl;
		return;
	}
	row = 0; col = 0;
	while (f_lines.good())
	{
		row++;
		Line newLine(paramReg->numScenarios, paramReg->numPeriods);
		string existing_or_potential;
		string s;
		if (!getline(f_lines, s)) break;

		col = 0;
		istringstream ss(s);
		while (ss)
		{
			col++;
			string value;
			if (!getline(ss, value, ',')) break;
			if (row == 1 && col == 2)
			{
				_model.numLine = atoi(value.c_str());
			}
			if (row > 2 && row <= 2 + _model.numLine)
			{
				switch (col)
				{
				case 1:
					newLine.from = atoi(value.c_str());
					break;
				case 2:
					break;
				case 3:
					newLine.to = atoi(value.c_str());
					break;
				case 5:
					existing_or_potential = value.c_str();
					break;
				case 8:
					newLine.capacity = atof(value.c_str());
					break;
				case 9:
					newLine.volt = atof(value.c_str());
					break;
				case 10:
					newLine.imped = atof(value.c_str());
					break;
				case 11:
					newLine.varCost = atof(value.c_str()) * 1000.0;
					break;
            case 12:
               newLine.length = atof(value.c_str());
               break;
				default:
					break;
				}
			}
		}
		if (row > 2 && row <= 2 + _model.numLine)
		{
			newLine.id = _id++;
			newLine.lossRatio = paramReg->lossPower * newLine.length;
			newLine.currentCapacity = (existing_or_potential == "E") ? newLine.capacity : 0;
			_model.lines.push_back(newLine);
		}
	}
	f_lines.close();
}

void ReadWrite::readStorages(Model &_model)
{
   const auto paramReg = ParamRegistry::instance();
	int col = 0;
	int row = 0;
	ifstream f_stg(stgFile.c_str());//open the csv file
	if (!f_stg.good())
	{
		cerr << "CANNOT OPEN THE Storage FILE!" << endl;
		return;
	}
	row = 0; col = 0;
	while (f_stg.good())
	{
		row++;
		Storage newStg(paramReg->numScenarios, paramReg->numPeriods);
		string existing_or_potential;
		string s;
		if (!getline(f_stg, s)) break;

		col = 0;
		istringstream ss(s);
		while (ss)
		{
			col++;
			string value;
			if (!getline(ss, value, ',')) break;
			if (row == 1 && col == 2)
			{
				_model.numStg = atoi(value.c_str());
			}
			if (row > 2 && row <= 2 + _model.numStg)
			{
				switch (col)
				{
				case 1:
					newStg.id = atoi(value.c_str());
					break;
				case 3:
					newStg.loc = atoi(value.c_str());
					break;
				case 4:
					existing_or_potential = value.c_str();
				case 5:
					newStg.capacity = atof(value.c_str());
					break;
				case 6:
					newStg.varCost = atof(value.c_str()) * 1000.0;
					break;
				default:
					break;
				}
			}
		}
		if (row > 2 && row <= 2 + _model.numStg)
		{
			newStg.currentCapacity = (existing_or_potential == "E") ? newStg.capacity : 0;
			_model.storages.push_back(newStg);
		}
	}
	f_stg.close();
}

void ReadWrite::dataGeneration(Model &_model)
{
   const auto paramReg = ParamRegistry::instance();

   const auto numNodes = _model.numNodes;
   const auto numDmdScenarios = paramReg->numDmdScenarios;
   const auto numGen = _model.numGen;
   const auto numLine = _model.numLine;
   const auto numStg = _model.numStg;
   const auto numBus = _model.numBus;
   const auto numRnGen = _model.numRnGen;
   const auto ucPeriods = _model.ucPeriods;
	int index = 0;
	std::default_random_engine generator;

	/* ************************** 1. Generate upper level scenario data ***************************** */
	for (int n = 0; n < numNodes; n++) {
	   const auto &node = _model.scenarioNodes[n];
	   for (auto &gen : _model.generators) {
	      const auto factor = gen.is_renewable ? paramReg->renewableCapex[n % paramReg->numRealizations] :
	            paramReg->fossilCapex[n % paramReg->numRealizations];
	      gen.cost[n] = (n == 0) ? int(gen.capacity * gen.varCost) : int(gen.cost[node.parent] * factor);

	      const auto fuelPriceFactor = paramReg->fuelPriceChange[n % paramReg->numRealizations];
	      gen.fuelPriceRate[n] = ((n == 0) || (gen.varCost != 900.0 * 1000.0)) ? 1.0 :
	            (gen.fuelPriceRate[node.parent] * fuelPriceFactor);   //TODO: generator type distinguished between coal and natural gas
	   }
	   for (auto &line : _model.lines) {
	      const auto factor = paramReg->lineCapex[n%paramReg->numRealizations];
	      line.cost[n] = (n == 0) ? int(line.capacity * line.varCost) : int(line.cost[node.parent] * factor);
	   }
	   for (auto &stg : _model.storages) {
	      const auto factor = paramReg->storageCapex[n%paramReg->numRealizations];
	      stg.cost[n] = (n == 0) ? int(stg.capacity * stg.varCost) : int(stg.cost[node.parent] * factor);
	   }
	}

	/* ************************** 2. Generate lower level load scenario data ***************************** */
	int * season_hourly_load_scenarios = new int[numNodes * Season::numSeasons * numDmdScenarios * ucPeriods * numBus];
	if (_model.is_master_proc()) 
	{	
	   index = 0;
		for (int n = 0; n < numNodes; n++) {
			/* ************************** Generate lower level scenario data ************************** */
			for (int season = 0; season != Season::numSeasons; season++) {
			   for (int d = 0; d < paramReg->numDmdScenarios; d++) {
					for (int h = 0; h < ucPeriods; h++) {
						for (auto & bus : _model.buses) {
						   if(paramReg->numDmdScenarios == 1){
						      season_hourly_load_scenarios[index++] = bus.season_hourly_loads[n][season][h];
						   }else{
						      if (bus.season_hourly_loads[n][season][h] == 0) {
						         season_hourly_load_scenarios[index++] = 0;
						      }
						      else {
						         const int mean_peak_load = bus.season_hourly_loads[n][season][h];
						         std::normal_distribution<double> distribution(double(mean_peak_load), stddev * mean_peak_load);
						         double sample = 0.;
						         do{
						            sample = distribution(generator);
						         }while((sample < mean_peak_load * (1 - stddev)) || (sample > mean_peak_load * (1 + stddev)));
						         season_hourly_load_scenarios[index++] = int(sample);
						      }
						   }
						}
					}
				}
			}			
		}		
	}
	MPI_Bcast(season_hourly_load_scenarios, numNodes * Season::numSeasons * numDmdScenarios * ucPeriods * numBus, MPI_INT, _model.get_master(), MPI_COMM_WORLD);
	index = 0;
	for (int n = 0; n < numNodes; n++) {
		/* ************************** Generate lower level scenario data ************************** */
		for (int season = 0; season != Season::numSeasons; season++) {
		   for (int d = 0; d < paramReg->numDmdScenarios; d++) {
				for (int h = 0; h < ucPeriods; h++) {
					for (auto & bus : _model.buses) {
						bus.season_hourly_load_scenarios.push_back(season_hourly_load_scenarios[index++]);
					}
				}
			}
		}
	}

	delete[] season_hourly_load_scenarios;

	/* ************************** 3. Generate uncertain supply data for renewable generators ***************************** */
	int * uncertain_supply = new int[Season::numSeasons * numDmdScenarios * ucPeriods * numRnGen];
	for (auto &gen : _model.generators) {
		gen.supply_mean = vector2int(Season::numSeasons);
		for (int s = 0; s < Season::numSeasons; s++) {
			gen.supply_mean[s] = std::vector<int>(_model.ucPeriods, 0);
		}
	}
	if (_model.is_master_proc())
	{
		index = 0;
		for (int season = 0; season != Season::numSeasons; season++) {
		   for (int d = 0; d < paramReg->numDmdScenarios; d++) {
				for (int h = 0; h < ucPeriods; h++) {
				   std::map<int, double> busWindSpeed;
				   std::map<int, double> busSolar;
					for (auto & gen : _model.generators) {
						if (!gen.is_renewable) continue;							
						if (gen.get_rn_type() == RenewableType::WIND) {
						   double wind_speed = 0;
						   const auto itrWind = busWindSpeed.find(gen.loc);
						   if(itrWind == busWindSpeed.end()){
						      // wind speed follows weibull distribution
						      std::weibull_distribution<double> distribution(gen.distn_params[season][0], gen.distn_params[season][1]);
						      wind_speed = distribution(generator);
						      busWindSpeed.emplace(gen.loc, wind_speed);
						   }else{
						      wind_speed = itrWind->second;
						   }
							const double avg_wind_speed = gen.distn_params[season][2];
							// power generation is dependent upon wind speed
							int supply = 0;
							if (wind_speed <= POWER_WIND_MIN) {
								supply = 0;
							}
							else if (wind_speed >= POWER_WIND_MAX) {
								supply = int(gen.capacity);
							}
							else {
								const double coef = gen.capacity / (POWER_WIND_MAX - POWER_WIND_MIN);
								supply = int(coef * (wind_speed - POWER_WIND_MIN));
							}
							uncertain_supply[index] = supply;

							int avg_supply = 0;
							if (avg_wind_speed <= POWER_WIND_MIN) {
								avg_supply = 0;
							}
							else if (avg_wind_speed >= POWER_WIND_MAX) {
								avg_supply = int(gen.capacity);
							}
							else {
								const double coef = gen.capacity / (POWER_WIND_MAX - POWER_WIND_MIN);
								avg_supply = int(coef * (avg_wind_speed - POWER_WIND_MIN));
							}
							if (d == 0) gen.supply_mean[season][h] = avg_supply;
						}
						else {
							const double supply_ratio = gen.distn_params[season][h%HOURSPERDAY];
							if (supply_ratio == 0.0) {
								uncertain_supply[index] = 0;
							}
							else {
							   double solar_power = 0.0;
							   const auto itrSolar = busSolar.find(gen.loc);
							   if(itrSolar == busSolar.end()){
							      std::normal_distribution<double> distribution(supply_ratio * gen.capacity, stddev * supply_ratio * gen.capacity);
							      do{
							         solar_power = distribution(generator);
							      }while((solar_power < supply_ratio * gen.capacity * (1 - stddev)) || (solar_power > supply_ratio * gen.capacity * (1 + stddev)));
							      busSolar.emplace(gen.loc, solar_power);
							   }else{
							      solar_power = itrSolar->second;
							   }

							   uncertain_supply[index] = std::min<int>(int(solar_power), int(gen.capacity));
							}
							if(d == 0) gen.supply_mean[season][h] = int(supply_ratio * gen.capacity);
						}
						if(paramReg->numDmdScenarios == 1){
						   uncertain_supply[index] = gen.supply_mean[season][h];
						}
						index++;
					}
				}
			}
		}
	}
	MPI_Bcast(uncertain_supply, Season::numSeasons * numDmdScenarios * ucPeriods * numRnGen, MPI_INT, _model.get_master(), MPI_COMM_WORLD);	
	for (int n = 0; n < numNodes; n++) {	
		index = 0;
		for (int season = 0; season != Season::numSeasons; season++) {
		   for (int d = 0; d < paramReg->numDmdScenarios; d++) {
				for (int h = 0; h < ucPeriods; h++) {
					for (auto & gen : _model.generators) {
						if (!gen.is_renewable) continue;
						gen.supply_scenarios.push_back(uncertain_supply[index++]);
					}
				}
			}
		}
	}
	delete[] uncertain_supply;
}

// write bounds
void ReadWrite::printBounds(int iteration, double ub, double* reducedCost, int numNodes, double lb, double gap)
{
	if (boundsFile.rfind(outputDirectory, 0) != 0) {
		boundsFile = outputDirectory + boundsFile;
	}
	std::ifstream bounds_file(boundsFile);
	std::fstream fbounds;

	if (bounds_file && reducedCost != nullptr)
	{
		fbounds.open(boundsFile.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		fbounds << iteration << "\t" << ub << "\t";
		for (int n = 0; n < numNodes; n++) {
			fbounds << reducedCost[n] << "\t";
		}
		fbounds << lb << "\t" << gap << std::endl;
	}
	else {
		fbounds.open(boundsFile.c_str(), std::fstream::in | std::fstream::out | std::fstream::trunc);
		fbounds << "Itn\tUB\t";
		for (int n = 0; n < numNodes; n++) {
			fbounds << n << "\t";
		}
		fbounds << "LB\tGap";
		fbounds << std::endl;
	}

	fbounds.close();
}

void ReadWrite::printSolution(Model &_model)
{
	if (!_model.is_master_proc()) return;
	PRINT_SECTION("Printing The Optimal Solution");
#ifdef LPMODEL
	char buf[100];
	sprintf(buf, "%s/direct.lp", outputDirectory.c_str());
	_model.MPmodel.write(buf);
#endif

	const auto paramReg = ParamRegistry::instance();
	char filename[133];
	char fullpath[177];

	cout << "Computation time <seconds>: -> " << _model.cmp_time << endl;

	switch (_model.get_algorithm())
	{
	case Algorithm::DIRECT:
	   sprintf(filename, "DIRECT_%dnodes_%dperiods_%ddemand", _model.numNodes, paramReg->numPeriods, paramReg->numDmdScenarios);
		break;
	case Algorithm::EXPECT:
	   sprintf(filename, "EXPECT_%dnodes_%dperiods_%ddemand", _model.numNodes, paramReg->numPeriods, paramReg->numDmdScenarios);
		break;
	case Algorithm::COLGEN:
	   sprintf(filename, "COLGEN_%dnodes_%dperiods_%ddemand", _model.numNodes, paramReg->numPeriods, paramReg->numDmdScenarios);
		break;
	case Algorithm::BENDERS:
	   sprintf(filename, "BENDERS_%dnodes_%dperiods_%ddemand", _model.numNodes, paramReg->numPeriods, paramReg->numDmdScenarios);
		break;
	case Algorithm::NESTED:
	   sprintf(filename, "NESTED_%dnodes_%dperiods_%ddemand", _model.numNodes, paramReg->numPeriods, paramReg->numDmdScenarios);
		break;
	default:
		break;
	}

	sprintf(fullpath, "%s/%s.csv", outputDirectory.c_str(), filename);
	remove(fullpath);
	FILE * outFile;
	outFile = fopen(fullpath, "a");
	fclose(outFile);
	fstream output(fullpath, fstream::in | fstream::out);
	output << "Total Run Time <seconds>:," << _model.cmp_time << std::endl;
	output << "Serial Run Time <seconds>:," << _model.totalSolutionTime << std::endl;
	output << "Parallel Computing Efficiency <%>," << _model.totalSolutionTime / (_model.cmp_time * _model.num_ranks)  * 100 << std::endl;
	output << "No. Iterations," << _model.iteration << std::endl;
	output << "Solution information:, " << endl;

	if(_model.MPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE || _model.MPmodel.get(GRB_IntAttr_Status) == GRB_UNBOUNDED)
		return;

	output << "IP Objective:," << _model.MpUB << endl;
	output << "LP objective:," << _model.MpLpObj << endl;
	output << "IP gap:," << double((_model.MpUB - _model.MpLpObj) / _model.MpLpObj);
	output << "Investment:," << _model.investmentCost << endl;
	output << "Operation:,";
	for (int n = 0; n < _model.numNodes; n++) {
		output << _model.operationalCost[n] << ",";
	}
	output << endl;
	for (int g = 0; g < _model.numGen; g++)
	{
		if (_model.generators[g].exists()) continue;
		output << _model.generators[g].name << "@bus" << _model.generators[g].loc << ",";
		for (int t = 0; t < paramReg->numPeriods; t++) {
			output << t << ",";		
		}
		output << endl;
		for (int s = 0; s < paramReg->numScenarios; s++) {
			output << s << ",";
			for (int t = 0; t < paramReg->numPeriods; t++) {
				const auto node = _model.getNode(s, t);
				output << _model.sol_xg[node->id][g] << ",";
			}
			output << endl;
		}
	}
	for (int l = 0; l < _model.numLine; l++)
	{
		if (_model.lines[l].exists()) continue;
		output << "Line" << _model.lines[l].id << "(" << _model.lines[l].from << "-" << _model.lines[l].to << "),";
		for (int t = 0; t < paramReg->numPeriods; t++) {
			output << t << ",";
		}
		output << endl;
		for (int s = 0; s < paramReg->numScenarios; s++) {
			output << s << ",";
			for (int t = 0; t < paramReg->numPeriods; t++) {
				const auto node = _model.getNode(s, t);
				output << _model.sol_xl[node->id][l] << ",";
			}
			output << std::endl;
		}
	}
	for (int k = 0; k < _model.numStg; k++)
	{
		if (_model.storages[k].exists()) continue;
		output << "Stg" << _model.storages[k].id << "@" << _model.storages[k].loc;
		for (int t = 0; t < paramReg->numPeriods; t++) {
			output << t << ",";
		}
		output << std::endl;
		for (int s = 0; s < paramReg->numScenarios; s++) {
			output << s << ",";
			for (int t = 0; t < paramReg->numPeriods; t++) {
				const auto node = _model.getNode(s, t);
				output << _model.sol_xs[node->id][k] << ",";
			}
			output << std::endl;
		}
	}
	output.close();
}

void ReadWrite::writeGeneratedData(Model &_model)
{
	if (!_model.is_master_proc())
		return;
	char fullpath[177];

	/* **************** Write out cost scenario at upper level ************************** */
	sprintf(fullpath, "%s/scenarios/cost.csv", outputDirectory.c_str());
	remove(fullpath);
	FILE * outFile;
	outFile = fopen(fullpath, "a");
	fclose(outFile);
	fstream output(fullpath, fstream::in | fstream::out);
	output << "ScenarioNode,";
	for (const auto &gen : _model.generators) {
		output << "Gen" << gen.id << ",";
	}
	for (const auto &line : _model.lines) {
		output << "Line" << line.id << ",";
	}
	for (const auto &stg : _model.storages) {
		output << "Stg" << stg.id << ",";
	}
	output << std::endl;

	for (const auto &node : _model.scenarioNodes) {
		output << node.id << ",";
		for (const auto &gen : _model.generators) {
			output << gen.cost[node.id] << ",";
		}
		for (const auto &line : _model.lines) {
			output << line.cost[node.id] << ",";
		}
		for (const auto &stg : _model.storages) {
			output << stg.cost[node.id] << ",";
		}
		output << std::endl;
	}
	/* **************** Write out demand scenario at lower level for node 0 ************************** */
	for (const auto &bus : _model.buses) {
		std::string str_season = "";
		for (int season = 0; season != Season::numSeasons; season++) {
			switch (season) {
			case Season::SpringWkdy:
				str_season = "SpringWkdy";
				break;
			case Season::SummerWkdy:
				str_season = "SummerWkdy";
				break;
			case Season::FallWkdy:
				str_season = "FallWkdy";
				break;
			case Season::WinterWkdy:
				str_season = "WinterWkdy";
				break;
			case Season::SpringWknd:
				str_season = "SpringWknd";
				break;
			case Season::SummerWknd:
				str_season = "SummerWknd";
				break;
			case Season::FallWknd:
				str_season = "FallWknd";
				break;
			case Season::WinterWknd:
				str_season = "WinterWknd";
				break;
			default:
				break;
			}

			sprintf(fullpath, "%s/scenarios/bus%d_%s.csv", outputDirectory.c_str(), bus.id, str_season.c_str());
			remove(fullpath);
			FILE * outFile;
			outFile = fopen(fullpath, "a");
			fclose(outFile);
			fstream output(fullpath, fstream::in | fstream::out);
			for (int d = 0; d < ParamRegistry::instance()->numDmdScenarios; d++) {
				for (int h = 0; h < _model.ucPeriods; h++) {
					output << bus.season_hourly_load_scenarios[_model.get2ndStageIndex(0, season, d, h)] << ",";
				}
				output << std::endl;
			}
		}
	}
}

/* Utility Functions */
std::pair<double, double> ReadWrite::split_range(const string & s) {
	const char delim = '-';
	std::string left, right;
	auto pos = s.find(delim);
	if(pos != string::npos){
	   left = s.substr(0, pos++);
	   right = s.substr(pos, s.length());
	}else{
	   left = s;
	   right = s;
	}

	return std::make_pair(atof(left.c_str()), atof(right.c_str()));
}

template <typename T>
T ReadWrite::NumGen(T minNum, T maxNum)
//generate a number between minNum and maxNum (including both ends)
{
   if(minNum == maxNum) return T(minNum);

   std::random_device rd;
   std::mt19937 mt(rd());
   std::uniform_real_distribution<double> dist(minNum, maxNum);

   return T(dist(mt));

}

std::vector<std::string> ReadWrite::split(const std::string &s) const
{
	std::vector<std::string> v;
	std::istringstream iss(s);
	std::string x;
	while (iss >> x) {
		v.push_back(x);
	}
	return v;
}

