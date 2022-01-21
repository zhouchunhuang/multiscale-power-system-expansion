#include "model.h"

Result Model::direct()
{
	PRINT_SECTION("Solving the problem directly with GUROBI");
	MPmodel.set(GRB_IntParam_OutputFlag, 1);
	MPmodel.set(GRB_DoubleParam_MIPGap, paramReg->mp_gap);
	MPmodel.set(GRB_DoubleParam_TimeLimit, paramReg->time_lmt);
	char buf[100];

#ifdef GetModelSize
	// number of variables
	unsigned long numIntVars = 0;
	unsigned long numContVar = 0;
	for (int n = 0; n < numNodes; n++) {
		for (int g = 0; g < numGen; g++) {
			if (generators[g].exists()) continue;
			numIntVars++;
		}
		for (int l = 0; l < numLine; l++) {
			if (lines[l].exists()) continue;
			numIntVars++;
		}
		for (int k = 0; k < numStg; k++) {
			if (storages[k].exists()) continue;
			numIntVars++;
		}
	}
	numIntVars += numNodes * Season::numSeasons * numDmdScenarios * ucPeriods * numGen * 2;
	numContVar += numNodes * Season::numSeasons * numDmdScenarios * ucPeriods * (numGen + numLine * 2 + numStg * 3);
	std::cout << "Number of discrete decision variables: " << numIntVars << std::endl;
	std::cout << "Number of continuous decision variables: " << numContVar << std::endl;

	// number of constraints
	unsigned long  numConstrs = 0;
	// ---- upper level constraints ------ //
	for (int n = 0; n < numNodes; n++) {
		const auto &node = scenarioNodes[n];
		if (node.period != numPeriods - 1) continue;
		// upper level constraints
		for (int g = 0; g < numGen; g++) {
			numConstrs++;
		}
		for (int l = 0; l < numLine; l++) {
			numConstrs++;
		}
		for (int k = 0; k < numStg; k++) {
			numConstrs++;
		}
	}
	// ---- lower level constraints ------ //
	for (int n = 0; n < numNodes; n++) {
		const auto &node = scenarioNodes[n];
		for (int season = 0; season != Season::numSeasons; season++) {
			for (int d = 0; d < numDmdScenarios; d++) {
				for (int h = 0; h < ucPeriods; h++) {
					const int index = getLowerIndex(n, season, d, h);
					for (int g = 0; g < numGen; g++) {
						const auto &gen = generators[g];
						numConstrs++;
						numConstrs++;
						if (gen.is_renewable) {
							numConstrs++;
						}
						numConstrs++;
						// minimum on and off
						if (h) {
							for (int tau = h; tau < ucPeriods && tau < h + generators[g].minUpTime; tau++) {
								numConstrs++;
							}
							for (int tau = h; tau < ucPeriods && tau < h + generators[g].minDownTime; tau++) {
								numConstrs++;
							}
						}
						// ramp up and down						
						if (h) {
							numConstrs++;
							numConstrs++;
						}
					}
					for (int l = 0; l < numLine; l++) {
						numConstrs++;
					}
					for (int k = 0; k < numStg; k++) {
						numConstrs++;
						numConstrs++;														//(21)
						numConstrs++;
					}
					for (int b = 0; b < numBus; b++) {
						numConstrs++;
					}
				}
			}
		}
	}

	std::cout << "Number of constraints: " << numConstrs << std::endl;
	exit(EXIT_SUCCESS);
#endif

	try {
		PRINT_SUBSECTION("initializing decision variables");
		/* *********************** Upper level expansion decisions *********************** */
		var_xg = new GRBVar*[numNodes];
		var_xl = new GRBVar*[numNodes];
		var_xs = new GRBVar*[numNodes];
		for (int n = 0; n < numNodes; n++) {
			var_xg[n] = MPmodel.addVars(numGen, GRB_BINARY);
			var_xl[n] = MPmodel.addVars(numLine, GRB_BINARY);
			var_xs[n] = MPmodel.addVars(numStg, GRB_BINARY);
#ifdef LPMODEL	
			for (int g = 0; g < numGen; g++) {
				std::sprintf(buf, "var_xg(%d,%d)", n, g);
				var_xg[n][g].set(GRB_StringAttr_VarName, buf);
			}
			for (int l = 0; l < numLine; l++) {
				std::sprintf(buf, "var_xl(%d,%d)", n, l);
				var_xl[n][l].set(GRB_StringAttr_VarName, buf);
			}
			for (int k = 0; k < numStg; k++) {
				std::sprintf(buf, "var_xs(%d,%d)", n, k);
				var_xs[n][k].set(GRB_StringAttr_VarName, buf);
			}
#endif
		}

		/* *********************** Lower level operational decisions: two stage UC problem *********************** */
		const int numLowerlevelIndices = numNodes * Season::numSeasons * paramReg->numDmdScenarios * ucPeriods;
		varGenStatus = new GRBVar *[numLowerlevelIndices];											//binary variable for on/off status of plants alpha[n][d][k][g]
		varStartUp = new GRBVar *[numLowerlevelIndices];												//binary variable for start-up/shut-down of plants alpha[n][d][k][g]
		varPowerGen = new GRBVar *[numLowerlevelIndices];												//amount of power generated from generator p[n][d][k][g]
		varPowerFlow = new GRBVar **[numLowerlevelIndices];												//power flow between bus i and bus j f[n][d][k][l]
		varPowerWithdrawal = new GRBVar *[numLowerlevelIndices];										//total power withdrawn from storage facilities of bus u[n][d][k][sg]
		varPowerInject = new GRBVar *[numLowerlevelIndices];											//total power injected into storage facilities of bus v[n][d][k][sg]
		varPowerRemaining = new GRBVar *[numLowerlevelIndices];										//total remaining power in storage facilities of bus at the beginning of k v[n][d][k][sg]
		for (int n = 0; n < numNodes; n++) {
			for (int season = 0; season != Season::numSeasons; season++) {
			   for (int d = 0; d < paramReg->numDmdScenarios; d++) {
					for (int h = 0; h < ucPeriods; h++) {
						const int index = get2ndStageIndex(n, season, d, h);
						varGenStatus[index] = MPmodel.addVars(numGen, GRB_BINARY);
						varStartUp[index] = MPmodel.addVars(numGen, GRB_BINARY);
						varPowerGen[index] = MPmodel.addVars(numGen);
						varPowerFlow[index] = new GRBVar *[numLine];
						varPowerWithdrawal[index] = MPmodel.addVars(numStg);
						varPowerInject[index] = MPmodel.addVars(numStg);
						varPowerRemaining[index] = MPmodel.addVars(numStg);
						for (int l = 0; l < numLine; l++) {
							varPowerFlow[index][l] = MPmodel.addVars(Direction::numDirections);
						}
#ifdef LPMODEL
						for (int g = 0; g < numGen; g++) {
							std::sprintf(buf, "varGenStatus(%d,%d,%d,%d,%d)", n, season, d, h, g);
							varGenStatus[index][g].set(GRB_StringAttr_VarName, buf);
							std::sprintf(buf, "varStartUp(%d,%d,%d,%d,%d)", n, season, d, h, g);
							varStartUp[index][g].set(GRB_StringAttr_VarName, buf);
							std::sprintf(buf, "varPowerGen(%d,%d,%d,%d,%d)", n, season, d, h, g);
							varPowerGen[index][g].set(GRB_StringAttr_VarName, buf);
						}
						for (int l = 0; l < numLine; l++) {
							std::sprintf(buf, "varPowerFlow(%d,%d,%d,%d,%d,+)", n, season, d, h, l);
							varPowerFlow[index][l][Direction::Forward].set(GRB_StringAttr_VarName, buf);
							std::sprintf(buf, "varPowerFlow(%d,%d,%d,%d,%d,-)", n, season, d, h, l);
							varPowerFlow[index][l][Direction::Backward].set(GRB_StringAttr_VarName, buf);
						}
						for (int k = 0; k < numStg; k++) {
							std::sprintf(buf, "varPowerWithdrawal(%d,%d,%d,%d,%d)", n, season, d, h, k);
							varPowerWithdrawal[index][k].set(GRB_StringAttr_VarName, buf);
							std::sprintf(buf, "varPowerInject(%d,%d,%d,%d,%d)", n, season, d, h, k);
							varPowerInject[index][k].set(GRB_StringAttr_VarName, buf);
							std::sprintf(buf, "varPowerRemaining(%d,%d,%d,%d,%d)", n, season, d, h, k);
							varPowerRemaining[index][k].set(GRB_StringAttr_VarName, buf);
						}
#endif
					}
				}
			}
		}

		PRINT_SUBSECTION("initializing objective function");
		/* *********************** Objective *********************** */
		GRBLinExpr xpr = 0;
		for (int n = 0; n < numNodes; n++) {
			const auto &node = scenarioNodes[n];
			for (int g = 0; g < numGen; g++) {
				const auto &gen = generators[g];
				xpr += node.probability * gen.cost[n] * var_xg[n][g];
			}
			for (int l = 0; l < numLine; l++) {
				const auto &line = lines[l];
				xpr += node.probability * line.cost[n] * var_xl[n][l];
			}
			for (int k = 0; k < numStg; k++) {
				const auto &stg = storages[k];
				xpr += node.probability * stg.cost[n] * var_xs[n][k];
			}
			for (int season = 0; season != Season::numSeasons; season++) {
				for (int h = 0; h < ucPeriods; h++) {
					for (int g = 0; g < numGen; g++) {
						const auto &gen = generators[g];
						for (int d = 0; d < paramReg->numDmdScenarios; d++) {
						   xpr += node.probability * paramReg->prob_dmd_scenario* get_opn_cost_factor(season) * (gen.startUpCost * varStartUp[getLowerIndex(n, season, d, h)][g]
								+ gen.fuelCostA * varGenStatus[getLowerIndex(n, season, d, h)][g]
								+ gen.fuelCostB * gen.fuelPriceRate[n] * varPowerGen[getLowerIndex(n, season, d, h)][g]);
						}
					}
				}
			}
		}
		MPmodel.setObjective(xpr, GRB_MINIMIZE);

		PRINT_SUBSECTION("initializing constraints");
		/* *********************** Constraints *********************** */
		// ---- upper level constraints ------ //
		for (int n = 0; n < numNodes; n++) {
			const auto &node = scenarioNodes[n];
			if (node.period != paramReg->numPeriods - 1) continue;
			// upper level constraints
			for (int g = 0; g < numGen; g++) {
				GRBLinExpr tempExpr = 0;
				for (const int i : node.nodesInPath) {
					tempExpr += var_xg[i][g];
				}
				MPmodel.addConstr(tempExpr, GRB_LESS_EQUAL, 1 - static_cast<int>(generators[g].exists()), "UpperLevel");
			}
			for (int l = 0; l < numLine; l++) {
				GRBLinExpr tempExpr = 0;
				for (const int i : node.nodesInPath) {
					tempExpr += var_xl[i][l];
				}
				MPmodel.addConstr(tempExpr, GRB_LESS_EQUAL, 1 - static_cast<int>(lines[l].exists()), "UpperLevel");
			}
			for (int k = 0; k < numStg; k++) {
				GRBLinExpr tempExpr = 0;
				for (const int i : node.nodesInPath) {
					tempExpr += var_xs[i][k];
				}
				MPmodel.addConstr(tempExpr, GRB_LESS_EQUAL, 1 - static_cast<int>(storages[k].exists()), "UpperLevel");
			}
		}

		// ---- lower level constraints ------ //
		for (int n = 0; n < numNodes; n++) {
			const auto &node = scenarioNodes[n];
			for (int season = 0; season != Season::numSeasons; season++) {
			   for (int d = 0; d < paramReg->numDmdScenarios; d++) {
					for (int h = 0; h < ucPeriods; h++) {
						const int index = getLowerIndex(n, season, d, h);
						for (int g = 0; g < numGen; g++) {
							const auto &gen = generators[g];

							GRBLinExpr tempExpr = int(generators[g].exists());
							for (const int i : node.nodesInPath) {
								tempExpr += var_xg[i][g];
							}
							std::sprintf(buf, "GenExist(%d,%d,%d,%d)", n, season, h, g);																							// Coupling constraint			
							MPmodel.addConstr(varGenStatus[index][g], GRB_LESS_EQUAL, tempExpr, buf);

							std::sprintf(buf, "PowUB(%d,%d,%d,%d,%d)", n, season, d, h, g);
							MPmodel.addConstr(varPowerGen[index][g] - gen.capacity * varGenStatus[index][g], GRB_LESS_EQUAL, 0, buf);

							if (gen.is_renewable) {
								std::sprintf(buf, "RnGen(%d,%d,%d,%d,%d)", n, season, d, h, g);
								MPmodel.addConstr(varPowerGen[index][g], GRB_LESS_EQUAL, gen.supply_scenarios[index], buf);
							}

							std::sprintf(buf, "GenAct(%d,%d,%d,%d)", n, season, h, g);
							if (h == 0) {																														// weekly UC cconstraints
								MPmodel.addConstr(varStartUp[index][g] - varGenStatus[index][g], GRB_GREATER_EQUAL, 0, buf);
							}
							else
							{
								MPmodel.addConstr(varStartUp[index][g] - varGenStatus[index][g] + varGenStatus[getLowerIndex(n, season, d, h - 1)][g], GRB_GREATER_EQUAL, 0, buf);
							}

							// minimum on and off
							if (h) {
								for (int tau = h; tau < ucPeriods && tau < h + generators[g].minUpTime; tau++) {
									std::sprintf(buf, "MinUp(%d,%d,%d,%d,%d)", n, season, h, tau, g);
									MPmodel.addConstr(varGenStatus[index][g] - varGenStatus[getLowerIndex(n, season, d, h - 1)][g],
										GRB_LESS_EQUAL, varGenStatus[getLowerIndex(n, season, d, tau)][g]);
								}
								for (int tau = h; tau < ucPeriods && tau < h + generators[g].minDownTime; tau++) {
									std::sprintf(buf, "MinDown(%d,%d,%d,%d,%d)", n, season, h, tau, g);
									MPmodel.addConstr(varGenStatus[index][g] - varGenStatus[getLowerIndex(n, season, d, h - 1)][g],
										GRB_GREATER_EQUAL, varGenStatus[getLowerIndex(n, season, d, tau)][g] - 1);
								}
							}

							// ramp up and down						
							if (h) {
								std::sprintf(buf, "RampUp(%d,%d,%d,%d,%d)", n, season, d, h, g);
								MPmodel.addConstr(varPowerGen[index][g] - varPowerGen[getLowerIndex(n, season, d, h - 1)][g], GRB_LESS_EQUAL, gen.rampUp, buf);
								std::sprintf(buf, "RampDown(%d,%d,%d,%d,%d)", n, season, d, h, g);
								MPmodel.addConstr(varPowerGen[getLowerIndex(n, season, d, h - 1)][g] - varPowerGen[index][g], GRB_LESS_EQUAL, gen.rampDown, buf);
							}
						}
						for (int l = 0; l < numLine; l++) {
							const auto &line = lines[l];
							GRBLinExpr tempExpr = 0;
							if (!line.exists()) {
								for (const int i : node.nodesInPath) {
									tempExpr += var_xl[i][l];
								}
							}
							tempExpr *= line.capacity;
							std::sprintf(buf, "Lne(%d,%d,%d,%d,%d)", n, season, d, h, l);
							MPmodel.addConstr(varPowerFlow[index][l][Direction::Forward] + varPowerFlow[index][l][Direction::Backward] - tempExpr,
								GRB_LESS_EQUAL, line.currentCapacity, buf);
						}
						for (int k = 0; k < numStg; k++) {
							const auto &storage = storages[k];
							GRBLinExpr tempExpr = 0;
							if (!storage.exists()) {
								for (const int i : node.nodesInPath) {
									tempExpr += var_xs[i][k];
								}
							}
							tempExpr *= storage.capacity;
							std::sprintf(buf, "Stg(%d,%d,%d,%d,%d)", n, season, d, h, k);																							//Constraint 23
							MPmodel.addConstr(varPowerRemaining[index][k] - tempExpr, GRB_LESS_EQUAL, storage.currentCapacity, buf);

							std::sprintf(buf, "StgLmt(%d,%d,%d,%d,%d)", n, season, d, h, k);
							MPmodel.addConstr(varPowerWithdrawal[index][k] - varPowerRemaining[index][k], GRB_LESS_EQUAL, 0, buf);														//(21)
							// std::sprintf(buf, "StgLmt(%d,%d,%d,%d)", n, d, h, k);
							// MPmodel.addConstr(varPowerRemaining[n][d][h][k], GRB_LESS_EQUAL, storage.capacity, buf);
							std::sprintf(buf, "StgBal(%d,%d,%d,%d,%d)", n, season, d, h, k);
							if (h == 0) {																														//in the 1st lower level period, reset the storage facility
								MPmodel.addConstr(varPowerRemaining[index][k], GRB_EQUAL, 0, buf);
							}
							else {
								MPmodel.addConstr(varPowerRemaining[index][k] - varPowerRemaining[get2ndStageIndex(n, season, d, h - 1)][k]
									- varPowerInject[get2ndStageIndex(n, season, d, h - 1)][k] + varPowerWithdrawal[get2ndStageIndex(n, season, d, h - 1)][k], GRB_EQUAL, 0, buf);							//Constraint 22
							}
						}
						for (int b = 0; b < numBus; b++) {
							std::sprintf(buf, "Dmd(%d,%d,%d,%d,%d)", n, season, d, h, b);
							GRBLinExpr tempExpr = 0;
							for (int l = 0; l < numLine; l++) {
							   if (lines[l].from == b)	tempExpr -= varPowerFlow[index][l][Direction::Forward] - varPowerFlow[index][l][Direction::Backward] * (1 - lines[l].lossRatio);
							   if (lines[l].to == b) tempExpr += varPowerFlow[index][l][Direction::Forward] * (1 - lines[l].lossRatio) - varPowerFlow[index][l][Direction::Backward];
							}
							for (int g = 0; g < numGen; g++) {
								if (generators[g].loc == b) tempExpr += varPowerGen[index][g];
							}
							for (int k = 0; k < numStg; k++) {
								if (storages[k].loc == b) tempExpr += varPowerWithdrawal[index][k] * paramReg->stgCeff - varPowerInject[index][k];
							}
							MPmodel.addConstr(tempExpr, GRB_GREATER_EQUAL, buses[b].season_hourly_load_scenarios[index], buf);
						}
					}
				}
			}
		}

		PRINT_SUBSECTION("optimizing");
		clock_t _start = clock();
		MPmodel.optimize();
		clock_t _end = clock();
		cmp_time = double((_end - _start)) / CLOCKS_PER_SEC;

		if (MPmodel.get(GRB_IntAttr_Status) != GRB_OPTIMAL)
		{
			return static_cast<Result>(MPmodel.get(GRB_IntAttr_Status));
		}

		MpUB = MPmodel.get(GRB_DoubleAttr_ObjVal);
		getExpansionSolution();

		operationalCost = new double[numNodes];
		for (int n = 0; n < numNodes; n++) {
			const auto &node = scenarioNodes[n];
			operationalCost[n] = 0;
			for (int season = 0; season != Season::numSeasons; season++) {
				for (int h = 0; h < ucPeriods; h++) {
					for (int g = 0; g < numGen; g++) {
						const auto &gen = generators[g];
						operationalCost[n] += node.probability * get_opn_cost_factor(season) * (gen.startUpCost * varStartUp[get1stStageIndex(n, season, h)][g].get(GRB_DoubleAttr_X) +
							gen.fuelCostA * varGenStatus[get1stStageIndex(n, season, h)][g].get(GRB_DoubleAttr_X));
						for (int d = 0; d < paramReg->numDmdScenarios; d++) {
						   operationalCost[n] += node.probability * paramReg->prob_dmd_scenario * get_opn_cost_factor(season) * gen.fuelCostB * gen.fuelPriceRate[n]* varPowerGen[get2ndStageIndex(n, season, d, h)][g].get(GRB_DoubleAttr_X);
						}
					}
				}
			}
		}
	}
	catch (GRBException e) {
		std::cerr << e.getMessage() << std::endl;
		exit(EXIT_FAILURE);
	}

	return Result::OPTIMAL;
}

Result Model::expectedValForm()
{
	PRINT_SECTION("Solving the counterpart with deterministic lower level with GUROBI");
	MPmodel.set(GRB_DoubleParam_MIPGap, paramReg->mp_gap);
	MPmodel.set(GRB_DoubleParam_TimeLimit, paramReg->time_lmt);
	char buf[100];
	/* *********************** Upper level expansion decisions *********************** */
	var_xg = new GRBVar*[numNodes];
	var_xl = new GRBVar*[numNodes];
	var_xs = new GRBVar*[numNodes];
	for (int n = 0; n < numNodes; n++) {
		var_xg[n] = MPmodel.addVars(numGen, GRB_BINARY);
		var_xl[n] = MPmodel.addVars(numLine, GRB_BINARY);
		var_xs[n] = MPmodel.addVars(numStg, GRB_BINARY);
#ifdef LPMODEL	
		for (int g = 0; g < numGen; g++) {
			std::sprintf(buf, "var_xg(%d,%d)", n, g);
			var_xg[n][g].set(GRB_StringAttr_VarName, buf);
		}
		for (int l = 0; l < numLine; l++) {
			std::sprintf(buf, "var_xl(%d,%d)", n, l);
			var_xl[n][l].set(GRB_StringAttr_VarName, buf);
		}
		for (int k = 0; k < numStg; k++) {
			std::sprintf(buf, "var_xs(%d,%d)", n, k);
			var_xs[n][k].set(GRB_StringAttr_VarName, buf);
	}
#endif
}

	/* *********************** Lower level operational decisions *********************** */
	const int numLowerIndices = numNodes * Season::numSeasons * ucPeriods;
	varGenStatus = new GRBVar *[numLowerIndices];											//binary variable for on/off status of plants alpha[n][d][k][g]
	varStartUp = new GRBVar *[numLowerIndices];												//binary variable for start-up/shut-down of plants alpha[n][d][k][g]
	varPowerGen = new GRBVar *[numLowerIndices];												//amount of power generated from generator p[n][d][k][g]
	varPowerFlow = new GRBVar **[numLowerIndices];												//power flow between bus i and bus j f[n][d][k][l]
	varPowerWithdrawal = new GRBVar *[numLowerIndices];										//total power withdrawn from storage facilities of bus u[n][d][k][sg]
	varPowerInject = new GRBVar *[numLowerIndices];											//total power injected into storage facilities of bus v[n][d][k][sg]
	varPowerRemaining = new GRBVar *[numLowerIndices];										//total remaining power in storage facilities of bus at the beginning of k v[n][d][k][sg]

	for (int n = 0; n < numNodes; n++) {
		for (int season = 0; season != Season::numSeasons; season++) {
			for (int h = 0; h < ucPeriods; h++) {
				const int index = getLowerIndexExpectedForm(n, season, h);
				varGenStatus[index] = MPmodel.addVars(numGen, GRB_BINARY);
				varStartUp[index] = MPmodel.addVars(numGen, GRB_BINARY);
				varPowerGen[index] = MPmodel.addVars(numGen);
				varPowerFlow[index] = new GRBVar *[numLine];
				varPowerWithdrawal[index] = MPmodel.addVars(numStg);
				varPowerInject[index] = MPmodel.addVars(numStg);
				varPowerRemaining[index] = MPmodel.addVars(numStg);
				for (int l = 0; l < numLine; l++) {
					varPowerFlow[index][l] = MPmodel.addVars(Direction::numDirections);
				}
#ifdef LPMODEL
				for (int g = 0; g < numGen; g++) {
					std::sprintf(buf, "varGenStatus(%d,%d,%d,%d)", n, season, h, g);
					varGenStatus[index][g].set(GRB_StringAttr_VarName, buf);
					std::sprintf(buf, "varStartUp(%d,%d,%d,%d)", n, season, h, g);
					varStartUp[index][g].set(GRB_StringAttr_VarName, buf);
					std::sprintf(buf, "varPowerGen(%d,%d,%d,%d)", n, season, h, g);
					varPowerGen[index][g].set(GRB_StringAttr_VarName, buf);
				}
				for (int l = 0; l < numLine; l++) {
					std::sprintf(buf, "varPowerFlow(%d,%d,%d,%d,+)", n, season, h, l);
					varPowerFlow[index][l][Direction::Forward].set(GRB_StringAttr_VarName, buf);
					std::sprintf(buf, "varPowerFlow(%d,%d,%d,%d,-)", n, season, h, l);
					varPowerFlow[index][l][Direction::Backward].set(GRB_StringAttr_VarName, buf);
				}
				for (int k = 0; k < numStg; k++) {
					std::sprintf(buf, "varPowerWithdrawal(%d,%d,%d,%d)", n, season, h, k);
					varPowerWithdrawal[index][k].set(GRB_StringAttr_VarName, buf);
					std::sprintf(buf, "varPowerInject(%d,%d,%d,%d)", n, season, h, k);
					varPowerInject[index][k].set(GRB_StringAttr_VarName, buf);
					std::sprintf(buf, "varPowerRemaining(%d,%d,%d,%d)", n, season, h, k);
					varPowerRemaining[index][k].set(GRB_StringAttr_VarName, buf);
			}
#endif			
		}
	}
	}

	/* *********************** Objective *********************** */
	GRBLinExpr xpr = 0;
	for (int n = 0; n < numNodes; n++) {
		const auto &node = scenarioNodes[n];
		for (int g = 0; g < numGen; g++) {
			const auto &gen = generators[g];
			xpr += node.probability * gen.cost[n] * var_xg[n][g];
		}
		for (int l = 0; l < numLine; l++) {
			const auto &line = lines[l];
			xpr += node.probability * line.cost[n] * var_xl[n][l];
		}
		for (int k = 0; k < numStg; k++) {
			const auto &stg = storages[k];
			xpr += node.probability * stg.cost[n] * var_xs[n][k];
		}
		for (int season = 0; season != Season::numSeasons; season++) {
			for (int h = 0; h < ucPeriods; h++) {
				const int lower_index = getLowerIndexExpectedForm(n, season, h);
				for (int g = 0; g < numGen; g++) {
					const auto &gen = generators[g];
					xpr += node.probability * get_opn_cost_factor(season) * (gen.startUpCost * varStartUp[lower_index][g] + gen.fuelCostB * gen.fuelPriceRate[n] * varPowerGen[lower_index][g]
						+ gen.fuelCostA * varGenStatus[lower_index][g]);
				}
			}
		}
	}
	MPmodel.setObjective(xpr, GRB_MINIMIZE);

	/* *********************** Constraints *********************** */
	// upper level constraints
	for (int n = 0; n < numNodes; n++) {
		const auto &node = scenarioNodes[n];
		if (node.period != paramReg->numPeriods - 1) continue;
		// upper level constraints
		for (int g = 0; g < numGen; g++) {
			if (generators[g].exists()) continue;
			GRBLinExpr tempExpr = 0;
			for (const int i : node.nodesInPath) {
				tempExpr += var_xg[i][g];
			}
			MPmodel.addConstr(tempExpr, GRB_LESS_EQUAL, 1 - static_cast<int>(generators[g].exists()), "UpperLevel");
		}
		for (int l = 0; l < numLine; l++) {
			GRBLinExpr tempExpr = 0;
			for (const int i : node.nodesInPath) {
				tempExpr += var_xl[i][l];
			}
			MPmodel.addConstr(tempExpr, GRB_LESS_EQUAL, 1 - static_cast<int>(lines[l].exists()), "UpperLevel");
		}
		for (int k = 0; k < numStg; k++) {
			GRBLinExpr tempExpr = 0;
			for (const int i : node.nodesInPath) {
				tempExpr += var_xs[i][k];
			}
			MPmodel.addConstr(tempExpr, GRB_LESS_EQUAL, 1 - static_cast<int>(storages[k].exists()), "UpperLevel");
		}
	}

	// lower level constraints
	for (int n = 0; n < numNodes; n++) {
		const auto &node = scenarioNodes[n];
		for (int season = 0; season != Season::numSeasons; season++) {
			for (int h = 0; h < ucPeriods; h++) {
				const int lower_index = getLowerIndexExpectedForm(n, season, h);
				for (int g = 0; g < numGen; g++) {
					const auto &gen = generators[g];
					GRBLinExpr tempExpr = int(gen.exists());
					for (const int i : node.nodesInPath) {
						tempExpr += var_xg[i][g];
					}

					std::sprintf(buf, "Gen(%d,%d,%d,%d)", n, season, h, g);																							// Coupling constraint			
					MPmodel.addConstr(varGenStatus[lower_index][g], GRB_LESS_EQUAL, tempExpr, buf);

					if (gen.is_renewable) {
						std::sprintf(buf, "RnGen(%d,%d,%d,%d)", n, season, h, g);
						MPmodel.addConstr(varPowerGen[lower_index][g], GRB_LESS_EQUAL, gen.supply_mean[season][h%HOURSPERDAY], buf);
					}

					std::sprintf(buf, "PowUB(%d,%d,%d,%d)", n, season, h, g);
					MPmodel.addConstr(varPowerGen[lower_index][g] - gen.capacity * varGenStatus[lower_index][g], GRB_LESS_EQUAL, 0, buf);

					std::sprintf(buf, "GenAct(%d,%d,%d,%d)", n, season, h, g);
					if (h == 0) {																														// weekly UC cconstraints
						MPmodel.addConstr(varStartUp[lower_index][g] - varGenStatus[lower_index][g], GRB_GREATER_EQUAL, 0, buf);
					}
					else
					{
						MPmodel.addConstr(varStartUp[lower_index][g] - varGenStatus[lower_index][g] + varGenStatus[getLowerIndexExpectedForm(n, season, h - 1)][g], GRB_GREATER_EQUAL, 0, buf);
					}

					if (h) {
						// min up and down
						for (int tau = h; tau < ucPeriods && tau < h + generators[g].minUpTime; tau++) {
							std::sprintf(buf, "MinUp(%d,%d,%d,%d,%d)", n, season, h, tau, g);
							MPmodel.addConstr(varGenStatus[lower_index][g] - varGenStatus[getLowerIndexExpectedForm(n, season, h - 1)][g], GRB_LESS_EQUAL,
								varGenStatus[getLowerIndexExpectedForm(n, season, tau)][g], buf);
						}
						for (int tau = h; tau < ucPeriods && tau < h + generators[g].minDownTime; tau++) {
							std::sprintf(buf, "MinDown(%d,%d,%d,%d,%d)", n, season, h, tau, g);
							MPmodel.addConstr(varGenStatus[lower_index][g] - varGenStatus[getLowerIndexExpectedForm(n, season, h - 1)][g], GRB_GREATER_EQUAL,
								varGenStatus[getLowerIndexExpectedForm(n, season, tau)][g] - 1, buf);
						}
						// ramp up and down
						std::sprintf(buf, "RampUp(%d,%d,%d,%d)", n, season, h, g);
						MPmodel.addConstr(varPowerGen[lower_index][g] - varPowerGen[getLowerIndexExpectedForm(n, season, h - 1)][g], GRB_LESS_EQUAL, generators[g].rampUp);
						std::sprintf(buf, "RampDown(%d,%d,%d,%d)", n, season, h, g);
						MPmodel.addConstr(varPowerGen[getLowerIndexExpectedForm(n, season, h - 1)][g] - varPowerGen[lower_index][g], GRB_LESS_EQUAL, generators[g].rampDown);
					}
				}
				for (int l = 0; l < numLine; l++) {
					const auto &line = lines[l];
					GRBLinExpr tempExpr = 0;
					if (!line.exists()) {
						for (const int i : node.nodesInPath) {
							tempExpr += var_xl[i][l];
						}
					}
					tempExpr *= line.capacity;
					std::sprintf(buf, "Lne(%d,%d,%d,%d)", n, season, h, l);
					MPmodel.addConstr(varPowerFlow[lower_index][l][Direction::Forward] + varPowerFlow[lower_index][l][Direction::Backward] - tempExpr,
						GRB_LESS_EQUAL, line.currentCapacity, buf);
				}
				for (int k = 0; k < numStg; k++) {
					const auto &storage = storages[k];
					GRBLinExpr tempExpr = 0;
					if (!storage.exists()) {
						for (const int i : node.nodesInPath) {
							tempExpr += var_xs[i][k];
						}
					}
					tempExpr *= storage.capacity;
					std::sprintf(buf, "Stg(%d,%d,%d,%d)", n, season, h, k);																							//Constraint 23
					MPmodel.addConstr(varPowerRemaining[lower_index][k] - tempExpr, GRB_LESS_EQUAL, storage.currentCapacity, buf);

					std::sprintf(buf, "StgLmt(%d,%d,%d,%d)", n, season, h, k);
					MPmodel.addConstr(varPowerWithdrawal[lower_index][k] - varPowerRemaining[lower_index][k], GRB_LESS_EQUAL, 0, buf);														//(21)
					std::sprintf(buf, "StgBal(%d,%d,%d,%d)", n, season, h, k);
					if (h == 0) {																														//in the 1st lower level period, reset the storage facility
						MPmodel.addConstr(varPowerRemaining[lower_index][k], GRB_EQUAL, 0, buf);
					}
					else {
						MPmodel.addConstr(varPowerRemaining[lower_index][k] - varPowerRemaining[getLowerIndexExpectedForm(n, season, h - 1)][k]
							- varPowerInject[getLowerIndexExpectedForm(n, season, h - 1)][k] + varPowerWithdrawal[getLowerIndexExpectedForm(n, season, h - 1)][k], GRB_EQUAL, 0, buf);							//Constraint 22
					}
				}
				for (int b = 0; b < numBus; b++) {
					std::sprintf(buf, "Dmd(%d,%d,%d,%d)", n, season, h, b);
					GRBLinExpr tempExpr = 0;
					for (int l = 0; l < numLine; l++) {
					   if (lines[l].from == b)	tempExpr -= varPowerFlow[lower_index][l][Direction::Forward] - varPowerFlow[lower_index][l][Direction::Backward] * (1 - lines[l].lossRatio);
					   if (lines[l].to == b) tempExpr += varPowerFlow[lower_index][l][Direction::Forward] * (1 - lines[l].lossRatio) - varPowerFlow[lower_index][l][Direction::Backward];
					}
					for (int g = 0; g < numGen; g++) {
						if (generators[g].loc == b) tempExpr += varPowerGen[lower_index][g];
					}
					for (int k = 0; k < numStg; k++) {
						if (storages[k].loc == b) tempExpr += varPowerWithdrawal[lower_index][k] * paramReg->stgCeff - varPowerInject[lower_index][k];
					}
					MPmodel.addConstr(tempExpr, GRB_GREATER_EQUAL, buses[b].season_hourly_loads[n][season][h], buf);
				}
			}			
		}
	}

#ifdef LPMODEL
	sprintf(buf, "%s/expected_form.lp", output_directory.c_str());
	MPmodel.write(buf);
#endif

	clock_t _start = clock();
	auto _result = solveMasterProblem();
	clock_t _end = clock();
	cmp_time = double((_end - _start)) / double(CLOCKS_PER_SEC);

	/*if (_result != Result::INFEASIBLE)
	{
		computeStochasticValue();
	}*/

	return _result;
}

void Model::computeStochasticValue()
{
	if (!is_master_proc())	return;

	PRINT_SUBSECTION("Compute Operational Costs");

	GRBModel grbModel = GRBModel(*env);
	grbModel.set(GRB_IntParam_OutputFlag, 0);
	// Decision Variables
	GRBVar** varGenStatus = new GRBVar *[ucPeriods];
	GRBVar** varStartUp = new GRBVar *[ucPeriods];
	GRBVar** varPowerGen = new GRBVar *[ucPeriods];
	GRBVar*** varPowerFlow = new GRBVar **[ucPeriods];
	GRBVar** varPowerWithdrawal = new GRBVar *[ucPeriods];
	GRBVar** varPowerInject = new GRBVar *[ucPeriods];
	GRBVar** varPowerRemaining = new GRBVar *[ucPeriods];

	for (int h = 0; h < ucPeriods; h++) {
		varGenStatus[h] = grbModel.addVars(numGen, GRB_BINARY);
		varStartUp[h] = grbModel.addVars(numGen, GRB_BINARY);
		varPowerGen[h] = grbModel.addVars(numGen);
		varPowerFlow[h] = new GRBVar *[numLine];
		varPowerWithdrawal[h] = grbModel.addVars(numStg);
		varPowerInject[h] = grbModel.addVars(numStg);
		varPowerRemaining[h] = grbModel.addVars(numStg);

		for (int l = 0; l < numLine; l++) {
			varPowerFlow[h][l] = grbModel.addVars(Direction::numDirections);
		}
	}

	GRBVar* unmetDemand = grbModel.addVars(ucPeriods * numBus);

	/* ************************* Objective **************************** */
	GRBLinExpr expr = 0;
	for (int h = 0; h < ucPeriods; h++) {
		for (int g = 0; g < numGen; g++) {
			const auto &gen = generators[g];
			expr += gen.startUpCost * varStartUp[h][g] + gen.fuelCostB * varPowerGen[h][g] + gen.fuelCostA * varGenStatus[h][g];
		}
		for (int b = 0; b < numBus; b++) {
			expr += paramReg->bigM * unmetDemand[h * numBus + b];
		}
	}
	grbModel.setObjective(expr, GRB_MINIMIZE);

	/* ************************* Constraints **************************** */
	GRBConstr* GenExistConstrs = new GRBConstr[ucPeriods * numGen];
	GRBConstr* GenLimitConstrs = new GRBConstr[ucPeriods * numGen];
	GRBConstr* RnGenConstrs = new GRBConstr[ucPeriods * numRnGen];
	GRBConstr* LineLimitConstrs = new GRBConstr[ucPeriods * numLine];
	GRBConstr* StgLimitConstrs = new GRBConstr[ucPeriods * numStg];
	GRBConstr* DemandConstrs = new GRBConstr[ucPeriods * numBus];
	for (int h = 0; h < ucPeriods; h++) {
		for (int g = 0; g < numGen; g++) {
			const auto &gen = generators[g];
			GenExistConstrs[h * numGen + g] =
				grbModel.addConstr(varGenStatus[h][g], GRB_LESS_EQUAL, 0);				//RHS will be updated later in the loop
			GenLimitConstrs[h * numGen + g] =
				grbModel.addConstr(varPowerGen[h][g], GRB_LESS_EQUAL, gen.capacity);
			
			grbModel.addConstr(varPowerGen[h][g] - gen.capacity * varGenStatus[h][g], GRB_LESS_EQUAL, 0);

			if (h == 0) {
				grbModel.addConstr(varStartUp[h][g] - varGenStatus[h][g], GRB_GREATER_EQUAL, 0);
			}
			else
			{
				grbModel.addConstr(varStartUp[h][g] - varGenStatus[h][g] + varGenStatus[h - 1][g], GRB_GREATER_EQUAL, 0);
			}		
		}
		for (int l = 0; l < numLine; l++) {
			const auto &line = lines[l];
			LineLimitConstrs[h * numLine + l] =
				grbModel.addConstr(varPowerFlow[h][l][Direction::Forward] + varPowerFlow[h][l][Direction::Backward], GRB_LESS_EQUAL, line.currentCapacity);						//RHS will be updated in benders algorithm											
			// grbModel.addConstr(varPowerFlow[index][l], GRB_LESS_EQUAL, line.capacity);
		}
		for (int k = 0; k < numStg; k++) {
			const auto &storage = storages[k];
			StgLimitConstrs[h * numStg + k] =
				grbModel.addConstr(varPowerRemaining[h][k], GRB_LESS_EQUAL, storage.currentCapacity);		//RHS will be updated in benders algorithm

			grbModel.addConstr(varPowerWithdrawal[h][k] - varPowerRemaining[h][k], GRB_LESS_EQUAL, 0);
			if (h == 0) {
				grbModel.addConstr(varPowerRemaining[h][k], GRB_EQUAL, 0);
			}
			else {
				grbModel.addConstr(varPowerRemaining[h][k] - varPowerRemaining[h - 1][k]
					- varPowerInject[h - 1][k] + varPowerWithdrawal[h - 1][k], GRB_EQUAL, 0);
			}
		}
		for (int b = 0; b < numBus; b++) {
			GRBLinExpr tempExpr = 0;
			for (int l = 0; l < numLine; l++) {
			   if (lines[l].from == b)	tempExpr -= varPowerFlow[h][l][Direction::Forward] - varPowerFlow[h][l][Direction::Backward] * (1 - lines[l].lossRatio);
			   if (lines[l].to == b) tempExpr += varPowerFlow[h][l][Direction::Forward] * (1 - lines[l].lossRatio) - varPowerFlow[h][l][Direction::Backward];
			}
			for (int g = 0; g < numGen; g++) {
				if (generators[g].loc == b) tempExpr += varPowerGen[h][g];
			}
			for (int k = 0; k < numStg; k++) {
				if (storages[k].loc == b) tempExpr += varPowerWithdrawal[h][k] * paramReg->stgCeff - varPowerInject[h][k];
			}
			DemandConstrs[h * numBus + b] =
				grbModel.addConstr(tempExpr + unmetDemand[h * numBus + b], GRB_GREATER_EQUAL, 0);
		}
	}
	
	// for every scenario node and every lower-level scenario, solve the corresponding UC problem and get operational costs and unmet demands	
	for (int n = 0; n < numNodes; n++) {
		operationalCost[n] = 0;
		for (int s = 0; s != Season::numSeasons; s++) {
		   for (int d = 0; d < paramReg->numDmdScenarios; d++) {
				// update and solve subproblem here
				for (int h = 0; h < ucPeriods; h++) {
					for (int g = 0; g < numGen; g++) {
						const auto &gen = generators[g];
						int _exist = gen.exists();
						for (const int i : scenarioNodes[n].nodesInPath) {
							_exist += int(sol_xg[i][g]);
						}
						GenExistConstrs[h * numGen + g].set(GRB_DoubleAttr_RHS, _exist);

						GenLimitConstrs[h * numGen + g].set(GRB_DoubleAttr_RHS, gen.is_renewable ? gen.supply_scenarios[get2ndStageIndex(n, s, d, h)] : gen.capacity);
					}
					for (int l = 0; l < numLine; l++) {
						const auto &line = lines[l];
						double cap = line.currentCapacity;
						for (const int i : scenarioNodes[n].nodesInPath) {
							cap += sol_xl[i][l] * line.capacity;
						}
						LineLimitConstrs[h * numLine + l].set(GRB_DoubleAttr_RHS, cap);
					}
					for (int k = 0; k < numStg; k++) {
						const auto &storage = storages[k];
						double cap = storage.currentCapacity;
						for (const int i : scenarioNodes[n].nodesInPath) {
							cap += sol_xs[i][k] * storage.capacity;
						}
						StgLimitConstrs[h * numStg + k].set(GRB_DoubleAttr_RHS, cap);
					}
					for (int b = 0; b < numBus; b++) {
						DemandConstrs[h * numBus + b].set(GRB_DoubleAttr_RHS, buses[b].season_hourly_load_scenarios[get2ndStageIndex(n, s, d, h)]);
					}
				}

				grbModel.optimize();

				operationalCost[n] += scenarioNodes[n].probability * paramReg->prob_dmd_scenario * grbModel.get(GRB_DoubleAttr_ObjVal);
			}
		}
	}
}


