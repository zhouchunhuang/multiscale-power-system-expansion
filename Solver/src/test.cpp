#include "model.h"

// This is used for testing some algorithm in a single machine
void Model::test()
{
	PRINT_SECTION("Testing");
	std::cout.precision(4);

	rank = 0; master = 0; 
	color = 0; sub_master = 0; sub_rank = 0;
	double mpRunTime = 0;

	initialize(true);

	clock_t _start = clock(), _end = 0;

	createMasterProblem();

	createSecondaryMP();
	SMPmodel.set(GRB_IntParam_OutputFlag, 0);

	createSubproblem(true);
	SPmodel.set(GRB_IntParam_OutputFlag, 0);

	_start = clock();
	/* *********************************************************** Nested Decomposition Algorithm********************************************************** */
	MPmodel.optimize();																					//solve master problem
	mpRunTime = MPmodel.get(GRB_DoubleAttr_Runtime);
	MpUB = MPmodel.get(GRB_DoubleAttr_ObjVal);
	getMpDuals();

	updateSecondaryMP();					// update secondary MP objective; constraints will not be changed; benders cuts are valid in future upper-level iterations. 

	PRINT_SECTION("Phase I");
	std::cout << "No.\t\tNUB\t\tUB\t\tLB\t\tGap\t\tTime" << std::endl;
	sub_iteration = 0;
	while (true) {
		sub_iteration++;

		solveSecondaryMP();

		testSolveSPsPhaseI();

		getSecondaryMpBounds();

		time_t _end = clock();
		const double smp_gap = fabs((SmpUB - SmpLB) / SmpUB);
		const double test_time = double((double(_end - _start)) / CLOCKS_PER_SEC);

		std::cout << sub_iteration << "\t\t" << SmpNewUB << "\t\t" << SmpUB << "\t\t" << SmpLB << "\t\t" << smp_gap << "\t\t" << test_time << std::endl;
		
		if (smp_gap <= paramReg->smp_gap) break;
	}

	PRINT_SECTION("Phase II");
	std::cout << "No.\t\tNUB\t\tUB\t\tLB\t\tGap\t\tTime" << std::endl;
	sub_iteration = 0;
	SmpUB = +paramReg->bigM; SmpLB = -paramReg->bigM; SmpNewUB = paramReg->bigM;
	sub_message = 0;
	sub_iteration = 0;
	while (true) {
		sub_iteration++;

		solveSecondaryMP();

		testSolveSPsPhaseII();

		getSecondaryMpBounds();

		time_t _end = clock();
		const double smp_gap = fabs((SmpUB - SmpLB) / SmpUB);
		const double test_time = double((double(_end - _start)) / CLOCKS_PER_SEC);

		std::cout << sub_iteration << "\t\t" << SmpNewUB << "\t\t" << SmpUB << "\t\t" << SmpLB << "\t\t" << smp_gap << "\t\t" << test_time << std::endl;

		if (smp_gap <= paramReg->smp_gap) break;
	}

	system("pause");

	exit(EXIT_SUCCESS);
}

void Model::testSolveSPsPhaseI()
{
	for (int i = 0; i < numSPsPerNode; i++) {
		const int season = get_season_of_subproblem(i);
		// update subproblem
		for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
			for (int h = 0; h < ucPeriods; h++) {
				for (int g = 0; g < numGen; g++) {
					GenLimitConstrs[d * ucPeriods * numGen + h * numGen + g].set(GRB_DoubleAttr_RHS, CapGen[g]);

					SPmodel.chgCoeff(GenLimit2Constrs[d * ucPeriods * numGen + h * numGen + g], varGenStatus[getLowerIndexInSP(d, h)][g],
						generators[g].is_renewable ? (-generators[g].supply_scenarios[getLowerIndex(color, season, get_dmd_scenario(d, i), h)]) : (-generators[g].capacity));

					varGenStatus[getLowerIndexInSP(d, h)][g].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
					varGenStatus[getLowerIndexInSP(d, h)][g].set(GRB_DoubleAttr_UB, GRB_INFINITY);
				}
				for (int l = 0; l < numLine; l++) {
					LineLimitConstrs[d * ucPeriods * numLine + h * numLine + l].set(GRB_DoubleAttr_RHS, CapLine[l]);
				}
				for (int k = 0; k < numStg; k++) {
					StgLimitConstrs[d * ucPeriods * numStg + h * numStg + k].set(GRB_DoubleAttr_RHS, CapStg[k]);
				}
				for (int b = 0; b < numBus; b++) {
					DemandConstrs[d * ucPeriods * numBus + h * numBus + b].set(GRB_DoubleAttr_RHS, buses[b].season_hourly_load_scenarios[
							get2ndStageIndex(color, season, get_dmd_scenario(d, i), h)]);
				}
			}
		}
		// objective
		GRBLinExpr spobj_exp = 0;
		int index = 0;
		for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
			for (int h = 0; h < ucPeriods; h++, index++) {
				for (int g = 0; g < numGen; g++) {
					const auto &gen = generators[g];
					spobj_exp += scenarioNodes[color].probability * paramReg->prob_dmd_scenario * (gen.startUpCost * varStartUp[index][g] + gen.fuelCostA * varGenStatus[index][g]) +
					      scenarioNodes[color].probability * paramReg->prob_dmd_scenario * gen.fuelCostB * varPowerGen[index][g];
				}
			}
		}
		SPmodel.setObjective(get_opn_cost_factor(season) * spobj_exp, GRB_MINIMIZE);

		try {			
			SPmodel.optimize();
			
			//get Duals or Extreme rays
			const bool infeasible = (SPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE);
			GRBLinExpr cutLhs = 0;
			GRBLinExpr cutRhs = 0;

			if (!infeasible) {
				cutLhs = var_theta[i];
				SPofv[i] = SPmodel.get(GRB_DoubleAttr_ObjVal);
				SPstatus[i] = SpStatus::LRoptimal;			
			}
			else {
				SPstatus[i] = SpStatus::LRinfeasible;
			}

			double dual = 0;	

			for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
				for (int h = 0; h < ucPeriods; h++) {
					for (int g = 0; g < numGen; g++) {
						const auto &gen = generators[g];

						dual = infeasible ? -abs(GenLimitConstrs[d * ucPeriods * numGen + h * numGen + g].get(GRB_DoubleAttr_FarkasDual)) :
							GenLimitConstrs[d * ucPeriods * numGen + h * numGen + g].get(GRB_DoubleAttr_Pi);
						cutRhs += dual * gen.currentCapacity;
						cutRhs += dual * var_zg[g] * gen.capacity;

						dual = dual = infeasible ? -abs(RampUpConstrs[d * ucPeriods * numGen + h * numGen + g].get(GRB_DoubleAttr_FarkasDual)) : 
							RampUpConstrs[d * ucPeriods * numGen + h * numGen + g].get(GRB_DoubleAttr_Pi);
						cutRhs += dual * ((h == 0) ? gen.capacity : gen.rampUp);

						dual = infeasible ? -abs(RampDownConstrs[d * ucPeriods * numGen + h * numGen + g].get(GRB_DoubleAttr_FarkasDual)) :
							RampDownConstrs[d * ucPeriods * numGen + h * numGen + g].get(GRB_DoubleAttr_Pi);
						cutRhs += dual * ((h == 0) ? gen.capacity : gen.rampDown);
					}
					for (int l = 0; l < numLine; l++) {
						dual = infeasible ? -abs(LineLimitConstrs[d * ucPeriods * numLine + h * numLine + l].get(GRB_DoubleAttr_FarkasDual)) : 
							LineLimitConstrs[d * ucPeriods * numLine + h * numLine + l].get(GRB_DoubleAttr_Pi);
						cutRhs += dual * lines[l].currentCapacity;
						cutRhs += dual * var_zl[l] * lines[l].capacity;
					}
					for (int k = 0; k < numStg; k++) {
						dual = infeasible ? -abs(StgLimitConstrs[d * ucPeriods * numStg + h * numStg + k].get(GRB_DoubleAttr_FarkasDual)) : 
							StgLimitConstrs[d * ucPeriods * numStg + h * numStg + k].get(GRB_DoubleAttr_Pi);
						cutRhs += dual * storages[k].currentCapacity;
						cutRhs += dual * var_zs[k] * storages[k].capacity;
					}
					for (int b = 0; b < numBus; b++) {
						dual = infeasible ? abs(DemandConstrs[d * ucPeriods * numBus + h * numBus + b].get(GRB_DoubleAttr_FarkasDual)) : 
							DemandConstrs[d * ucPeriods * numBus + h * numBus + b].get(GRB_DoubleAttr_Pi);
						cutRhs += dual * buses[b].season_hourly_load_scenarios[get2ndStageIndex(color, get_season_of_subproblem(i), get_dmd_scenario(d, i), h)];
					}
				}
			}

			for (const auto &it : MinDownConstrs) {
				dual = infeasible ? -abs(it.get(GRB_DoubleAttr_FarkasDual)) : it.get(GRB_DoubleAttr_Pi);
				cutRhs += dual;	// RHS of min down constraints are always 1
			}
			SMPmodel.addConstr(cutLhs, GRB_GREATER_EQUAL, cutRhs, "BendersCut");
		}
		catch (GRBException &e) {
			std::cerr << "Error while solving subproblem " << i << "! error code = " << e.getErrorCode() << "; error message: " << e.getMessage() << std::endl;
			std::cout << "The obj of the subproblem = " << SPofv[i] << std::endl;
			std::cout << "The status of the subproblem = " << SPmodel.get(GRB_IntAttr_Status) << std::endl;
			char buf[100];
			std::sprintf(buf, "%s/sp-%d-%d.lp", output_directory.c_str(), season, i);
			SPmodel.write(buf);
			exit(EXIT_FAILURE);
		}
	}
}

void Model::testSolveSPsPhaseII()
{
	for (int i = 0; i < numSPsPerNode; i++) {
		recFuncLB[i] = SPofv[i];
	}

	for (int i = 0; i < numSPsPerNode; i++) {
		const int season = get_season_of_subproblem(i);
		// update subproblem
		for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
			for (int h = 0; h < ucPeriods; h++) {
				for (int g = 0; g < numGen; g++) {
					GenLimitConstrs[d * ucPeriods * numGen + h * numGen + g].set(GRB_DoubleAttr_RHS, CapGen[g]);

					SPmodel.chgCoeff(GenLimit2Constrs[d * ucPeriods * numGen + h * numGen + g], varGenStatus[getLowerIndexInSP(d, h)][g],
						generators[g].is_renewable ? (-generators[g].supply_scenarios[getLowerIndex(color, season, get_dmd_scenario(d, i), h)]) : (-generators[g].capacity));

					varGenStatus[getLowerIndexInSP(d, h)][g].set(GRB_CharAttr_VType, GRB_BINARY);
				}
				for (int l = 0; l < numLine; l++) {
					LineLimitConstrs[d * ucPeriods * numLine + h * numLine + l].set(GRB_DoubleAttr_RHS, CapLine[l]);
				}
				for (int k = 0; k < numStg; k++) {
					StgLimitConstrs[d * ucPeriods * numStg + h * numStg + k].set(GRB_DoubleAttr_RHS, CapStg[k]);
				}
				for (int b = 0; b < numBus; b++) {
					DemandConstrs[d * ucPeriods * numBus + h * numBus + b].set(GRB_DoubleAttr_RHS, buses[b].season_hourly_load_scenarios[
						get2ndStageIndex(color, season, get_dmd_scenario(d, i), h)]);
				}
			}
		}
		// objective
		GRBLinExpr spobj_exp = 0;
		int index = 0;
		for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
			for (int h = 0; h < ucPeriods; h++, index++) {
				for (int g = 0; g < numGen; g++) {
					const auto &gen = generators[g];
					spobj_exp += scenarioNodes[color].probability * paramReg->prob_dmd_scenario * (gen.startUpCost * varStartUp[index][g] + gen.fuelCostA * varGenStatus[index][g]) +
					      scenarioNodes[color].probability * paramReg->prob_dmd_scenario * gen.fuelCostB * varPowerGen[index][g];
				}
			}
		}
		SPmodel.setObjective(get_opn_cost_factor(season) * spobj_exp, GRB_MINIMIZE);

		try {
			// solve the subproblem 
			SPmodel.optimize();

			if (SPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
				std::cerr << "SP " << i << "is infeasible" << std::endl;
				exit(EXIT_FAILURE);
			}

			SPofv[i] = SPmodel.get(GRB_DoubleAttr_ObjVal);

			const double recourseValue = SPmodel.get(GRB_DoubleAttr_ObjVal);
			zEqOne = 0;
			combnXpr = 0;
            for (int g = 0; g < numGen; g++) {
                if (var_zg[g].get(GRB_DoubleAttr_X) > 0.5) {
                    zEqOne++;
                    combnXpr += var_zg[g];
                }
                else {
                    combnXpr -= var_zg[g];
                }
            }
            for (int l = 0; l < numLine; l++) {
                if (var_zl[l].get(GRB_DoubleAttr_X) > 0.5) {
                    zEqOne++;
                    combnXpr += var_zl[l];
                }
                else {
                    combnXpr -= var_zl[l];
                }
            }
            for (int k = 0; k < numStg; k++) {
                if (var_zs[k].get(GRB_DoubleAttr_X) > 0.5) {
                    zEqOne++;
                    combnXpr += var_zs[k];
                }
                else {
                    combnXpr -= var_zs[k];
                }
            }

			GRBLinExpr cut2Rhs = (recourseValue - recFuncLB[i]) * (combnXpr - zEqOne + 1) + recFuncLB[i];
			SMPmodel.addConstr(var_theta[i], GRB_GREATER_EQUAL, cut2Rhs, "IntegerLshapedOptCut");
		}
		catch (GRBException &e) {
			std::cerr << "Error while solving subproblem " << i << "! error code = " << e.getErrorCode() << "; error message: " << e.getMessage() << std::endl;
			std::cout << "The obj of the subproblem = " << SPofv[i] << std::endl;
			std::cout << "The status of the subproblem = " << SPmodel.get(GRB_IntAttr_Status) << std::endl;
			char buf[100];
			std::sprintf(buf, "%s/sp-%d-%d.lp", output_directory.c_str(), season, i);
			SPmodel.write(buf);
			exit(EXIT_FAILURE);
		}
	}
}