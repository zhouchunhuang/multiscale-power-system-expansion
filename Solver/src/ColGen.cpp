#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

#include "model.h"


using namespace std;

Result Model::columnGeneration()
{
	if (is_master_proc()) PRINT_SECTION("Solving the Problem with Column Generation")
	ReadWrite rw;
	rw.set_outputDirectory(output_directory);
	double mpRunTime = 0;

	initialize();
	
	auto _start = std::chrono::steady_clock::now();

	createMasterProblem();

	createUpperSubproblem();

	if(is_master_proc())
	{																
		rw.printBounds(0, MpUB, nullptr, numNodes, MpLB, MpUB);
#ifdef StreamAlgInfo
		std::cout << "UB\t\t";
		for (int n = 0; n < numNodes; n++) {
			std::cout << n << "\t\t";
		}
		std::cout << "LB" << std::endl;
#endif
	}																		
		
	///==================================Std Column Generation Algorithm======================================///
	iteration = 0;
	while(true){
		iteration++;
		//................Master processor..................//
		if(is_master_proc())																
		{
            MPmodel.optimize();
            if (MPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
#ifdef LPMODEL
                char buf[100];
				std::sprintf(buf, "%s/mp_infeasible.lp", output_directory.c_str());
				MPmodel.write(buf);
#endif // LPMODEL
                std::cerr << "Primary Master Problem is infeasible!" << std::endl;
                exit(EXIT_FAILURE);
            }
            else if (MPmodel.get(GRB_IntAttr_Status) == GRB_NUMERIC) {
                MPmodel.set(GRB_IntParam_Crossover, 1);
                MPmodel.optimize();
                MPmodel.set(GRB_IntParam_Crossover, 0);
            }
			mpRunTime = MPmodel.get(GRB_DoubleAttr_Runtime);
			MpUB = MPmodel.get(GRB_DoubleAttr_ObjVal);
			getMpDuals();
		}	
		for (int n = 0; n < numNodes; n++) {
			const auto &node = scenarioNodes[n];
			MPI_Bcast(DualGen[n], numGen, MPI_DOUBLE, master, MPI_COMM_WORLD);
			MPI_Bcast(DualLine[n], numLine, MPI_DOUBLE, master, MPI_COMM_WORLD);
			MPI_Bcast(DualStg[n], numStg, MPI_DOUBLE, master, MPI_COMM_WORLD);
			MPI_Bcast(DualConvex, numNodes, MPI_DOUBLE, master, MPI_COMM_WORLD);
		}
		
		solveUpperSubproblem();													/* ******** Solve subproblem and Get a feasible expansion plan as new columns ******** */
		
		addColumns();															/* ******** Add columns to master problem ******** */

		getMPbounds();
		
		MPI_Bcast(&message, 1, MPI_INT, master, MPI_COMM_WORLD);

		if (is_master_proc()){
			if (gap <= paramReg->mp_gap) message = 1;															// generate columns until the gap is relatively small
			auto _end = std::chrono::steady_clock::now();
			cmp_time = double(std::chrono::duration_cast<std::chrono::seconds>(_end - _start).count());
			if (cmp_time >= paramReg->time_lmt)
				message = 1;

#ifdef StreamAlgInfo
			std::cout << mpRunTime << "\t\t";
			for (int n = 0; n < numNodes; n++) {
				std::cout << PspTime[n] << "(" << SmpGap[n] << ")" << "\t\t";
			}
			std::cout << std::endl;

			std::cout << MpUB << "\t\t";
			for (int n = 0; n < numNodes; n++) {
				std::cout << ReducedCost[n] << "\t\t";
			}
				
			std::cout << MpLB << std::endl;
#endif
			rw.printBounds(iteration, MpUB, ReducedCost, numNodes, MpLB, gap);
		}
		MPI_Bcast(&message,1,MPI_INT,master,MPI_COMM_WORLD);
		if(message) break;		
	}

	solveMasterProblem();
	
	MPI_Finalize();                         // terminate MPI
	return Result::OPTIMAL;
}

void Model::initialize(bool nested)
{
	// --------------Feasible Expansion Plan ------------------//
	ColExpPlanGen = new double*[numNodes];
	ColExpPlanLine = new double*[numNodes];
	ColExpPlanStg = new double*[numNodes];
	// --------------Dual Solution From Master Problem ------------------//
	DualGen = new double*[numNodes];
	DualLine = new double*[numNodes];
	DualStg = new double*[numNodes];
	DualConvex = new double[numNodes];
	// --------------Store Reduced Cost obtained by Solving Subproblems ------------------//
	ReducedCost = new double[numNodes];
	LambdaCoef = new double[numNodes];
	SmpOFV = new double[numNodes];
	SmpGap = new double[numNodes];
   PspTime = new double[numNodes];
   PspSmpTime = new double[numNodes];
   PspPhase2Time = new double[numNodes];

	for (int n = 0; n < numNodes; n++) {
		ColExpPlanGen[n] = new double[numGen];
		ColExpPlanLine[n] = new double[numLine];
		ColExpPlanStg[n] = new double[numStg];

		DualGen[n] = new double[numGen];
		DualLine[n] = new double[numLine];
		DualStg[n] = new double[numStg];

		ReducedCost[n] = paramReg->bigM;
		LambdaCoef[n] = paramReg->bigM;
		SmpOFV[n] = paramReg->bigM;
		SmpGap[n] = 1;
      PspTime[n] = 0;
      PspSmpTime[n] = 0;
      PspPhase2Time[n] = 0;
	}
	//-------------------------------------------------------//
    // new arrays to store columns
    col_gen_expansion = new double**[numNodes];
    col_line_expansion = new double**[numNodes];
    col_stg_expansion = new double**[numNodes];
    col_obj = new double*[numNodes];
    for(int n = 0; n < numNodes; n++) {
        const auto& node = scenarioNodes[n];
        const int num_solutions = NUM_COLUMNS;
        col_obj[n] = new double[num_solutions];
        col_gen_expansion[n] = new double *[num_solutions];
        col_line_expansion[n] = new double *[num_solutions];
        col_stg_expansion[n] = new double *[num_solutions];
        for (int i = 0; i < num_solutions; i++) {
            col_obj[n][i] = -1;
            col_gen_expansion[n][i] = new double[numGen];
            col_line_expansion[n][i] = new double[numLine];
            col_stg_expansion[n][i] = new double[numStg];
        }
    }
	/* ----------------------------------------------------- */
	message = 0;

	MpUB = paramReg->bigM;
	MpLB = -paramReg->bigM;
	gap = paramReg->bigM;
	MpNewLB = 0;

	//-------------------------------------------------------//
	if (nested) {
	   numDmdScenariosPerSubproblem = int(paramReg->numDmdScenarios / paramReg->numSPsPerSeason);
        numSPsPerNode = Season::numSeasons * paramReg->numSPsPerSeason;

		if (alg == Algorithm::TEST) {
			numCoresPerNode = 1;
			numSPsPerCore = numSPsPerNode;
		}
		else {
			if (num_ranks % numNodes) {
				std::cerr << "Incorrect parallel computing setting: numCoresPerNode not integer!" << std::endl;
				exit(EXIT_FAILURE);
			}
			numCoresPerNode = int(num_ranks / numNodes);

			if (numSPsPerNode % numCoresPerNode) {
				std::cerr << "Incorrect parallel computing setting: numSPsPerCore not integer!" << std::endl;
				exit(EXIT_FAILURE);
			}			
			numSPsPerCore = int(numSPsPerNode / numCoresPerNode);

			color = int(rank / numCoresPerNode);
			MPI_Comm_split(MPI_COMM_WORLD, color, rank, &sub_comm);			// split the world to a set of communicators
			MPI_Comm_rank(sub_comm, &sub_rank);								// get sub-world rank

			sub_master = 0;
			sub_message = 0;
		}
	
		ColExpPlanGenTemp = new double*[numNodes];
		ColExpPlanLineTemp = new double*[numNodes];
		ColExpPlanStgTemp = new double*[numNodes];
		for (int n = 0; n < numNodes; n++) {
			ColExpPlanGenTemp[n] = new double[numGen];
			ColExpPlanLineTemp[n] = new double[numLine];
			ColExpPlanStgTemp[n] = new double[numStg];
		}

		BendersInfCuts.clear();
	}
}

void Model::createMasterProblem()
{
	if (!is_master_proc()) return;

	PRINT_SUBSECTION("create Master Problem");
	char buf[100];

#ifdef SolverStreamOff
	MPmodel.set(GRB_IntParam_OutputFlag, 0);
#endif
	
	MPmodel.set(GRB_IntParam_Method, GRB_METHOD_BARRIER);
	MPmodel.set(GRB_IntParam_Crossover, 0);
	try {
        /* ********************** Decision variables ************************** */
        var_xg = new GRBVar *[numNodes];
        var_xl = new GRBVar *[numNodes];
        var_xs = new GRBVar *[numNodes];
        columns = std::vector<std::vector<GRBVar> >(numNodes);
        for (int n = 0; n < numNodes; n++) {
            var_xg[n] = MPmodel.addVars(numGen);
            var_xl[n] = MPmodel.addVars(numLine);
            var_xs[n] = MPmodel.addVars(numStg);
            columns[n].clear();
            columns[n].push_back(MPmodel.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS));

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

            std::sprintf(buf, "column(%d,0)", n);
            columns[n].back().set(GRB_StringAttr_VarName, buf);
            // var_lambda[n].set(GRB_DBL_ATTR_UB, 1);
#endif
        }

        /* ************************* Objective **************************** */
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
            xpr += node.probability * BIGM * columns[n].back();
        }
        MPmodel.setObjective(xpr, GRB_MINIMIZE);

        /* ********************** Constraints ************************** */
        SplitGen = new GRBConstr *[numNodes];                                            // multi-dimensional constraints for splitting variables
        SplitLine = new GRBConstr *[numNodes];                                            // multi-dimensional constraints for splitting variables
        SplitStg = new GRBConstr *[numNodes];                                            // multi-dimensional constraints for splitting variables
        // Convex = IloAdd(MP,IloRangeArray(mpEnv,N,1,1));									// convexity constraints in primary master problem
        Convex = new GRBConstr[numNodes];

        // the expansion in a scenario path should be at most once
        for (int n = 0; n < numNodes; n++) {
            const auto &node = scenarioNodes[n];
            if (node.period != paramReg->numPeriods - 1) continue;
            // upper level constraints
            for (int g = 0; g < numGen; g++) {
                GRBLinExpr tempExpr = 0;
                for (const int i : node.nodesInPath) {
                    tempExpr += var_xg[i][g];
                }
                MPmodel.addConstr(tempExpr, GRB_LESS_EQUAL, generators[g].exists() ? 0 : 1, "UpperLevel");
            }
            for (int l = 0; l < numLine; l++) {
                GRBLinExpr tempExpr = 0;
                for (const int i : node.nodesInPath) {
                    tempExpr += var_xl[i][l];
                }
                MPmodel.addConstr(tempExpr, GRB_LESS_EQUAL, lines[l].exists() ? 0 : 1, "UpperLevel");
            }
            for (int k = 0; k < numStg; k++) {
                GRBLinExpr tempExpr = 0;
                for (const int i : node.nodesInPath) {
                    tempExpr += var_xs[i][k];
                }
                MPmodel.addConstr(tempExpr, GRB_LESS_EQUAL, storages[k].exists() ? 0 : 1, "UpperLevel");
            }
        }

        // in one location only one potential generator can be built
#ifdef RestrictSameSiteExpansion
        for (int n = 0; n < numNodes; n++) {
            const auto &node = scenarioNodes[n];
            if (node.period != numPeriods - 1) continue;
            for (int b = 0; b < numBus; b++) {
                GRBLinExpr numGensInSite = 0;
                bool hasCandidateGens = false;
                for (int g = 0; g < numGen; g++) {
                    const auto &gen = generators[g];
                    if (gen.exists() || (gen.loc != b)) continue;
                    for (const int i : node.nodesInPath) {
                        numGensInSite += var_xg[i][g];
                    }
                    hasCandidateGens = true;
                }
                if (hasCandidateGens) {
                    MPmodel.addConstr(numGensInSite, GRB_LESS_EQUAL, 1);
                }
            }
        }
#endif

        /*	Formulate objecitve and constraints for primary master problem (initial: all expansions are made at node 0 and opn cost is bigM) */
        for (int n = 0; n < numNodes; n++) {
            const Node &node = scenarioNodes[n];
            auto &var_lambda = columns[n].back();
            std::sprintf(buf, "Convexity(%d)", n);
            Convex[n] = MPmodel.addConstr(var_lambda, GRB_EQUAL, 1, buf);
            SplitGen[n] = new GRBConstr[numGen];
            SplitLine[n] = new GRBConstr[numLine];
            SplitStg[n] = new GRBConstr[numStg];
            for (int g = 0; g < numGen; g++) {
                GRBLinExpr tempExpr = 0;
                for (const int i : node.nodesInPath) {
                    tempExpr += var_xg[i][g];
                }
                const int coef = 1 - int(generators[g].exists());
                std::sprintf(buf, "Gen(%d,%d)", n, g);
                SplitGen[n][g] = MPmodel.addConstr(tempExpr - coef * var_lambda, GRB_GREATER_EQUAL, 0, buf);
            }
            for (int l = 0; l < numLine; l++) {
                GRBLinExpr tempExpr = 0;
                for (const int i : node.nodesInPath) {
                    tempExpr += var_xl[i][l];
                }
                const int coef = 1 - int(lines[l].exists());
                std::sprintf(buf, "Line(%d,%d)", n, l);
                SplitLine[n][l] = MPmodel.addConstr(tempExpr - coef * var_lambda, GRB_GREATER_EQUAL, 0, buf);
            }
            for (int k = 0; k < numStg; k++) {
                GRBLinExpr tempExpr = 0;
                for (const int i : node.nodesInPath) {
                    tempExpr += var_xs[i][k];
                }
                const int coef = 1 - int(storages[k].exists());
                std::sprintf(buf, "Stg(%d,%d)", n, k);
                SplitStg[n][k] = MPmodel.addConstr(tempExpr - coef * var_lambda, GRB_GREATER_EQUAL, 0, buf);
            }

        }
    }
    catch (GRBException &e) {
        std::cerr << "Error while creating master problem; Error message: " << e.getMessage() << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Model::createUpperSubproblem()
{
	char buf[100];

	if (is_master_proc()) PRINT_SUBSECTION("create subproblems");

#ifdef SolverStreamOff
	SPmodel.set(GRB_IntParam_OutputFlag, 0);
#else
	if (!is_master_proc()) {
		SPmodel.set(GRB_IntParam_OutputFlag, 0);
	}
#endif

	const int n = rank;
	const Node &node = scenarioNodes[n];
	//SPmodel.set(GRB_DoubleParam_TimeLimit, paramReg->smp_timlmt);
	SPmodel.set(GRB_DoubleParam_MIPGap, paramReg->smp_gap);

	try {
        /* ********************** Decision variables ************************** */
        // capacity expansion
        var_zg = SPmodel.addVars(numGen, GRB_BINARY);
        var_zl = SPmodel.addVars(numLine, GRB_BINARY);
        var_zs = SPmodel.addVars(numStg, GRB_BINARY);
#ifdef LPMODEL
        for (int g = 0; g < numGen; g++) {
            std::sprintf(buf, "zg(%d,%d)", color, g);
            var_zg[g].set(GRB_StringAttr_VarName, buf);
        }
        for (int l = 0; l < numLine; l++) {
            std::sprintf(buf, "zl(%d,%d)", color, l);
            var_zl[l].set(GRB_StringAttr_VarName, buf);
        }
        for (int k = 0; k < numStg; k++) {
            std::sprintf(buf, "zs(%d,%d)", color, k);
            var_zs[k].set(GRB_StringAttr_VarName, buf);
        }
#endif // LPMODEL
        // unit commitment variables
        const int numLowerLevelIndices = Season::numSeasons * paramReg->numDmdScenarios * ucPeriods;
        varGenStatus = new GRBVar *[numLowerLevelIndices];                                            //binary variable for on/off status of plants alpha[n][d][k][g]
        varStartUp = new GRBVar *[numLowerLevelIndices];
        varPowerGen = new GRBVar *[numLowerLevelIndices];                                                //amount of power generated from generator p[n][d][k][g]
        varPowerFlow = new GRBVar **[numLowerLevelIndices];                                                //power flow between bus i and bus j f[n][d][k][l]
        varPowerWithdrawal = new GRBVar *[numLowerLevelIndices];                                        //total power withdrawn from storage facilities of bus u[n][d][k][sg]
        varPowerInject = new GRBVar *[numLowerLevelIndices];                                            //total power injected into storage facilities of bus v[n][d][k][sg]
        varPowerRemaining = new GRBVar *[numLowerLevelIndices];                                        //total remaining power in storage facilities of bus at the beginning of k v[n][d][k][sg]
        for (int season = 0; season != Season::numSeasons; ++season) {
           for (int d = 0; d < paramReg->numDmdScenarios; d++) {
                for (int h = 0; h < ucPeriods; h++) {
                    const int index = get2ndStageIndexInNode(season, d, h);
#ifdef PhaseII
                    varGenStatus[index] = SPmodel.addVars(numGen, GRB_BINARY);
#else
                    varGenStatus[index] = SPmodel.addVars(numGen);
#endif
                    varStartUp[index] = SPmodel.addVars(numGen);
                    varPowerGen[index] = SPmodel.addVars(numGen);
                    varPowerFlow[index] = new GRBVar *[numLine];
                    varPowerWithdrawal[index] = SPmodel.addVars(numStg);
                    varPowerInject[index] = SPmodel.addVars(numStg);
                    varPowerRemaining[index] = SPmodel.addVars(numStg);
                    for (int l = 0; l < numLine; l++) {
                        varPowerFlow[index][l] = SPmodel.addVars(Direction::numDirections);
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

        /* ************************* Objective **************************** */
        GRBLinExpr spobj_exp = 0;
        for (int season = 0; season != Season::numSeasons; ++season) {
           for (int d = 0; d < paramReg->numDmdScenarios; d++) {
                GRBLinExpr expr = 0;
                for (int h = 0; h < ucPeriods; h++) {
                    const int index = get2ndStageIndexInNode(season, d, h);
                    for (int g = 0; g < numGen; g++) {
                        const auto &gen = generators[g];
                        expr += gen.startUpCost * varStartUp[index][g] + gen.fuelCostB * gen.fuelPriceRate[n] * varPowerGen[index][g] +
                                gen.fuelCostA * varGenStatus[index][g];
                    }
                }
                spobj_exp += expr * paramReg->prob_dmd_scenario * get_opn_cost_factor(season);
            }
        }
        spobj_exp *= node.probability;
        SPmodel.setObjective(spobj_exp, GRB_MINIMIZE);

        // *********** Constraints in SP <<<--->>> UC & capacity limit constraints (GenAvlbty will be updated in the main function during the CG algorithm)
        for (int g = 0; g < numGen; g++) {
            SPmodel.addConstr(var_zg[g], GRB_LESS_EQUAL, 1 -
                                                         int(generators[g].exists()));                // at most one expansion along the scenario path
        }
        for (int l = 0; l < numLine; l++) {
            SPmodel.addConstr(var_zl[l], GRB_LESS_EQUAL, 1 - int(lines[l].exists()));
        }
        for (int k = 0; k < numStg; k++) {
            SPmodel.addConstr(var_zs[k], GRB_LESS_EQUAL, 1 - int(storages[k].exists()));

            // build energy storage only if there is a renewable energy plant in the same site
            GRBLinExpr capGenInSite = 0;
            for (int g = 0; g < numGen; g++) {
                if ((generators[g].loc != storages[k].loc) || !generators[g].is_renewable) continue;
                capGenInSite += int(generators[g].exists()) + var_zg[g];
            }
            SPmodel.addConstr(var_zs[k], GRB_LESS_EQUAL, capGenInSite);
        }

        // total capacity after expansion must exceed total demand
        GRBLinExpr totalCapacity = 0;
        for (int g = 0; g < numGen; g++) {
            const auto &gen = generators[g];
            totalCapacity += gen.capacity * (int(gen.exists()) + var_zg[g]);
        }
        for (int k = 0; k < numStg; k++) {
            const auto &stg = storages[k];
            totalCapacity += stg.capacity * (int(stg.exists()) + var_zs[k]);
        }

        double totalDemand = 0;
        for (int _s = 0; _s != Season::numSeasons; _s++) {
           for (int d = 0; d < paramReg->numDmdScenarios; d++) {
                for (int h = 0; h < ucPeriods; h++) {
                    double temp = 0;
                    for (int b = 0; b < numBus; b++) {
                        temp += buses[b].season_hourly_load_scenarios[get2ndStageIndex(n, _s, d, h)];
                    }
                    totalDemand = (temp > totalDemand) ? temp : totalDemand;
                }
            }
        }
        SPmodel.addConstr(totalCapacity, GRB_GREATER_EQUAL, totalDemand);

        for (int season = 0; season != Season::numSeasons; season++) {
           for (int d = 0; d < paramReg->numDmdScenarios; d++) {
                for (int h = 0; h < ucPeriods; h++) {
                    const int index = get2ndStageIndexInNode(season, d, h);
                    for (int g = 0; g < numGen; g++) {
                        const auto &gen = generators[g];
                        std::sprintf(buf, "GenExist(%d,%d,%d,%d)", season, d, h, g);
                        SPmodel.addConstr(varGenStatus[index][g] - int(gen.exists()) - var_zg[g], GRB_LESS_EQUAL, 0,
                                          buf);

                        std::sprintf(buf, "GenAct(%d,%d,%d)", season, h, g);
                        if (h == 0) {
                            SPmodel.addConstr(varStartUp[index][g] - varGenStatus[index][g], GRB_GREATER_EQUAL, 0, buf);
                        } else {
                            SPmodel.addConstr(varStartUp[index][g] - varGenStatus[index][g] +
                                              varGenStatus[get2ndStageIndexInNode(season, d, h - 1)][g],
                                              GRB_GREATER_EQUAL, 0, buf);
                        }

                        if (h) {
                            // min up and down
                            for (int tau = h; tau < ucPeriods && tau < h + generators[g].minUpTime; tau++) {
                                std::sprintf(buf, "MinUp(%d,%d,%d,%d)", season, h, tau, g);
                                SPmodel.addConstr(varGenStatus[index][g] -
                                                  varGenStatus[get2ndStageIndexInNode(season, d, h - 1)][g],
                                                  GRB_LESS_EQUAL,
                                                  varGenStatus[get2ndStageIndexInNode(season, d, tau)][g]);
                            }
                            for (int tau = h; tau < ucPeriods && tau < h + generators[g].minDownTime; tau++) {
                                std::sprintf(buf, "MinDown(%d,%d,%d,%d)", season, h, tau, g);
                                SPmodel.addConstr(varGenStatus[index][g] -
                                                  varGenStatus[get2ndStageIndexInNode(season, d, h - 1)][g],
                                                  GRB_GREATER_EQUAL,
                                                  varGenStatus[get2ndStageIndexInNode(season, d, tau)][g] - 1);
                            }
                        }

                        std::sprintf(buf, "PowUB(%d,%d,%d,%d)", season, d, h, g);
                        SPmodel.addConstr(varPowerGen[index][g] - gen.capacity * varGenStatus[index][g], GRB_LESS_EQUAL,
                                          0, buf);

                        if (gen.is_renewable) {
                            std::sprintf(buf, "RnGen(%d,%d,%d,%d)", n, season, h, g);
                            SPmodel.addConstr(varPowerGen[index][g], GRB_LESS_EQUAL,
                                              gen.supply_scenarios[get2ndStageIndex(n, season, d, h)], buf);
                        }

                        if (h) {
                            std::sprintf(buf, "RampUp(%d,%d,%d,%d)", n, season, h, g);
                            SPmodel.addConstr(
                                    varPowerGen[index][g] - varPowerGen[get2ndStageIndexInNode(season, d, h - 1)][g],
                                    GRB_LESS_EQUAL, generators[g].rampUp);
                            std::sprintf(buf, "RampDown(%d,%d,%d,%d)", n, season, h, g);
                            SPmodel.addConstr(
                                    varPowerGen[get2ndStageIndexInNode(season, d, h - 1)][g] - varPowerGen[index][g],
                                    GRB_LESS_EQUAL, generators[g].rampDown);
                        }
                    }
                    for (int l = 0; l < numLine; l++) {
                        const auto &line = lines[l];
                        std::sprintf(buf, "FlowLimit(%d,%d,%d,%d)", season, d, h, l);
                        SPmodel.addConstr(varPowerFlow[index][l][Direction::Forward] +
                                          varPowerFlow[index][l][Direction::Backward] - line.capacity * var_zl[l],
                                          GRB_LESS_EQUAL, line.currentCapacity, buf);
                    }
                    for (int k = 0; k < numStg; k++) {
                        const auto &storage = storages[k];
                        std::sprintf(buf, "Stg(%d,%d,%d,%d)", season, d, h,
                                     k);                                                                                            //Constraint 23
                        SPmodel.addConstr(varPowerRemaining[index][k] - storage.capacity * var_zs[k], GRB_LESS_EQUAL,
                                          storage.currentCapacity, buf);

                        std::sprintf(buf, "StgLmt(%d,%d,%d,%d)", season, d, h, k);
                        SPmodel.addConstr(varPowerWithdrawal[index][k] - varPowerRemaining[index][k], GRB_LESS_EQUAL, 0,
                                          buf);                                                        //(21)
                        std::sprintf(buf, "StgLmt(%d,%d,%d,%d)", season, d, h, k);
                        SPmodel.addConstr(varPowerRemaining[index][k], GRB_LESS_EQUAL, storage.capacity, buf);
                        std::sprintf(buf, "StgBal(%d,%d,%d,%d)", season, d, h, k);
                        if (h ==
                            0) {                                                                                                                        //in the 1st lower level period, reset the storage facility
                            SPmodel.addConstr(varPowerRemaining[index][k], GRB_EQUAL, 0, buf);
                        } else {
                            SPmodel.addConstr(varPowerRemaining[index][k] -
                                              varPowerRemaining[get2ndStageIndexInNode(season, d, h - 1)][k]
                                              - varPowerInject[get2ndStageIndexInNode(season, d, h - 1)][k] +
                                              varPowerWithdrawal[get2ndStageIndexInNode(season, d, h - 1)][k],
                                              GRB_EQUAL, 0, buf);                            //Constraint 22
                        }
                    }
                    for (int b = 0; b < numBus; b++) {
                        std::sprintf(buf, "Dmd(%d,%d,%d,%d)", season, d, h, b);
                        GRBLinExpr tempExpr = 0;
                        for (int l = 0; l < numLine; l++) {
                            if (lines[l].from == b) tempExpr -= varPowerFlow[index][l][Direction::Forward] -
                                  varPowerFlow[index][l][Direction::Backward] * (1 - lines[l].lossRatio);
                            if (lines[l].to == b) tempExpr += varPowerFlow[index][l][Direction::Forward] * (1 - lines[l].lossRatio) -
                                                               varPowerFlow[index][l][Direction::Backward];
                        }
                        for (int g = 0; g < numGen; g++) {
                            if (generators[g].loc == b) tempExpr += varPowerGen[index][g];
                        }
                        for (int k = 0; k < numStg; k++) {
                            if (storages[k].loc == b) tempExpr += varPowerWithdrawal[index][k] * paramReg->stgCeff -
                                                                  varPowerInject[index][k];
                        }
                        SPmodel.addConstr(tempExpr, GRB_GREATER_EQUAL,
                                          buses[b].season_hourly_load_scenarios[get2ndStageIndex(n, season, d, h)],
                                          buf);
                    }
                }
            }
        }
    }catch (GRBException &e) {
        std::cerr << "Error while creating subproblem:" << e.getMessage() << "---" << e.getErrorCode() << std::endl;
        exit(EXIT_FAILURE);
    }catch (const std::out_of_range& oor) {
        std::cerr << "Out of Range error: " << oor.what() << '\n';
    }catch(const std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Model::getMpDuals()
{	
	for (int n = 0; n < numNodes; n++) {														//obtain dual solution
		for (int tau = 0; tau <= scenarioNodes[n].period; tau++) {			
			for (int g = 0; g < numGen; g++)	DualGen[n][g] = SplitGen[n][g].get(GRB_DoubleAttr_Pi);
			for (int l = 0; l < numLine; l++)	DualLine[n][l] = SplitLine[n][l].get(GRB_DoubleAttr_Pi);
			for (int k = 0; k < numStg; k++)	DualStg[n][k] = SplitStg[n][k].get(GRB_DoubleAttr_Pi);
		}
		DualConvex[n] = Convex[n].get(GRB_DoubleAttr_Pi);
	}
}

Result Model::solveUpperSubproblem()
{
	//................All processors: solve subproblems.....................//
    for (int g = 0; g < numGen; g++)
        var_zg[g].set(GRB_DoubleAttr_Obj, DualGen[rank][g]);
    for (int l = 0; l < numLine; l++)
        var_zl[l].set(GRB_DoubleAttr_Obj, DualLine[rank][l]);
    for (int k = 0; k < numStg; k++)
        var_zs[k].set(GRB_DoubleAttr_Obj, DualStg[rank][k]);
#ifdef LPMODEL
	char buf[100];
	sprintf(buf, "%s/sp.lp", output_directory.c_str());
	SPmodel.write(buf);
#endif
	try{
		SPmodel.optimize();

		if (SPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
			std::cerr << "The subproblem " << rank << " is infeasible!" << std::endl;
			exit(EXIT_FAILURE);
		}

		ReducedCost[rank] = SPmodel.get(GRB_DoubleAttr_ObjBound) - DualConvex[rank];											// reduced cost
		SmpOFV[rank] = SPmodel.get(GRB_DoubleAttr_ObjVal);
		LambdaCoef[rank] = SPmodel.get(GRB_DoubleAttr_ObjVal);																							// operation cost -> coefficient of lambda in obj of MP
		SmpGap[rank] = SPmodel.get(GRB_DoubleAttr_MIPGap);
      PspTime[rank] = SPmodel.get(GRB_DoubleAttr_Runtime);

        // get columns from subproblems
        for (int g = 0; g < numGen; g++) {
            ColExpPlanGen[rank][g] = var_zg[g].get(GRB_DoubleAttr_X);
            LambdaCoef[rank] -= ColExpPlanGen[rank][g] * DualGen[rank][g];
        }
        for (int l = 0; l < numLine; l++){
            ColExpPlanLine[rank][l] = var_zl[l].get(GRB_DoubleAttr_X);
            LambdaCoef[rank] -= ColExpPlanLine[rank][l] * DualLine[rank][l];
        }
        for (int k = 0; k < numStg; k++) {
            ColExpPlanStg[rank][k] = var_zs[k].get(GRB_DoubleAttr_X);
            LambdaCoef[rank] -= ColExpPlanStg[rank][k] * DualStg[rank][k];
        }
	}
	catch (GRBException &e) {
		std::cerr << "Error while solving Subproblem " << color << "; Error message: " << e.getMessage() << std::endl;
		exit(EXIT_FAILURE);
	}

	return Result::OPTIMAL;
}

void Model::addColumns(bool nested)
{
	const int node_number = nested ? color : rank;
	const bool is_sender = nested ? ((!is_master_proc()) && is_sub_master_proc()) : (!is_master_proc());

	if (is_sender)
	{
		MPI_Send(ColExpPlanGen[node_number], numGen, MPI_DOUBLE, master, 0, MPI_COMM_WORLD);						// send columns to master processor
		MPI_Send(ColExpPlanLine[node_number], numLine, MPI_DOUBLE, master, 1, MPI_COMM_WORLD);
		MPI_Send(ColExpPlanStg[node_number], numStg, MPI_DOUBLE, master, 2, MPI_COMM_WORLD);

		MPI_Send(&ReducedCost[node_number], 1, MPI_DOUBLE, master, 3, MPI_COMM_WORLD);									// send reduced cost to master processor
		MPI_Send(&LambdaCoef[node_number], 1, MPI_DOUBLE, master, 4, MPI_COMM_WORLD);
		MPI_Send(&SmpOFV[node_number], 1, MPI_DOUBLE, master, 5, MPI_COMM_WORLD);
		MPI_Send(&SmpGap[node_number], 1, MPI_DOUBLE, master, 6, MPI_COMM_WORLD);
		MPI_Send(&PspTime[node_number], 1, MPI_DOUBLE, master, 7, MPI_COMM_WORLD);
		if(nested){
		   MPI_Send(&PspPhase2Time[node_number], 1, MPI_DOUBLE, master, 8, MPI_COMM_WORLD);
		   MPI_Send(&PspSmpTime[node_number], 1, MPI_DOUBLE, master, 9, MPI_COMM_WORLD);
		}
	}

	if(is_master_proc()) {
		for (int n = 0; n < numNodes; n++)
		{
			if (n == master) continue;

			const int sender = nested ? (n * numCoresPerNode) : n;

			MPI_Recv(ColExpPlanGen[n], numGen, MPI_DOUBLE, sender, 0, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);					// receive columns from workder processors
			MPI_Recv(ColExpPlanLine[n], numLine, MPI_DOUBLE, sender, 1, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(ColExpPlanStg[n], numStg, MPI_DOUBLE, sender, 2, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

			MPI_Recv(&ReducedCost[n], 1, MPI_DOUBLE, sender, 3, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(&LambdaCoef[n], 1, MPI_DOUBLE, sender, 4, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(&SmpOFV[n], 1, MPI_DOUBLE, sender, 5, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(&SmpGap[n], 1, MPI_DOUBLE, sender, 6, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			MPI_Recv(&PspTime[n], 1, MPI_DOUBLE, sender, 7, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

			if(nested){
			   MPI_Recv(&PspPhase2Time[n], 1, MPI_DOUBLE, sender, 8, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			   MPI_Recv(&PspSmpTime[n], 1, MPI_DOUBLE, sender, 9, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			}
		}

		// update master problem
		auto add_column_start = std::chrono::steady_clock::now();
		try {
            for (int n = 0; n < numNodes; n++) {
                columns[n].push_back(MPmodel.addVar(0, GRB_INFINITY, LambdaCoef[n], GRB_CONTINUOUS));
                auto &newVar = columns[n].back();
                // newVars[n].set(GRB_DoubleAttr_Obj, LambdaCoef[n]);
#ifdef LPMODEL
                char name_buf[30];
                std::sprintf(name_buf, "column(%d,%d)", n, static_cast<int>(columns[n].size() - 1));
                newVar.set(GRB_StringAttr_VarName, name_buf);
#endif
                MPmodel.chgCoeff(Convex[n], newVar, 1);
                for (int g = 0; g < numGen; g++) {
                    MPmodel.chgCoeff(SplitGen[n][g], newVar, -ColExpPlanGen[n][g]);
                }
                for (int l = 0; l < numLine; l++) {
                    MPmodel.chgCoeff(SplitLine[n][l], newVar,-ColExpPlanLine[n][l]);
                }
                for (int k = 0; k < numStg; k++) {
                    MPmodel.chgCoeff(SplitStg[n][k], newVar, -ColExpPlanStg[n][k]);
                }
            }
        }
        catch (GRBException &e) {
            std::cerr << "Error while adding columns! Error message: " << e.getMessage() << std::endl;
            exit(EXIT_FAILURE);
        }
        auto add_column_end = std::chrono::steady_clock::now();
        rankSolutionTime = double(std::chrono::duration_cast<std::chrono::seconds>(add_column_end - add_column_start).count());
	}
}

void Model::getMPbounds()
{
	if (!is_master_proc()) return;
	MpNewLB = MpUB;
	for (int n = 0; n < numNodes; n++) {
		if (ReducedCost[n] < 0) {
			MpNewLB += ReducedCost[n];
		}
	}
	if (MpNewLB > MpLB)	MpLB = MpNewLB;
	gap = abs((MpUB - MpLB) / MpLB);
}

Result Model::solveMasterProblem(bool nested)
{
	if (!is_master_proc()) return Result::OPTIMAL;
	PRINT_SECTION("Solve the Master Problem in the last Column Generation Iteration");
#ifdef LPMODEL
	char buf[100];
	const std::string cg_or_nested = nested ? "NestedMP" : "MP";
	std::sprintf(buf, "%s/%s.lp", output_directory.c_str(), cg_or_nested.c_str());
	MPmodel.write(buf);
#endif
	MPmodel.set(GRB_IntParam_Method, GRB_METHOD_BARRIER);
	MPmodel.set(GRB_IntParam_Crossover, 1);

	MPmodel.optimize();
	if (MPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
	{
	   return Result::INFEASIBLE;
	}
	MpLpObj = MPmodel.get(GRB_DoubleAttr_ObjVal);
	rankSolutionTime += MPmodel.get(GRB_DoubleAttr_Runtime);

	// enforce binary constraints
	for (int n = 0; n < numNodes; n++) {
		for (int g = 0; g < numGen; g++)	var_xg[n][g].set(GRB_CharAttr_VType, GRB_BINARY);
		for (int l = 0; l < numLine; l++)	var_xl[n][l].set(GRB_CharAttr_VType, GRB_BINARY);
		for (int k = 0; k < numStg; k++)	var_xs[n][k].set(GRB_CharAttr_VType, GRB_BINARY);
	}

	for(auto &nodeColumns : columns){
	   for(auto &itCol : nodeColumns){
	      itCol.set(GRB_CharAttr_VType, GRB_BINARY);
	   }
	}

	MPmodel.optimize();
	if (MPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
	{
		return Result::INFEASIBLE;
	}

	MpUB = MPmodel.get(GRB_DoubleAttr_ObjVal);
	rankSolutionTime += MPmodel.get(GRB_DoubleAttr_Runtime);
	getExpansionSolution();

	operationalCost = new double[numNodes];
	for (int n = 0; n < numNodes; n++) {
		operationalCost[n] = 0;
		if (alg == Algorithm::COLGEN || alg == Algorithm::NESTED) {
			for (const auto &column : columns[n]) {
				operationalCost[n] += column.get(GRB_DoubleAttr_X) * column.get(GRB_DoubleAttr_Obj);
			}
		}
	}

	return Result::OPTIMAL;
}




