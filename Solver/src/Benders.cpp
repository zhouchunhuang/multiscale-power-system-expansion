#include "model.h"


Result Model::bendersDecomposition() 
{
    if (is_master_proc())
    {
        PRINT_SECTION("Solving the Problem with Benders Decomposition in Parallel");
        PRINT_SUBSECTION("Initializing the algorithm");
        std::cout.precision(4);
    }

	initializeBenders();

	createBendersMasterProblem();

	createSubproblem();

   bendersLoop(0);

   bendersLoop(1);

#ifdef PhaseII
	message = 0;

	MpUB = paramReg->bigM;
	MpLB = -paramReg->bigM;
	MpNewUB = 0;
	getRecFuncLBInNodes();

	SPmodel.set(GRB_DoubleParam_MIPGap, paramReg->sp_gap);
	bendersLoop(2);
#endif

	solveMasterProblem();

	return Result::OPTIMAL;
}

void Model::bendersLoop(int phase) {
	if (is_master_proc()) {
		PRINT_SUBSECTION("Benders Decomposition Phase " + std::to_string(phase));
#ifdef StreamAlgInfo
        std::cout << "#\t\tLB\t";
        for (int n = 0; n < numNodes; n++)  std::cout << n << "\t";
        std::cout << "NewUB\tUB\tGap/Time" << std::endl;
#endif
	}

	message = 0;

	MpUB = paramReg->bigM;
	MpLB = -paramReg->bigM;
	MpNewUB = 0;

	if(is_master_proc()){
	   if(phase == 0){
	      for (int n = 0; n < numNodes; n++) {
	         for (int g = 0; g < numGen; g++) {
	            var_xg[n][g].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
	         }
	         for (int l = 0; l < numLine; l++) {
	            var_xl[n][l].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
	         }
	         for (int k = 0; k < numStg; k++) {
	            var_xs[n][k].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
	         }
	      }
	   }else{
	      for (int n = 0; n < numNodes; n++) {
	         for (int g = 0; g < numGen; g++) {
	            var_xg[n][g].set(GRB_CharAttr_VType, GRB_BINARY);
	         }
	         for (int l = 0; l < numLine; l++) {
	            var_xl[n][l].set(GRB_CharAttr_VType, GRB_BINARY);
	         }
	         for (int k = 0; k < numStg; k++) {
	            var_xs[n][k].set(GRB_CharAttr_VType, GRB_BINARY);
	         }
	      }
	   }
	}

	auto _start = std::chrono::steady_clock::now();
	iteration = 0;
	while (true) {
		iteration++;

		solveBendersMP();

		if (MPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
			std::cerr << "No feasible solution exists for the problem!" << std::endl;
			exit(EXIT_FAILURE);
		}

		for (int n = 0; n < numNodes; n++) {
			MPI_Bcast(CapGenInNodes[n], numGen, MPI_DOUBLE, master, MPI_COMM_WORLD);
			MPI_Bcast(CapLineInNodes[n], numLine, MPI_DOUBLE, master, MPI_COMM_WORLD);
			MPI_Bcast(CapStgInNodes[n], numStg, MPI_DOUBLE, master, MPI_COMM_WORLD);
		}

		for (int sp = 0; sp < numSPsPerCore; sp++) {	// every CPU core solves a set of subproblems in series;
			solveBendersSubproblem(sp, phase);
		}

		dataTransfer();

        /*if (!is_master_proc()) {
            for (int sp = 0; sp < numSPsPerCore; sp++) {
                const int noSpInNode = get_subproblem_no_in_node(sp);
                MPI_Send(SecDualGenInNodes[color][noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, master, 10000 * noSpInNode + 0, MPI_COMM_WORLD);
                MPI_Send(SecDualRampUpInNodes[color][noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, master, 10000 * noSpInNode + 1, MPI_COMM_WORLD);
                MPI_Send(SecDualRampDownInNodes[color][noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, master, 10000 * noSpInNode + 2, MPI_COMM_WORLD);
                MPI_Send(SecDualLineInNodes[color][noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numLine, MPI_DOUBLE, master, 10000 * noSpInNode + 3, MPI_COMM_WORLD);
                MPI_Send(SecDualStgInNodes[color][noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numStg, MPI_DOUBLE, master, 10000 * noSpInNode + 4, MPI_COMM_WORLD);
                MPI_Send(SecDualDmdInNodes[color][noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numBus, MPI_DOUBLE, master, 10000 * noSpInNode + 5, MPI_COMM_WORLD);
                MPI_Send(SecDualMinDownInNodes[color][noSpInNode], int(MinDownConstrs.size()), MPI_DOUBLE, master, 10000 * noSpInNode + 6, MPI_COMM_WORLD);
                MPI_Send(&SPstatusInNodes[color][noSpInNode], 1, MPI_INT, master, 10000 * noSpInNode + 7, MPI_COMM_WORLD);
                if (SPstatusInNodes[color][noSpInNode] == SpStatus::Optimal)
                {
                    MPI_Send(&recFuncValInNodes[color][noSpInNode], 1, MPI_DOUBLE, master, 10000 * noSpInNode + 8, MPI_COMM_WORLD);
                    MPI_Send(&SPofvInNodes[color][noSpInNode], 1, MPI_DOUBLE, master, 10000 * noSpInNode + 9, MPI_COMM_WORLD);
                }
            }
        }
        else {
            for (int n = 0; n < numNodes; n++) {
                for (int i = 0; i < numSPsPerNode; i++) {
                    const int sender = n * numCoresPerNode + int(i / numSPsPerCore);
                    if (sender == master) {
                        continue;
                    }
                    MPI_Recv(SecDualGenInNodes[n][i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE,
                             sender, 10000 * i + 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(SecDualRampUpInNodes[n][i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE,
                             sender, 10000 * i + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(SecDualRampDownInNodes[n][i], numDmdScenariosPerSubproblem * ucPeriods * numGen,
                             MPI_DOUBLE, sender, 10000 * i + 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(SecDualLineInNodes[n][i], numDmdScenariosPerSubproblem * ucPeriods * numLine, MPI_DOUBLE,
                             sender, 10000 * i + 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(SecDualStgInNodes[n][i], numDmdScenariosPerSubproblem * ucPeriods * numStg, MPI_DOUBLE,
                             sender, 10000 * i + 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(SecDualDmdInNodes[n][i], numDmdScenariosPerSubproblem * ucPeriods * numBus, MPI_DOUBLE,
                             sender, 10000 * i + 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(SecDualMinDownInNodes[n][i], int(MinDownConstrs.size()), MPI_DOUBLE, sender, 10000 * i + 6,
                             MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&SPstatusInNodes[n][i], 1, MPI_INT, sender, 10000 * i + 7, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                    if (SPstatusInNodes[n][i] == SpStatus::Optimal) {
                        MPI_Recv(&recFuncValInNodes[n][i], 1, MPI_DOUBLE, sender, 10000 * i + 8, MPI_COMM_WORLD,
                                 MPI_STATUS_IGNORE);
                        MPI_Recv(&SPofvInNodes[n][i], 1, MPI_DOUBLE, sender, 10000 * i + 9, MPI_COMM_WORLD,
                                 MPI_STATUS_IGNORE);
                    }
                }
            }
        }*/

		if(is_master_proc()){
			getBendersMpBounds();

			addBendersCuts(phase);

			auto _end = std::chrono::steady_clock::now();
			cmp_time = double(std::chrono::duration_cast<std::chrono::seconds>(_end - _start).count());
			const double mp_gap = std::abs((MpUB - MpLB) / MpUB);

			if ((mp_gap <= paramReg->mp_gap) || (cmp_time >= paramReg->time_lmt))
			{
				message = 1;
			}

			if((phase == 1) && (cmp_time > paramReg->time_lmt / 2) && (MpUB < BIGM)){
			   message = 1;
			}


#ifdef StreamAlgInfo
            std::cout << "-" << iteration << "-\t" << MpLB << "(" << MPmodel.get(GRB_DoubleAttr_Runtime) << ")\t";
            for (int n = 0; n < numNodes; n++){
                double opnCostInNode = 0;
                for(int i = 0; i < numSPsPerNode; i++){
                    opnCostInNode += SPofvInNodes[n][i];
                }
                std::cout << opnCostInNode << "\t";
            }
            std::cout << MpNewUB << "\t" << MpUB << "\t" << mp_gap << "\t" << cmp_time << std::endl;
#endif
		}
		/* ******************************* */
		MPI_Bcast(&message, 1, MPI_INT, master, MPI_COMM_WORLD);
		if (message == 1) break;
	}
}

void Model::initializeBenders() 
{
	message = 0;

	MpUB = BIGM;
	MpLB = -BIGM;
	MpNewUB = 0;

	numDmdScenariosPerSubproblem = int(paramReg->numDmdScenarios / paramReg->numSPsPerSeason);
    numSPsPerNode = Season::numSeasons * paramReg->numSPsPerSeason;

    if (num_ranks % numNodes) {
        std::cerr << "Incorrect parallel computing setting: numCoresPerNode not integer!" << std::endl;
        std::cerr << "Number of Ranks: " << num_ranks << "; number of scenario nodes: " << numNodes << std::endl;
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

#ifdef TEST_MODE_2
	if(is_master_proc()) {
        std::cout << "numDmdScenarios = " << numDmdScenarios << std::endl;
        std::cout << "numSPsPerSeason = " << paramReg->numSPsPerSeason << std::endl;
        std::cout << "numDmdScenariosPerSubproblem = " << numDmdScenariosPerSubproblem << std::endl;
        std::cout << "numSPsPerNode = " << numSPsPerNode << std::endl;
        std::cout << "numCoresPerNode = " << numCoresPerNode << std::endl;
        std::cout << "numSPsPerCore = " << numSPsPerCore << std::endl;
    }
    exit(EXIT_SUCCESS);
#endif
}

void Model::createBendersMasterProblem()
{	
	if (!is_master_proc()) return;		// create master problem only in the master processor
	PRINT_SUBSECTION("Create Master Problem");
#ifdef LPMODEL	
	char buf[100];
#endif

#ifdef SolverStreamOff
	MPmodel.set(GRB_IntParam_OutputFlag, 0);
#endif

	MPmodel.set(GRB_DoubleParam_MIPGap, paramReg->mp_gap);
	/* *********************** Variables *********************** */
	var_xg = new GRBVar*[numNodes];
	var_xl = new GRBVar*[numNodes];
	var_xs = new GRBVar*[numNodes];
	var_y = new GRBVar*[numNodes];		
	for (int n = 0; n < numNodes; n++) {
		var_xg[n] = MPmodel.addVars(numGen, GRB_BINARY);
		var_xl[n] = MPmodel.addVars(numLine, GRB_BINARY);
		var_xs[n] = MPmodel.addVars(numStg, GRB_BINARY);
		var_y[n] = MPmodel.addVars(numSPsPerNode);
#ifdef LPMODEL	
		for (int g = 0; g < numGen; g++) {
			std::sprintf(buf, "xg(%d,%d)", n, g);
			var_xg[n][g].set(GRB_StringAttr_VarName, buf);
		}
		for (int l = 0; l < numLine; l++) {
			std::sprintf(buf, "xl(%d,%d)", n, l);
			var_xl[n][l].set(GRB_StringAttr_VarName, buf);
		}
		for (int k = 0; k < numStg; k++) {
			std::sprintf(buf, "xs(%d,%d)", n, k);
			var_xs[n][k].set(GRB_StringAttr_VarName, buf);
		}
		for (int i = 0; i < numSPsPerNode; i++) {
			std::sprintf(buf, "y(%d)", n);
			var_y[n][i].set(GRB_StringAttr_VarName, buf);
		}
#endif
	}
		
	/* *********************** Objective *********************** */
	GRBLinExpr xpr = 0;
	for (int n = 0; n < numNodes; n++) {
		const auto &node = scenarioNodes[n];
		for (int g = 0; g < numGen; g++) {
			xpr += node.probability * generators[g].cost[n] * var_xg[n][g];
		}
		for (int l = 0; l < numLine; l++) {
			xpr += node.probability * lines[l].cost[n] * var_xl[n][l];
		}
		for (int k = 0; k < numStg; k++) {
			xpr += node.probability * storages[k].cost[n] * var_xs[n][k];
		}
		for (int i = 0; i < numSPsPerNode; i++) {
			xpr += var_y[n][i];
		}		
	}
	MPmodel.setObjective(xpr, GRB_MINIMIZE);

	/* *********************** Constraints *********************** */
	for (int n = 0; n < numNodes; n++) {
		const auto &node = scenarioNodes[n];
		if (node.period != paramReg->numPeriods - 1) continue;
		// upper level constraints
		for (int g = 0; g < numGen; g++) {
			GRBLinExpr tempXpr = 0;
			for (const int j : node.nodesInPath) {
			   tempXpr += var_xg[j][g];
			}
			MPmodel.addConstr(tempXpr, GRB_LESS_EQUAL, generators[g].exists() ? 0 : 1, "UpperLevel");
		}
		for (int l = 0; l < numLine; l++) {
		   GRBLinExpr tempXpr = 0;
			for (const int j : node.nodesInPath) {
			   tempXpr += var_xl[j][l];
			}
			MPmodel.addConstr(tempXpr, GRB_LESS_EQUAL, lines[l].exists() ? 0 : 1, "UpperLevel");
		}
		for (int k = 0; k < numStg; k++) {
		   GRBLinExpr tempXpr = 0;
			for (const int j : node.nodesInPath) {
			   tempXpr += var_xs[j][k];
			}
			MPmodel.addConstr(tempXpr, GRB_LESS_EQUAL, storages[k].exists() ? 0 : 1, "UpperLevel");
		}		
	}

	// total generation capacity must exceed the total demand in each demand scenario
	for (int n = 0; n < numNodes; n++) {
		const auto &node = scenarioNodes[n];
		GRBLinExpr totalGenCapacity = 0;
		double totalDemand = 0;
		for (int g = 0; g < numGen; g++) {
			totalGenCapacity += generators[g].currentCapacity;
			for (const int j : node.nodesInPath) {
				totalGenCapacity += var_xg[j][g] * generators[g].capacity;
			}
		}
		for (int k = 0; k < numStg; k++) {
			totalGenCapacity += storages[k].currentCapacity;
			for (const int j : node.nodesInPath) {
				totalGenCapacity += var_xs[j][k] * storages[k].capacity;
			}
		}
		for (int s = 0; s < Season::numSeasons; s++) {
		   for (int d = 0; d < paramReg->numDmdScenarios; d++) {
				for (int h = 0; h < ucPeriods; h++) {
					double temp = 0;
					for (int b = 0; b < numBus; b++) {
						temp += buses[b].season_hourly_load_scenarios[get2ndStageIndex(n, s, d, h)];
					}
					totalDemand = totalDemand > temp ? totalDemand : temp;
				}
			}
		}
		MPmodel.addConstr(totalGenCapacity, GRB_GREATER_EQUAL, totalDemand);
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
    for (int n = 0; n < numNodes; n++) {
	    const auto &node = scenarioNodes[n];
        for (int g = 0; g < numGen; g++) {
            GRBLinExpr cumulativeExpansion = 0;
            for (auto j : node.nodesInPath) {
                cumulativeExpansion += var_xg[j][g];
            }
            MPmodel.addConstr(cumulativeExpansion, GRB_LESS_EQUAL, 1 -
                                                                   int(generators[g].exists()));                // at most one expansion along the scenario path
        }
        for (int l = 0; l < numLine; l++) {
            GRBLinExpr cumulativeExpansion = 0;
            for (auto j : node.nodesInPath) {
                cumulativeExpansion += var_xl[j][l];
            }
            MPmodel.addConstr(cumulativeExpansion, GRB_LESS_EQUAL, 1 - int(lines[l].exists()));
        }
        for (int k = 0; k < numStg; k++) {
            GRBLinExpr cumulativeExpansion = 0;
            for (auto j : node.nodesInPath) {
                cumulativeExpansion += var_xs[j][k];
            }
            MPmodel.addConstr(cumulativeExpansion, GRB_LESS_EQUAL, 1 - int(storages[k].exists()));

            // build energy storage only if there is a renewable energy plant in the same site
            GRBLinExpr capGenInSite = 0;
            for (int g = 0; g < numGen; g++) {
                if ((generators[g].loc != storages[k].loc) || !generators[g].is_renewable) continue;
                capGenInSite += int(generators[g].exists());
                for (auto j : node.nodesInPath) {
                    capGenInSite += var_xg[j][g];
                }
            }
            MPmodel.addConstr(cumulativeExpansion, GRB_LESS_EQUAL, capGenInSite);
        }
    }

	// reserve margin requirement
	for (int n = 0; n < numNodes; n++){
	   const auto &node = scenarioNodes[n];
	   GRBLinExpr totalGenCap = 0;
	   for (int g = 0; g < numGen; g++) {
	      const auto &gen = generators[g];
	      totalGenCap += generators[g].currentCapacity;
         for (const int j : node.nodesInPath) {
            totalGenCap += var_xg[j][g] * generators[g].capacity;
         }
	   }
	   double rmReq = 0;
	   for (int b = 0; b < numBus; b++) {
	      rmReq += buses[b].peak_loads_in_scenario_nodes[n];
	   }
	   MPmodel.addConstr(totalGenCap, GRB_GREATER_EQUAL, rmReq * paramReg->reserveMarginRate);
	}
}

void Model::getRecFuncLBInNodes()
{
	// use maximum possible capacity
	for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
		for (int h = 0; h < ucPeriods; h++) {
			for (int g = 0; g < numGen; g++) {
				GenLimitConstrs[d * ucPeriods * numGen + h * numGen + g].set(GRB_DoubleAttr_RHS, generators[g].capacity);
			}
			for (int l = 0; l < numLine; l++) {
				LineLimitConstrs[d * ucPeriods * numLine + h * numLine + l].set(GRB_DoubleAttr_RHS, lines[l].capacity);
			}
			for (int k = 0; k < numStg; k++) {
				StgLimitConstrs[d * ucPeriods * numStg + h * numStg + k].set(GRB_DoubleAttr_RHS, storages[k].capacity);
			}
		}
	}

	SPmodel.optimize();
	if (SPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
		std::cerr << "SP with max capacity is infeasible!\n";
		return;
	}

	recFuncLBInNodes[color][sub_rank] = SPmodel.get(GRB_DoubleAttr_ObjVal);

	MPI_Bcast(&message, 1, MPI_INT, master, MPI_COMM_WORLD);
	if (!is_master_proc()) {
		for (int sp = 0; sp < numSPsPerCore; sp++) {
			const int noSpInNode = get_subproblem_no_in_node(sp);
			MPI_Send(&recFuncLBInNodes[color][noSpInNode], 1, MPI_DOUBLE, master, 100, MPI_COMM_WORLD);
		}
	}
	else {
		for (int n = 0; n < numNodes; n++) {
			for (int i = 0; i < numSPsPerNode; i++) {
				const int sender = n * numCoresPerNode + int(i / numSPsPerCore);
				if (sender == master) {
					continue;
				}			
				MPI_Recv(&recFuncLBInNodes[n][i], 1, MPI_DOUBLE, sender, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			}
		}
	}

}

Result Model::solveBendersMP()
{
	if (!is_master_proc())
		return Result::OPTIMAL;

    std::vector<GRBConstr> tempConstrs;

	if(iteration == 1){
        for(int n = 0; n < numNodes; n++) {
            const auto &node = scenarioNodes[n];
            const int t = node.period;
            for (int g = 0; g < numGen; g++) {
                tempConstrs.push_back(MPmodel.addConstr(var_xg[n][g], GRB_EQUAL,
                                                         int((t == 0) && !(generators[g].exists()))));
            }
            for (int l = 0; l < numLine; l++) {
                tempConstrs.push_back(
                        MPmodel.addConstr(var_xl[n][l], GRB_EQUAL, int((t == 0) && !(lines[l].exists()))));
                for (int k = 0; k < numStg; k++) {
                    tempConstrs.push_back(
                            MPmodel.addConstr(var_xs[n][k], GRB_EQUAL, int((t == 0) && !(storages[k].exists()))));
                }
            }
        }
	}

	MPmodel.optimize();

	if (MPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) { 
		char buf[100];
		std::sprintf(buf, "%s/benders_mp_infeasible.lp", output_directory.c_str());
		MPmodel.write(buf);
		std::cerr << "The problem is infeasible!";
		exit(EXIT_FAILURE);
	}

	if(iteration > 1) {
        MpLB = MPmodel.get(GRB_DoubleAttr_ObjBound) > MpLB ? MPmodel.get(GRB_DoubleAttr_ObjBound) : MpLB;
    }
	MpNewUB = MPmodel.get(GRB_DoubleAttr_ObjVal);
	for (int n = 0; n < numNodes; n++) {
		for (int i = 0; i < numSPsPerNode; i++) {
			MpNewUB -= var_y[n][i].get(GRB_DoubleAttr_X);		// will be calculated further(+recourse funciton values from subproblems)
		}
	}

	// optimal solution from current master problem
	for (int n = 0; n < numNodes; n++) {
		const auto &node = scenarioNodes[n];
		for (int g = 0; g < numGen; g++) {
			CapGenInNodes[n][g] = generators[g].currentCapacity;
			for (const auto j : node.nodesInPath) {
				CapGenInNodes[n][g] += var_xg[j][g].get(GRB_DoubleAttr_X) * generators[g].capacity;
			}			
		}

		for (int l = 0; l < numLine; l++) {
			CapLineInNodes[n][l] = lines[l].currentCapacity;
			for (const auto j : node.nodesInPath) {
				CapLineInNodes[n][l] += var_xl[j][l].get(GRB_DoubleAttr_X) * lines[l].capacity;
			}
		}
		for (int k = 0; k < numStg; k++) {
			CapStgInNodes[n][k] = storages[k].currentCapacity;
			for (const auto j : node.nodesInPath) {
				CapStgInNodes[n][k] += var_xs[j][k].get(GRB_DoubleAttr_X) * storages[k].capacity;
			}
		}
	}

	if(iteration == 1){
        /* Remove temp constraints */
        if (is_master_proc()) {
            for (auto it : tempConstrs) {
                MPmodel.remove(it);
            }
        }
        tempConstrs.clear();
	}

	return Result::OPTIMAL;
}

Result Model::solveBendersSubproblem(int sp/*no.SubproblemInCore*/, int phase)
{
	const int noSpInNode = get_subproblem_no_in_node(sp);		// no. subproblem in a node
	const int season = get_season_of_subproblem(noSpInNode);
    const auto& node = scenarioNodes[color];
	// update subproblem
	for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
		for (int h = 0; h < ucPeriods; h++) {
			for (int g = 0; g < numGen; g++) {
				const auto &gen = generators[g];
				GenLimitConstrs[d * ucPeriods * numGen + h * numGen + g].set(GRB_DoubleAttr_RHS, CapGenInNodes[color][g]);
				SPmodel.chgCoeff(GenLimit2Constrs[d * ucPeriods * numGen + h * numGen + g], varGenStatus[getLowerIndexInSP(d, h)][g],
					gen.is_renewable ? (-gen.supply_scenarios[getLowerIndex(color, season, get_dmd_scenario(d, noSpInNode), h)]) : (-gen.capacity));
			}
			for (int l = 0; l < numLine; l++) {
				LineLimitConstrs[d * ucPeriods * numLine + h * numLine + l].set(GRB_DoubleAttr_RHS, CapLineInNodes[color][l]);
			}
			for (int k = 0; k < numStg; k++) {
				StgLimitConstrs[d * ucPeriods * numStg + h * numStg + k].set(GRB_DoubleAttr_RHS, CapStgInNodes[color][k]);
			}
			for (int b = 0; b < numBus; b++) {
				DemandConstrs[d * ucPeriods * numBus + h * numBus + b].set(GRB_DoubleAttr_RHS,
					buses[b].season_hourly_load_scenarios[getLowerIndex(color, season, get_dmd_scenario(d, noSpInNode), h)]);
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
                spobj_exp += node.probability * paramReg->prob_dmd_scenario * (gen.startUpCost * varStartUp[index][g] + gen.fuelCostA * varGenStatus[index][g]) +
                      node.probability * paramReg->prob_dmd_scenario * gen.fuelCostB * gen.fuelPriceRate[color] * varPowerGen[index][g];
            }
            for (int b = 0; b < numBus; b++) {
               spobj_exp += node.probability * paramReg->prob_dmd_scenario * paramReg->unmetDmdPnlty * varUnmetDemand[index][b];
            }
        }
    }
    SPmodel.setObjective(get_opn_cost_factor(season) * spobj_exp, GRB_MINIMIZE);

	/* ************************** Solve SP linear relaxation ************************************** */
	for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
		for (int h = 0; h < ucPeriods; h++) {
			for (int g = 0; g < numGen; g++) {
				varGenStatus[getLowerIndexInSP(d, h)][g].set(GRB_CharAttr_VType, GRB_CONTINUOUS);
				varGenStatus[getLowerIndexInSP(d, h)][g].set(GRB_DoubleAttr_UB, GRB_INFINITY);
			}
		}
	}

	SPmodel.optimize();

	//get Duals or Extreme rays
	if (SPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
		SPstatusInNodes[color][noSpInNode] = SpStatus::LRinfeasible;

        for (int i = 0; i < numDmdScenariosPerSubproblem * ucPeriods * numGen; i++) {
            SecDualGenInNodes[color][noSpInNode][i] = -abs(GenLimitConstrs[i].get(GRB_DoubleAttr_FarkasDual));
            SecDualRampUpInNodes[color][noSpInNode][i] = -abs(RampUpConstrs[i].get(GRB_DoubleAttr_FarkasDual));
            SecDualRampDownInNodes[color][noSpInNode][i] = -abs(RampDownConstrs[i].get(GRB_DoubleAttr_FarkasDual));
        }
        for (int i = 0; i < numDmdScenariosPerSubproblem * ucPeriods * numLine; i++) {
            SecDualLineInNodes[color][noSpInNode][i] = -abs(LineLimitConstrs[i].get(GRB_DoubleAttr_FarkasDual));
        }
        for (int i = 0; i < numDmdScenariosPerSubproblem * ucPeriods * numStg; i++) {
            SecDualStgInNodes[color][noSpInNode][i] = -abs(StgLimitConstrs[i].get(GRB_DoubleAttr_FarkasDual));
        }
        for (int i = 0; i < numDmdScenariosPerSubproblem * ucPeriods * numBus; i++) {
            SecDualDmdInNodes[color][noSpInNode][i] = abs(DemandConstrs[i].get(GRB_DoubleAttr_FarkasDual));
        }
        int constr_index = 0;
        for (const auto &it : MinDownConstrs) {
            SecDualMinDownInNodes[color][noSpInNode][constr_index++] = -abs(it.get(GRB_DoubleAttr_FarkasDual));
        }
        return Result::INFEASIBLE; //if the LR is infeasible, no need to solve SP
	}

    for (int i = 0; i < numDmdScenariosPerSubproblem * ucPeriods * numGen; i++) {
        SecDualGenInNodes[color][noSpInNode][i] = GenLimitConstrs[i].get(GRB_DoubleAttr_Pi);
        SecDualRampUpInNodes[color][noSpInNode][i] = RampUpConstrs[i].get(GRB_DoubleAttr_Pi);
        SecDualRampDownInNodes[color][noSpInNode][i] = RampDownConstrs[i].get(GRB_DoubleAttr_Pi);
    }
    for (int i = 0; i < numDmdScenariosPerSubproblem * ucPeriods * numLine; i++) {
        SecDualLineInNodes[color][noSpInNode][i] = LineLimitConstrs[i].get(GRB_DoubleAttr_Pi);
    }
    for (int i = 0; i < numDmdScenariosPerSubproblem * ucPeriods * numStg; i++) {
        SecDualStgInNodes[color][noSpInNode][i] = StgLimitConstrs[i].get(GRB_DoubleAttr_Pi);
    }
    for (int i = 0; i < numDmdScenariosPerSubproblem * ucPeriods * numBus; i++) {
        SecDualDmdInNodes[color][noSpInNode][i] = DemandConstrs[i].get(GRB_DoubleAttr_Pi);
    }
    int constr_index = 0;
    for (const auto &it : MinDownConstrs) {
        SecDualMinDownInNodes[color][noSpInNode][constr_index++] = it.get(GRB_DoubleAttr_Pi);
    }

    SPstatusInNodes[color][noSpInNode] = SpStatus::LRoptimal;
    SPofvInNodes[color][noSpInNode] = SPmodel.get(GRB_DoubleAttr_ObjVal);
    recFuncValInNodes[color][noSpInNode] = SPmodel.get(GRB_DoubleAttr_ObjVal);
	
	if (phase > 1) {
		/* ************************** Solve SP MIP model ************************************** */
		for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
			for (int h = 0; h < ucPeriods; h++) {
				for (int g = 0; g < numGen; g++) {
					varGenStatus[getLowerIndexInSP(d, h)][g].set(GRB_CharAttr_VType, GRB_BINARY);
				}
			}
		}
		SPmodel.optimize();

		if (SPmodel.get(GRB_IntAttr_Status) == GRB_INFEASIBLE) {
			/*std::cerr << "The Subproblem " << color << " subrank " << sub_rank << " is infeasible while its LR is feasible!" << std::endl;
			char buf[100];
			std::sprintf(buf, "%s/sp%d-%d-%d.lp", output_directory.c_str(), color, season, noSpInNode);
			SPmodel.write(buf);
			exit(EXIT_FAILURE);*/
			SPstatusInNodes[color][noSpInNode] = SpStatus::MIPinfeasible;

			return Result::INFEASIBLE;
		}
		SPstatusInNodes[color][noSpInNode] = SpStatus::Optimal;
		SPofvInNodes[color][noSpInNode] = SPmodel.get(GRB_DoubleAttr_ObjVal);
		recFuncValInNodes[color][noSpInNode] = SPmodel.get(GRB_DoubleAttr_ObjVal);
	}

	return Result::OPTIMAL;
}

void Model::dataTransfer()
{
    // transfer data in sub communication world
    if (!is_sub_master_proc()) {
        for (int sp = 0; sp < numSPsPerCore; sp++) {
            const int noSpInNode = get_subproblem_no_in_node(sp);
            MPI_Send(SecDualGenInNodes[color][noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 0, sub_comm);
            MPI_Send(SecDualRampUpInNodes[color][noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 1, sub_comm);
            MPI_Send(SecDualRampDownInNodes[color][noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 2, sub_comm);
            MPI_Send(SecDualLineInNodes[color][noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numLine, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 3, sub_comm);
            MPI_Send(SecDualStgInNodes[color][noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numStg, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 4, sub_comm);
            MPI_Send(SecDualDmdInNodes[color][noSpInNode], numDmdScenariosPerSubproblem * ucPeriods * numBus, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 5, sub_comm);
            MPI_Send(SecDualMinDownInNodes[color][noSpInNode], int(MinDownConstrs.size()), MPI_DOUBLE, sub_master, 10000 * noSpInNode + 6, sub_comm);
            MPI_Send(&SPstatusInNodes[color][noSpInNode], 1, MPI_INT, sub_master, 10000 * noSpInNode + 7, sub_comm);
            if ((SPstatusInNodes[color][noSpInNode] == SpStatus::Optimal) || (SPstatusInNodes[color][noSpInNode] == SpStatus::LRoptimal))
            {
                MPI_Send(&recFuncValInNodes[color][noSpInNode], 1, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 8, sub_comm);
                MPI_Send(&SPofvInNodes[color][noSpInNode], 1, MPI_DOUBLE, sub_master, 10000 * noSpInNode + 9, sub_comm);
            }
        }
    }
    else {
        for (int i = 0; i < numSPsPerNode; i++) {
            const int sender = int(i / numSPsPerCore);
            if (sender == sub_master) {
                continue;
            }
            MPI_Recv(SecDualGenInNodes[color][i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, sender,
                     10000 * i + 0, sub_comm, MPI_STATUS_IGNORE);
            MPI_Recv(SecDualRampUpInNodes[color][i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE,
                     sender, 10000 * i + 1, sub_comm, MPI_STATUS_IGNORE);
            MPI_Recv(SecDualRampDownInNodes[color][i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE,
                     sender, 10000 * i + 2, sub_comm, MPI_STATUS_IGNORE);
            MPI_Recv(SecDualLineInNodes[color][i], numDmdScenariosPerSubproblem * ucPeriods * numLine, MPI_DOUBLE,
                     sender, 10000 * i + 3, sub_comm, MPI_STATUS_IGNORE);
            MPI_Recv(SecDualStgInNodes[color][i], numDmdScenariosPerSubproblem * ucPeriods * numStg, MPI_DOUBLE, sender,
                     10000 * i + 4, sub_comm, MPI_STATUS_IGNORE);
            MPI_Recv(SecDualDmdInNodes[color][i], numDmdScenariosPerSubproblem * ucPeriods * numBus, MPI_DOUBLE, sender,
                     10000 * i + 5, sub_comm, MPI_STATUS_IGNORE);
            MPI_Recv(SecDualMinDownInNodes[color][i], int(MinDownConstrs.size()), MPI_DOUBLE, sender, 10000 * i + 6,
                     sub_comm, MPI_STATUS_IGNORE);
            MPI_Recv(&SPstatusInNodes[color][i], 1, MPI_INT, sender, 10000 * i + 7, sub_comm, MPI_STATUS_IGNORE);
            if ((SPstatusInNodes[color][i] == SpStatus::Optimal) || (SPstatusInNodes[color][i] == SpStatus::LRoptimal)) {
                MPI_Recv(&recFuncValInNodes[color][i], 1, MPI_DOUBLE, sender, 10000 * i + 8, sub_comm,
                         MPI_STATUS_IGNORE);
                MPI_Recv(&SPofvInNodes[color][i], 1, MPI_DOUBLE, sender, 10000 * i + 9, sub_comm,
                         MPI_STATUS_IGNORE);
            }
        }
    }

    // transfer data to master rank
    if ((!is_master_proc()) && (is_sub_master_proc())) {
        for (int i = 0; i < numSPsPerNode; i++) {
            MPI_Send(SecDualGenInNodes[color][i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, master, 10000 * i + 0, MPI_COMM_WORLD);
            MPI_Send(SecDualRampUpInNodes[color][i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, master, 10000 * i + 1, MPI_COMM_WORLD);
            MPI_Send(SecDualRampDownInNodes[color][i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, master, 10000 * i + 2, MPI_COMM_WORLD);
            MPI_Send(SecDualLineInNodes[color][i], numDmdScenariosPerSubproblem * ucPeriods * numLine, MPI_DOUBLE, master, 10000 * i + 3, MPI_COMM_WORLD);
            MPI_Send(SecDualStgInNodes[color][i], numDmdScenariosPerSubproblem * ucPeriods * numStg, MPI_DOUBLE, master, 10000 * i + 4, MPI_COMM_WORLD);
            MPI_Send(SecDualDmdInNodes[color][i], numDmdScenariosPerSubproblem * ucPeriods * numBus, MPI_DOUBLE, master, 10000 * i + 5, MPI_COMM_WORLD);
            MPI_Send(SecDualMinDownInNodes[color][i], int(MinDownConstrs.size()), MPI_DOUBLE, master, 10000 * i + 6, MPI_COMM_WORLD);
            MPI_Send(&SPstatusInNodes[color][i], 1, MPI_INT, master, 10000 * i + 7, MPI_COMM_WORLD);
            if ((SPstatusInNodes[color][i] == SpStatus::Optimal) || (SPstatusInNodes[color][i] == SpStatus::LRoptimal))
            {
                MPI_Send(&recFuncValInNodes[color][i], 1, MPI_DOUBLE, master, 10000 * i + 8, MPI_COMM_WORLD);
                MPI_Send(&SPofvInNodes[color][i], 1, MPI_DOUBLE, master, 10000 * i + 9, MPI_COMM_WORLD);
            }
        }
    }
    else if(is_master_proc()){
        for(int n = 0; n < numNodes; n++){
            if(n == master) continue;
            const int sender = n * numCoresPerNode;
            for (int i = 0; i < numSPsPerNode; i++) {
                MPI_Recv(SecDualGenInNodes[n][i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE, sender,
                         10000 * i + 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(SecDualRampUpInNodes[n][i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE,
                         sender, 10000 * i + 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(SecDualRampDownInNodes[n][i], numDmdScenariosPerSubproblem * ucPeriods * numGen, MPI_DOUBLE,
                         sender, 10000 * i + 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(SecDualLineInNodes[n][i], numDmdScenariosPerSubproblem * ucPeriods * numLine, MPI_DOUBLE,
                         sender, 10000 * i + 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(SecDualStgInNodes[n][i], numDmdScenariosPerSubproblem * ucPeriods * numStg, MPI_DOUBLE, sender,
                         10000 * i + 4, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(SecDualDmdInNodes[n][i], numDmdScenariosPerSubproblem * ucPeriods * numBus, MPI_DOUBLE, sender,
                         10000 * i + 5, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(SecDualMinDownInNodes[n][i], int(MinDownConstrs.size()), MPI_DOUBLE, sender, 10000 * i + 6,
                         MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&SPstatusInNodes[n][i], 1, MPI_INT, sender, 10000 * i + 7, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if ((SPstatusInNodes[n][i] == SpStatus::Optimal) || (SPstatusInNodes[n][i] == SpStatus::LRoptimal)) {
                    MPI_Recv(&recFuncValInNodes[n][i], 1, MPI_DOUBLE, sender, 10000 * i + 8, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                    MPI_Recv(&SPofvInNodes[n][i], 1, MPI_DOUBLE, sender, 10000 * i + 9, MPI_COMM_WORLD,
                             MPI_STATUS_IGNORE);
                }
            }
        }
    }

}

void Model::getBendersMpBounds()
{
	if (!is_master_proc()) return;
	for (int n = 0; n < numNodes; n++) {
		for (int i = 0; i < numSPsPerNode; i++) {
			if ((SPstatusInNodes[n][i] == SpStatus::LRinfeasible) || (SPstatusInNodes[n][i] == SpStatus::MIPinfeasible)) {
				MpNewUB = +BIGM;
				break;
			}
			MpNewUB += SPofvInNodes[n][i];
		}
	}

	MpUB = (MpNewUB < MpUB) ? MpNewUB : MpUB;
}

void Model::addBendersCuts(int phase)
{
   // add combinatorial cut
   zEqOne = 0;
   combnXpr = 0;

   for (int n = 0; n < numNodes; n++) {
      for (int g = 0; g < numGen; g++) {
         if (var_xg[n][g].get(GRB_DoubleAttr_X) > 0.5) {
            zEqOne++;
            combnXpr += var_xg[n][g];
         }
         else {
            combnXpr -= var_xg[n][g];
         }
      }
      for (int l = 0; l < numLine; l++) {
         if (var_xl[n][l].get(GRB_DoubleAttr_X) > 0.5) {
            zEqOne++;
            combnXpr += var_xl[n][l];
         }
         else {
            combnXpr -= var_xl[n][l];
         }
      }
      for (int k = 0; k < numStg; k++) {
         if (var_xs[n][k].get(GRB_DoubleAttr_X) > 0.5) {
            zEqOne++;
            combnXpr += var_xs[n][k];
         }
         else {
            combnXpr -= var_xs[n][k];
         }
      }
   }

   bool addCombCut = false;
   for (int n = 0; n < numNodes; n++) {
      for (int i = 0; i < numSPsPerNode; i++) {
         if(SPstatusInNodes[n][i] == SpStatus::MIPinfeasible){
            MPmodel.addConstr(combnXpr, GRB_LESS_EQUAL, zEqOne - 1, "IntegerInfCut");
            addCombCut = true;
            break;
         }
      }
      if(addCombCut) break;
   }

   // benders cuts
	for (int n = 0; n < numNodes; n++) {
		const auto &node = scenarioNodes[n];
		for (int i = 0; i < numSPsPerNode; i++) {
            const int season = get_season_of_subproblem(i);
            GRBLinExpr cut1Lhs = var_y[n][i];
            GRBLinExpr cut1Rhs = 0;
            for (int d = 0; d < numDmdScenariosPerSubproblem; d++) {
               for (int h = 0; h < ucPeriods; h++) {
                  for (int g = 0; g < numGen; g++) {
                     const auto &gen = generators[g];
                     cut1Rhs += SecDualGenInNodes[n][i][d * ucPeriods * numGen + h * numGen + g] * gen.currentCapacity;
                     for (const auto j : node.nodesInPath) {
                        cut1Rhs += SecDualGenInNodes[n][i][d * ucPeriods * numGen + h * numGen + g] * var_xg[j][g] * gen.capacity;
                     }

                     cut1Rhs += SecDualRampUpInNodes[n][i][d * ucPeriods * numGen + h * numGen + g] *
                        ((h == 0) ? gen.capacity : gen.rampUp);

                     cut1Rhs += SecDualRampDownInNodes[n][i][d * ucPeriods * numGen + h * numGen + g] *
                        ((h == 0) ? gen.capacity : gen.rampDown);
                  }
                  for (int l = 0; l < numLine; l++) {
                     cut1Rhs += SecDualLineInNodes[n][i][d * ucPeriods * numLine + h * numLine + l] * lines[l].currentCapacity;
                     for (const auto j : node.nodesInPath) {
                        cut1Rhs += SecDualLineInNodes[n][i][d * ucPeriods * numLine + h * numLine + l] * var_xl[j][l] * lines[l].capacity;
                     }
                  }
                  for (int k = 0; k < numStg; k++) {
                     cut1Rhs += SecDualStgInNodes[n][i][d * ucPeriods * numStg + h * numStg + k] * storages[k].currentCapacity;
                     for (const auto j : node.nodesInPath) {
                        cut1Rhs += SecDualStgInNodes[n][i][d * ucPeriods * numStg + h * numStg + k] * var_xs[j][k] * storages[k].capacity;
                     }
                  }
                  for (int b = 0; b < numBus; b++) {
                     cut1Rhs += SecDualDmdInNodes[n][i][d * ucPeriods * numBus + h * numBus + b] *
                        buses[b].season_hourly_load_scenarios[get2ndStageIndex(n, season, get_dmd_scenario(d, i), h)];
                  }
               }
            }

            int constr_index = 0;
            for (const auto &it : MinDownConstrs) {
               cut1Rhs += SecDualMinDownInNodes[n][i][constr_index++];	// RHS of min down constraints are always 1
            }

            switch(SPstatusInNodes[n][i])
            {
               case SpStatus::LRinfeasible:
                  MPmodel.addConstr(0, GRB_GREATER_EQUAL, cut1Rhs, "BendersInfCut");
                  break;
               case SpStatus::LRoptimal:
                  MPmodel.addConstr(cut1Lhs, GRB_GREATER_EQUAL, cut1Rhs, "BendersOptCut");
                  break;
               case SpStatus::Optimal:
                  MPmodel.addConstr(cut1Lhs, GRB_GREATER_EQUAL, cut1Rhs, "BendersOptCut");
                  MPmodel.addConstr(var_y[n][i], GRB_GREATER_EQUAL,
                                    (recFuncValInNodes[n][i] - recFuncLBInNodes[n][i]) * (combnXpr - zEqOne + 1)
                                    + recFuncLBInNodes[n][i], "LShapedCut");
                  break;
               default:
                  break;
            }
		}
	}
}