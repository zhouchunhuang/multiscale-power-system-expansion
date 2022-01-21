#include <vector>
#include "model.h"

using namespace std;

int Model::optimize() 
{
	Result result = Result::OPTIMAL;
	switch (alg)
	{
	case Algorithm::DIRECT:
		result = direct();
		break;
	case Algorithm::EXPECT:
		result = expectedValForm();
		break;
	case Algorithm::COLGEN:
		result = columnGeneration();
		break;
	case Algorithm::NESTED:
		result = nestedDecomposition();
		break;
	case Algorithm::BENDERS:
		result = bendersDecomposition();
		break;
	case Algorithm::TEST:
		test();
		break;
	default:
		break;
	}
	return result;
}

void Model::constructScnStructure()
{
	scenarioNodes.clear();
	int node_id = 0;
	for (int t = 0; t < paramReg->numPeriods; t++) {
	   const int numNodesInPeriod = int(pow(paramReg->numRealizations, t));
		for (int n = 0; n < numNodesInPeriod; n++) {					
		   scenarioNodes.emplace_back(node_id++, t, int(paramReg->numScenarios / numNodesInPeriod) * n, int(paramReg->numScenarios / numNodesInPeriod) * (n + 1));
			// find the nodes in the scenario path
			auto &newNode = scenarioNodes.back();
			newNode.parent = (newNode.id == 0) ? -1 : getNode(int(paramReg->numScenarios / numNodesInPeriod) * n, t - 1)->id;
			newNode.nodesInPath.clear();
			for (int tau = 0; tau <= t; tau++) {
			   Node* nodeInPath = getNode(int(paramReg->numScenarios / numNodesInPeriod) * n, tau);
				newNode.nodesInPath.push_back(nodeInPath->id);
			}
		}
	}
	numNodes = static_cast<int>(scenarioNodes.size());
	const auto prob_scenario = double(1.0 / double(paramReg->numScenarios));

	for (auto &node : scenarioNodes) {
		node.probability = (node.end - node.start) * prob_scenario;
	}
}

void Model::getExpansionSolution()
{
	sol_xg = new double*[numNodes];
	sol_xl = new double*[numNodes];
	sol_xs = new double*[numNodes];
	for (int n = 0; n < numNodes; n++) {
		sol_xg[n] = MPmodel.get(GRB_DoubleAttr_X, var_xg[n], numGen);
		sol_xl[n] = MPmodel.get(GRB_DoubleAttr_X, var_xl[n], numLine);
		sol_xs[n] = MPmodel.get(GRB_DoubleAttr_X, var_xs[n], numStg);
	}

	investmentCost = 0;
	for (int n = 0; n < numNodes; n++) {
		const auto &node = scenarioNodes[n];
		for (int g = 0; g < numGen; g++) {
			const auto &gen = generators[g];
			investmentCost += node.probability * gen.cost[n] * sol_xg[n][g];
		}
		for (int l = 0; l < numLine; l++) {
			const auto &line = lines[l];
			investmentCost += node.probability * line.cost[n] * sol_xl[n][l];
		}
		for (int k = 0; k < numStg; k++) {
			const auto &stg = storages[k];
			investmentCost += node.probability * stg.cost[n] * sol_xs[n][k];
		}
	}

}

/*
void Model::outputSol()
{
	_end=clock();
	cmp_time=_end-_start;
	
#ifdef MODEL_EXPORT
	sprintf(buf,"%sColGenBS/cg_mp.opt.lp",outputDirectory.c_str());
	MPSolver.exportModel(buf);
#endif
#ifdef FILEOUTON
	MPSolver.solve();
	output<<"CPLEX computation time <seconds>:\t" << cmp_time/CLOCKS_PER_SEC <<endl;
	output << "Solution information..................................................."<<endl;
	output << "Solution status:\t" << MPSolver.getStatus() << endl;
	output << "OFV = \t" << MPSolver.getObjValue() << endl;
	output << "Gap = \t" << gap << endl;
	if((MPSolver.getStatus() != IloAlgorithm::Infeasible ) && (MPSolver.getStatus() != IloAlgorithm::InfeasibleOrUnbounded)) 
	{
		for(n=0;n<N;n++)
		{
			output<< "Solution of Node " << n <<endl;
			output<< "capacity expansion decision xg:\t";
			for(g=0;g<n_Gen;g++)
			{
				for(i=0;i<n_tech;i++)
					output<< MPSolver.getValue(v_xg[n][g][i])<<"\t";
			}
			output<<"\n";
					
			output<< "capacity expansion decision xl:\t";
			for(l=0;l<n_Line;l++)
			{
				for(i=0;i<n_tech;i++)
					output<< MPSolver.getValue(v_xl[n][l][i])<<"\t";
			}
			output<<"\n";

			output<< "capacity expansion decision xs:\t";
			for(sg=0;sg<n_Stg;sg++)
			{
				for(i=0;i<n_tech;i++)
					output<< MPSolver.getValue(v_xs[n][sg][i])<<"\t";
			}
			output<<"\n";
			output<<"\n";

		}			
	}
#endif	
}
*/



