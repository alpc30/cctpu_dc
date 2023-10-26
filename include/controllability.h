/*
 * controllability.h
 *
 *  Created on: 07/09/2015
 *      Author: jing
 */

#ifndef CONTROLLABILITY_H_
#define CONTROLLABILITY_H_

#include "distance_graph.h"

namespace CDS
{

	class Controllability
	{
		stnu stnu_;
		DistanceGraph distance_graph_;
		long long run_time_;
		vector<int> backprop_depth; //0-init, 1-terminate successfully, other-terminate with NC
		std::vector<double> v_result_;
		std::vector<int> v_neg_edge_id;
		vector<bool> being_back; //the node is being back propagated
		struct depth_edge
		{
			int depth;
			DistanceEdge edge;
			bool operator<(const depth_edge& a) const
			{
				return depth < a.depth;
			}
			depth_edge(int d, DistanceEdge e)
			{
				depth = d;
				edge = e;
			}
		};
		vector<vector<depth_edge> > v_q_neg_edge;
		vector<DistanceEdge> ng_path; // the path the ng (used for cubic dc checking)
		vector<Conflict> v_conf;
		vector<int> v_src_nodes; //for bellmanFord if the graph has more than one connected parts
	public:

		Controllability();

		//checking DC or SC
		void InitDistanceGraph(const stnu& a);
		stnu InitDistanceGraph(const stnu& a, const Candidate & b);
		bool Checking(stnu current_stnu);
		bool CheckingSC();
		bool CheckingDC();
		void edgeEdge(LINK a);

		//DC checking propagation and negative cycle
		std::vector<Conflict> Propagation(DistanceEdge edge); //return relaxable disjunctions
		std::vector<Conflict> PropagationDijkstra(DistanceEdge edge); // alternative function of propagation
		std::vector<Conflict> BellmanFord(); //return negative cycles
		int isTighterOrdinaryEdge(DistanceEdge edge);
		int isTighterUpperEdge(DistanceEdge edge);
		bool addOrdinaryEdge(DistanceEdge edge);
		bool addUpperEdge(DistanceEdge edge);
		bool containNegLoop(DistanceEdge edge);
		bool isLabelConsistent(const DistanceEdge &e1, const DistanceEdge &e2);

		//DC conflict extract
		std::vector<Conflict> DCConflict(stnu & s);
		std::vector<Conflict> DCConflict(const stnu& a, const Candidate & b);

		std::vector<Conflict> DCConflict1(const stnu& a, const Candidate & b);
//	std::vector<Conflict> DCConflict(const stnu &n, std::vector<double> v_relax, std::vector<int> v_relax_id);
		std::vector<Conflict> CheckingDCNegCycle();
		std::vector<Conflict> CheckingDCNegCycle1();
		Conflict getNegCycle(std::vector<DistanceEdge> v_e);
		Conflict getBlockedRelaxation(DistanceEdge e);

		//Cubic DC checking algorithim (TODO: varify)
		bool CubicCheckingDC();
		bool DCbackpropCheck(DistanceEdge e);

		//other functions
		long long get_run_time();
		void printResult(FILE *fp);

		~Controllability();

		// extract envelope based on the reduction rules
		void extractEnvelope(vector<Conflict>& v_neg_cycle);
		int DCbackprop(DistanceEdge e, vector<Conflict> &a,
		//vector<DistanceEdge > &b,
				const int &depth);
		void completeNegCycle(DistanceEdge e, vector<Conflict> &a,
				vector<DistanceEdge>&b, DistanceEdge& dis);

		bool dfs(int id, vector<DistanceEdge>& dis, vector<bool>& v,
				const DistanceEdge &source, vector<Conflict> &neg_cycle
				//, vector<DistanceEdge > &b
				);
		int get_id_backprop(const DistanceEdge & edge);

	};

} /* namespace CDS */
#endif /* CONTROLLABILITY_H_ */
