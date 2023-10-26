/*
 * cds.h
 *
 *  Created on: 21/09/2015
 *      Author: jing
 *    Function: main process of conflict directed search
 */

#ifndef CDS_H_
#define CDS_H_

#include "controllability.h"
#include "stnu.h"
//#include "lp.h"
#include "candidate.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

namespace CDS
{

	class SubSolution
	{
	public:
		//a class to describe the dynamic controllable solution for a STNU
		bool b_solved;
		bool b_solvable;
		//std::vector<std::vector<SimConflict> >v_sim_conflicts;
		std::vector<Conflict> v_conflicts;
		std::vector<pair<int, int> > v_value_dis_var;
		std::vector<pair<bool, double> > v_relaxed_bounds;
		void print();
		void simplifyConflicts(const stnu & a); //transfer v_conflicts to v_sim_conflicts
		bool isCovered(const SimConflict &a, const vector<SimConflict> &b);

		SubSolution()
		{
			b_solved = false;
			b_solvable = true;
			this->v_union_conj_conflicts.clear();
		}

		//a new structure of the sub-solution
		//if it works well, it will replace the old version
		//new structure is a union of conjunctions
		std::vector<std::vector<SimConflict> > v_union_conj_conflicts;
		std::vector<std::vector<SimConflict> > v_current_conflicts;
		std::vector<double> v_lb, v_ub; //prehistory + ctgi within[lbi,ubi]

	};

	class Cds
	{
		//the main class of search engine

		std::priority_queue<Candidate> q_candidate;
		std::vector<Conflict> v_conflicts;
		std::vector<Conflict> v_unresolved_conflicts;
		std::set<double> candidate_value;
		//std::vector<int> v_value_dis_var;
		std::vector<SubSolution> v_sub_solutions;
		stnu whole_stnu;
		int n_combines;
		int n_leaves;

		struct LBUB
		{
			double lb, ub;

			bool operator==(LBUB a) const
			{
				if (this->lb == a.lb && this->ub == a.ub)
					return true;
				return false;
			}
			bool operator<(LBUB a) const
			{
				if (this->lb != a.lb)
					return (this->lb) < a.lb;
				return this->ub < a.ub;
			}
		};

	public:
		std::vector<std::set<LBUB> > v_cover_ctg;
		std::vector<std::set<LBUB> > v_cover_rqm;
		// keep record of layer
		int mlayer, tlayer;

		//void init_v_cover();
		long long run_time;
		Cds();

		void setWholeSTNU(const stnu & t);
		void expandOnVariables(Candidate &a);
		//bool CombineSTNU(const stnu & a);
		virtual ~Cds();
		//int CDRU(stnu & a);
		bool DCCheckCSTNU(const stnu & a);
		void CDRU1(const stnu &a);
		void CDRU2(const Candidate & a, SubSolution & sol);
		void expandOnConflict(const Candidate & a,
				const std::vector<Conflict> c);
		void expandOnConflict1(const Candidate & a,
				const std::vector<Conflict> c);
		void CreateCandidateOnConflict(const Candidate & a, const Conflict & b);

		//Tree-search structure, return dc or not
		int DCCheckCSTNU1(const stnu & a);
		stnu getSTNU(const stnu &a, std::vector<pair<int, int> > assignments);
		std::vector<pair<int, int> > getNextAssignment(
				std::vector<pair<int, int> > v);
		bool combine2STNU(SubSolution& a, SubSolution b, int v_id);
		bool combineSubSol(SubSolution &a, SubSolution b, int v_id);
		//std::vector<SimConflict>  getUnion(std::vector<SimConflict> a, std::vector<SimConflict> b, int v);
		void updateSimConflict(SubSolution &a, int vid);
		bool updateSingleSimConflict(SimConflict &a, int vid);
		bool upOneSimConfbyCtg(SimConflict &a, int cid);
		bool upSimConfbyCtg(vector<vector<SimConflict> >&res,
				vector<vector<SimConflict> > a, int vid);
		//bool isInfeasible(const SimConflict &a);
		void removeRedundant(SubSolution &a, SubSolution &b);
		bool isRedundant(const SimConflict &a, const SubSolution &b);

		bool checkSubSolution(SubSolution &a);
		bool checkCombConflict(vector<SimConflict> &a);

		//A traditional CDRU which only expands on variables but not combining process
		//int DCCheckNoCombine(const stnu & a);

		//rewrite the dc checking process cause the former one is not correct
		Candidate CdsCCTPU(const stnu &a);

		std::vector<Conflict> resolveConflicts(Candidate & a);
		int resolved(const Conflict & c, Candidate &a);
		bool isComplete(const Candidate &a);

		string result;
		void generateResult(const Candidate &a);
		void generateResult(const int &x);

		bool tryCombine(SubSolution& a, SubSolution b, int var_id);
		bool checkPre(SubSolution b, int en);
	};

} /* namespace CDS */
#endif /* CDS_H_ */
