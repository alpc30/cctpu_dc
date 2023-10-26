/*
 * conflict.h
 *
 *  Created on: 01/06/2016
 *      Author: jing
 */
#include "head.h"
#include "stnu.h"
#include <unordered_set>
#ifndef CONFLICT_H_
#define CONFLICT_H_

namespace CDS
{

	class Conflict
	{
	public:

		int id;
		//std::vector<DistanceGraphEdge> v_conflict_edge; // to be change to int
		std::vector<int> v_conflict_edge_id;
		std::vector<std::pair<int, int> > v_assignments;
		std::vector<vector<int> > v_alter_blocked_moats;
		//bool b_var;
//	std::vector<double> v_relaxed_bounds; // keep a record of the relaxed bounds
		// it seems the relaxation here is not useful

		Conflict()
		{
			id = -1;
			v_conflict_edge_id.clear();
			v_assignments.clear();
			//	v_relaxed_bounds.clear();
		}
		;

		bool operator<(const Conflict a) const
		{
			if (v_conflict_edge_id.size() != a.v_conflict_edge_id.size())
				return v_conflict_edge_id.size() < a.v_conflict_edge_id.size();

			if (v_assignments.size() != a.v_assignments.size())
				return v_assignments.size() < a.v_assignments.size();
			for (int i = 0; i < v_conflict_edge_id.size(); i++)
			{
				if (v_conflict_edge_id[i] != a.v_conflict_edge_id[i])
					return v_conflict_edge_id[i] < a.v_conflict_edge_id[i];
			}
			return false;
		}
		;
		void print(const stnu & a) const;
		void print() const;

		bool isSame(const Conflict a) const
		{
			if (v_conflict_edge_id.size() != a.v_conflict_edge_id.size())
				return false;

			if (v_assignments.size() != a.v_assignments.size())
				return false;
			for (int i = 0; i < v_conflict_edge_id.size(); i++)
			{
				if (v_conflict_edge_id[i] != a.v_conflict_edge_id[i])
					return false;
			}
			return true;
		}

	};

	class SimEdge
	{
	public:
		std::string start, end;
		//int ctg_id;
		int edge_id;
		bool b_ub; //0-lb, 1-ub
		std::vector<int> dis_vars;
		bool operator<(const SimEdge & a) const
		{
			if (start != a.start)
				return start < a.start;
			return end < a.end;
		}
	};

#define SIM_CONF_U 1 //conflict: sum(u_i)<=x
#define SIM_CONF_L 2 //conflict: sum(l_i)>=y
#define SIM_CONF_D 3 //conflict: sum(u_i-l_i)>=z (D=U|L)
#define SIM_CONF_O 0 //conflict: others (just in case)

	class SimConflict // a simplified version of conflict
	{
	public:
		SimEdge sim_edge;

		std::vector<SimEdge> v_variable;
		double lb;
		double olb;
		bool b_sign; //0-less, 1-greater
		int conf_type; // SIM_CONF_U/L/D
		int start_node; // the start node
		//std::vector<double> v_lb, v_ub; //prehistory + ctgi within[lbi,ubi]

		int ret_sign() const
		{
			if (v_variable[0].b_ub)
				return -1;
			return 1;
		}

		void print() const
		{
			printf("type %d: ", this->conf_type);
			for (int i = 0; i < v_variable.size(); i++)
			{
				if (v_variable[i].b_ub)
					printf(" - (U)");
				else
					printf(" + (L)");
				printf("%s->%s", v_variable[i].start.c_str(),
						v_variable[i].end.c_str());
			}
			//if(b_sign)
			//	printf(" >= ");
			//else
			printf(" >= ");
			printf("%lf\n", lb);
		}

		int compare_vars(const SimConflict &a) const
		{
			// if no variable is in sim-conflict(a), return 1
			// otherwise, return 0;
			unordered_set<int> v;
			for (int i = 0; i < a.v_variable.size(); ++i)
			{
				v.insert(a.v_variable[i].edge_id);
			}

			for (int i = 0; i < v_variable.size(); ++i)
			{
				if (v.find(v_variable[i].edge_id) != v.end())
					return 0;
			}

			return 1;
		}

		int compare(const SimConflict &a) const
		{
			// return 0 or 1 if a cannot cover the current SimConflict

			// if the variables are different, return 0;
			// else if all variables are L/U, if lb<a.lb return 1, else return -1;
			// else if lb(L)<=lb(U), return -1, else return 1.

			if (a.v_variable.size() != v_variable.size())
				return 0;

			for (int i = 0; i < v_variable.size(); i--)
			{
				if (v_variable[i] < a.v_variable[i])
					return 0;
				if (a.v_variable[i] < v_variable[i])
					return 0;
			}		//broaden compare

			double lb1 = lb;
			double lb2 = a.lb;

			if (v_variable[0].b_ub != a.v_variable[0].b_ub)		//lb vs ub
			{
				if (ret_sign() < 0) //ub
				{
					if (lb * ret_sign() >= a.lb)
						return -1;
				}
				else
				{
					if (lb <= a.lb * a.ret_sign())
						return -1;
				}
				return 1;
				//if(a.lb*a.ret_sign() <= b.lb*b.ret_sign())
				//return true;
			}

			if (lb < a.lb)
				return 1;
			else
				return -1;
		}

		void removeCtg(int cid);
		bool isNew(const vector<SimConflict> a) const;
	};

}

#endif /* CONFLICT_H_ */
