#include "candidate.h"

namespace CDS
{

Candidate::Candidate(const stnu& a)
{
//	this->stnu_ = a;

	b_checked = false;

	for (int i = 0; i < a.n_vars; i++)
	{
		this->v_unassigned_vars.push(a.order_dv[i]);
		//value_preference += a.v_dis_var[i].m_utility;
	}

	for (int i = 0; i < 2 * (a.n_ctg + a.n_rqm); i++)
	{
		v_relaxed_bounds.push_back(make_pair(false, 0));
	}

	if (a.s_info.id_obj == O_MAX_DELAY)
	{
		setPreference(a);
		reward_vars = a.utility;
	}
	else
	{
		if(a.s_info.id_obj == O_RELAX_COST)
		{
			this->value_preference = a.utility;
			this->reward_vars = a.utility;
		}
		else
		{
			value_preference = 0;
			reward_vars = 0;
		}
	}
}
;

Candidate::Candidate(const stnu& a, vector<pair<int, int> > v_assignments) :
		Candidate(a)
{
	//this->value_preference = a.utility;
	//this->reward_vars = a.utility;

	for (int i = 0; i < 2 * (a.n_ctg + a.n_rqm); i++)
	{
		v_relaxed_bounds[i] = (make_pair(false, 0));
	}
	vec_ass = v_assignments;
}
;

void Candidate::print()
{
	printf("Candidate: %.3lf (%.3lf)\n", value_preference, reward_vars);
	printf("Assignments(%d):", this->v_assignments.size());
	map<int, int>::iterator it;
	for (it = v_assignments.begin(); it != v_assignments.end(); it++)
		printf(" v%d=%d", it->first, it->second);
	printf("\n");
	printf("Resolved Conflicts: %d\n", this->v_resolved_conflicts.size());
	for (int i = 0; i < this->v_resolved_conflicts.size(); i++)
	{
		v_resolved_conflicts[i].print();
	}
}

void Candidate::print(const stnu&a)
{
	printf("Candidate: %.3lf (%.3lf) %d\n", value_preference, reward_vars,
			b_checked);
	printf("Assignments(%d):", this->v_assignments.size());
	map<int, int>::iterator it;
	for (it = v_assignments.begin(); it != v_assignments.end(); it++)
		printf(" v%d=%d", it->first, it->second);
	printf("\n");
	LINK t = a.ctg[0];
	t.print();
	printf("Resolved Conflicts: %d\n", this->v_resolved_conflicts.size());
	for (int i = 0; i < this->v_resolved_conflicts.size(); i++)
	{
		v_resolved_conflicts[i].print(a);
	}
}

void Candidate::setPreference(const stnu& a)
{
	double x = g_inf;
	if (a.s_info.id_obj == O_MAX_DELAY)
	{
		for (int i = 0; i < a.n_ctg; i++)
		{
			double dis;

			dis = relaxedBound(a.ctg[i], true) - relaxedBound(a.ctg[i], false);
			//double dis = a.ctg[i].ub-a.ctg[i].lb;
			x = x < dis ? x : dis;
		}
	}
	else if (a.s_info.id_obj == O_MAX_EARLINESS)
	{
		for (int i = 0; i < a.n_ctg; i++)
		{
			double dis = a.ctg[i].ub - a.ctg[i].lb;
			x = x < dis ? x : dis;
		}
	}
	else if (a.s_info.id_obj == O_DC_CHECK)
	{
		x = reward_vars;
		double relax = 0;
		for (int i = 0; i < a.n_ctg; i++)
		{
			relax += getRelaxation(a.ctg[i]);
		}
		x -= relax;
	}
	else if (a.s_info.id_obj % 10== O_RELAX_COST)
	{
		if(a.s_info.id_obj == O_RELAX_COST_NR)
			x=0;
		else x = reward_vars;
		double relax = 0;
		for (int i = 0; i < a.n_ctg; i++)
		{
			relax += getRelaxation(a.ctg[i]);
		}
		for (int i = 0; i < a.n_rqm; i++)
		{
			relax -= getRelaxation(a.rqm[i]);
		}

		x -= relax;
	}

	value_preference = x;

}

double Candidate::getRelaxation(const LINK &a)
{
	int id_relax = a.id_lr;
	double res = 0;
	if (id_relax < v_relaxed_bounds.size())
	{
		if (v_relaxed_bounds[id_relax].first)
		{
			res += (v_relaxed_bounds[id_relax].second - a.lb) * a.lb_cost_r;
		}
	}
	else
	{
		printf("error: relax_id exceeds the size of v_relaxed_bounds.\n");
		exit(0);
	}

	id_relax = a.id_ur;
	if (id_relax < v_relaxed_bounds.size())
	{
		if (v_relaxed_bounds[id_relax].first)
		{
			res += (a.ub - v_relaxed_bounds[id_relax].second) * a.ub_cost_r;
		}
	}
	else
	{
		printf("error: relax_id exceeds the size of v_relaxed_bounds.\n");
		exit(0);
	}
	return res;
}

double Candidate::relaxedBound(const LINK &a, const bool &b_ub) const
{
	int id_relax = a.id_lr;
	double res = 0;
	if (b_ub)
	{
		id_relax = a.id_ur;
		if (id_relax < v_relaxed_bounds.size())
		{
			if (v_relaxed_bounds[id_relax].first)
			{
				res = v_relaxed_bounds[id_relax].second;
			}
			else
				res = a.ub;
		}
		else
		{
			printf("error: relax_id exceeds the size of v_relaxed_bounds.\n");
			exit(0);
		}
	}
	else
	{
		if (id_relax < v_relaxed_bounds.size())
		{
			if (v_relaxed_bounds[id_relax].first)
			{
				res = v_relaxed_bounds[id_relax].second;
			}
			else
				res = a.lb;
		}
		else
		{
			printf("error: relax_id exceeds the size of v_relaxed_bounds.\n");
			exit(0);
		}
	}
	return res;
}

void Candidate::setMaxDelay(const stnu&a, const double &x)
{
	for (int i = 0; i < a.n_ctg; i++)
	{
		setRelaxedBound(a.ctg[i], true, a.ctg[i].lb + x);
	}
}

void Candidate::setRelaxedBound(const LINK &a, const bool &b, const double &nb)
{
	if (b)
		v_relaxed_bounds[a.id_ur] = make_pair(true, nb);
	else
		v_relaxed_bounds[a.id_lr] = make_pair(true, nb);
}
}
