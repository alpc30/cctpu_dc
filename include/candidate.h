/*
 * candidate.h
 *
 *  Created on: 26/05/2016
 *      Author: jing
 */

#ifndef CANDIDATE_H_
#define CANDIDATE_H_

#include "conflict.h"
#include "stnu.h"

namespace CDS
{
	class Candidate
	{
	public:
		std::vector<Conflict> v_resolved_conflicts;
		std::queue<int> v_unassigned_vars;
		std::map<int, int> v_assignments;
		vector<pair<int, int> > vec_ass;
		bool b_checked;

		// keep a record of the relaxed bounds:
		// 1. if the boolean variable is true, the value of the second double is
		// the relaxed bound.
		// 2. the index is mapped from the original stnu, so it only takes O(1)
		// to achieve the relaxed bound.
		std::vector<pair<bool, double> > v_relaxed_bounds;
		double value_preference;
		double reward_vars;

		bool b_solvable;

		bool operator<(const Candidate & a) const
		{
			if (value_preference != a.value_preference)
				return value_preference < a.value_preference;
			return reward_vars < a.reward_vars;
		}

		void print();
		void print(const stnu&a);

		void setPreference(const stnu &a);

		double getRelaxation(const LINK &a);
		double relaxedBound(const LINK &a, const bool &b) const;
		void setRelaxedBound(const LINK &a, const bool &b, const double &newb);
		void setMaxDelay(const stnu &a, const double &x);

		Candidate()
		{
			b_checked = false;
		}
		;
		Candidate(const stnu &a);
		Candidate(const stnu &a, vector<pair<int, int> > v_assignments);
	};

}

#endif /* CANDIDATE_H_ */
