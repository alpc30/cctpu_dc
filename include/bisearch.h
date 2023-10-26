/*
 * bisearch.h
 *
 *  Created on: 16/05/2017
 *      Author: jing
 */

#ifndef BISEARCH_H_
#define BISEARCH_H_

#include "stnu.h"
#include "cds.h"

class BiSearch
{
public:
	/* do binary search */
	BiSearch()
	{
		p_stnu.n_ctg = -1, low = 0;
		high = 1e20;
		run_time = 0;
	}
	;
	~BiSearch()
	{
	}
	;
	CDS::stnu p_stnu;
	string result;
	double low, high;
	double run_time;

	BiSearch(CDS::stnu p) :
			BiSearch()
	{
		p_stnu = p;
	}
	;

	void bisearch();
	void init();

	bool check(double x);

	void print(double x);
};

#endif /* BISEARCH_H_ */
