/*
 * testgenerator.h
 *
 *  Created on: 16/05/2016
 *      Author: jing
 */

#ifndef TESTGENERATOR_H_
#define TESTGENERATOR_H_

#include "stnu.h"
#include "cds.h"
#include "controllability.h"


namespace CSTNU
{

class TestGenerator
{
	int n_stnu;
	int n_ori_ndc;
	int n_mod_dc;
	int n_lb;
	int n_ub;
public:
	TestGenerator();
	string createTest(CDS::stnu &n);
	void expandOnVariables(CDS::stnu &n);
	void createSTNU(const CDS::stnu &n, std::vector<int> va);
	virtual ~TestGenerator();
};

} /* namespace CSTNU */
#endif /* TESTGENERATOR_H_ */
