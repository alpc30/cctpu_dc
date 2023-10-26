/*
 * cstnu.h
 *
 *  Created on: 30/11/2015
 *      Author: jing
 */

#ifndef CSTNU_H_
#define CSTNU_H_

#include <vector>
#include <string>
#include "stnu.h"

namespace CSTNU
{

	class DisVar //Discrete Variables in CSTNU
	{
	public:
		//do disvars activate edges or nodes? (try edges)
		int id;
		std::vector<int> v_dc;
		std::vector<std::vector<int> > vv_id_act_edge;
	};

	typedef class CLINK
	{
		std::vector<int> v_dis_var; //-1=null
		int st, en;
		CLINK(int tst, int ten)
		{
			this->st = tst;
			this->en = ten;
			//TRANSPLAN::LINK (st, en);
			v_dis_var.clear();
		}
//	CLINK(int st, int en) : TRANSPLAN::LINK(st,en)
//	{
//		v_dis_var.clear();
//	}
	} CCTG, CRQM;

	class cstnu: public TRANSPLAN::stnu
	{

	public:

		std::vector<DisVar> v_dis_vars;

		cstnu();
		cstnu(std::string problem_file, TRANSPLAN::S_PARAM info);

		void readCSTNU(const std::string & problem_file);
		virtual ~cstnu();
	};

} /* namespace TRANSPLAN */
#endif /* CSTNU_H_ */
