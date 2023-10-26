/*
 * cstnu.cc
 *
 *  Created on: 30/11/2015
 *      Author: jing
 */

#include "../include/cstnu.h"
#include "reader.h"

namespace CSTNU
{

cstnu::cstnu()
{
	// TODO Auto-generated constructor stub

}

cstnu::cstnu(std::string problem_file, TRANSPLAN::S_PARAM s_info)
{
	stnu(problem_file, s_info);
}

void cstnu::readCSTNU(const std::string & problem_file)
{
	this->file_name = problem_file;
	reader = new TRANSPLAN::Reader(this);
	reader->read(problem_file);
}

cstnu::~cstnu()
{
	// TODO Auto-generated destructor stub
}

} /* namespace TRANSPLAN */
