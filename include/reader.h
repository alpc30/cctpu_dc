#ifndef READER_H_
#define READER_H_

//#include "constants.h"
//#include "structures.h"
#include "stnu.h"
#include "utils.h"

namespace CDS
{
 // class Utils;
 // class Resource;
 // class Activity;

  class Reader
  {

  public:
	stnu* t;

    Reader(){t=NULL;};
    Reader(stnu* p_stnu){t = p_stnu;/*printf("%d\n",t);*/};

    stnu * getOwner(){/*printf("%d\n",p_stnu);*/return t;};
  //  Transplan * getOwner(){return p_transplan;};

    //Read problem file and create state-variables, actions and transitions
    void read(std::string problem_file, stnu & s);
    void readout(std::string problem_file);
    void readcctp(std::string problem_file);
    void readpos(std::string problem_file);
    void readgrid(std::string problem_file);
    void readppstp(std::string problem_file);
    void readevacplan(std::string problem_file);
    void readcstnu(std::string problem_file);
    void readcst(std::string problem_file);
    void PrintDot(std::string problem_file);

    int getType(std::string x);

    std::string getName(std::string problem_file);

    // void read_preced_matrix(std::string preced_matrix_file);
  };
}
#endif
