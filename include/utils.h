#ifndef UTILS_H
#define UTILS_H

#include <vector>
#include <set>
#include <map>
#include <iostream>
#include<fstream>
#include<string>
#include<vector>
#include<cstdlib>
#include<cassert>

namespace TRANSPLAN
{
  class Utils
  {
  public:
    Utils();
    static void print_vector(std::vector<int>& l);

    static void print_set(std::set<int>& l);

    static void print_pair_set(std::set<std::pair<int,int> >& l);

    static void print_int_map(std::map<int,int>& m);

    static bool does_intersect(std::set<int>::iterator begin_set_1, std::set<int>::iterator end_set_1,
			       std::set<int>::iterator begin_set_2, std::set<int>::iterator end_set_2);

    static bool exists(int index, std::vector<int>& vlist);

    static void int_set_to_int_vec(std::set<int>& from, std::vector<int>& to);

    static void int_vec_to_int_set(std::vector<int>& from, std::set<int>& to);

    static void remove_from_vector(int index, std::vector<int>& vlist);
    
    static void trim(std::string& str);

    static void tokenize(std::string& str, std::vector<std::string>& tokens, const std::string& delim);

    static void get_int_interval(std::string& str, std::pair<int,int>& interval);

    //(,,,,)
    static void extract_from_first_bracket(std::string& str, std::vector<std::string>& tokens);

    //{,,,,}
    static void extract_from_second_bracket(std::string& str, std::vector<std::string>& tokens);

    //[,,,,]
    static void extract_from_third_bracket(std::string& str, std::vector<std::string>& tokens);

    ///<,,,,,>
    static void extract_from_angle_bracket(std::string& str, std::vector<std::string>& tokens);

    //Random number generations

    //use this first function to seed the random number generator,
    //call this before any of the other functions
    static void initrand();
    
    //generates a psuedo-random integer between 0 and 32767
    static int randint();
    
    //generates a psuedo-random integer between 0 and max
    static int randint(int max);
   
    //generates a psuedo-random integer between min and max
    static int randint(int min, int max);
    
    //generates a psuedo-random float between 0.0 and 0.999...
    static float randfloat();
    
    //generates a psuedo-random float between 0.0 and max
    static float randfloat(float max);

    //generates a psuedo-random float between min and max
    static float randfloat(float min, float max);
    
    //generates a psuedo-random double between 0.0 and 0.999...
    static double randdouble();
    
    //generates a psuedo-random double between 0.0 and max
    static double randdouble(double max);
    
    //generates a psuedo-random double between min and max
    static double randdouble(double min, double max);
    
    static char* randstr(char* str, unsigned long length);

    static void cartesian_product(int set_lvl,
				  std::vector<std::vector<int> >& sets, 
				  std::vector<std::vector<int> >& prods, 
				  std::vector<int> curr_pr = std::vector<int>());

    static void sort_set_ascending  (std::set<int>& set_to_sort, std::map<int,int>& keys, std::vector<int>& sorted_set);
    static void sort_set_descending (std::set<int>& set_to_sort, std::map<int,int>& keys, std::vector<int>& sorted_set);

    static void swap(std::vector<int>& vec, int from, int to);

    static std::vector<int> set_to_vector(const std::set<int>& from_set);
  };
}
#endif
