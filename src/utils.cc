#include "utils.h"

namespace TRANSPLAN
{
  Utils::Utils()
  {

  }


  bool Utils::exists(int val, std::vector<int>& vlist)
  {
    std::vector<int>::iterator itr = vlist.begin();
    for ( ; itr != vlist.end(); itr++ )
    {
      if ( *itr == val )
	return true;
    }

    return false;
  }


  void Utils::int_set_to_int_vec(std::set<int>& from, std::vector<int>& to)
  {
    
    std::set<int>::iterator itr = from.begin();
    for (; itr != from.end(); ++itr)
      to.push_back(*itr);
  }


  void Utils::int_vec_to_int_set(std::vector<int>& from, std::set<int>& to)
  {    
    for (int i = 0; i < from.size(); ++i)
      to.insert(i);
  }


  void Utils::print_vector(std::vector<int>& l)
  {
    std::cout << "{";
    for( std::vector<int>::iterator itr1 = l.begin(); itr1 != l.end(); itr1++ )
    {
      std::cout << *itr1 << ", ";
    }
    std::cout <<  "}" << std::endl;
  }


  void Utils::print_set(std::set<int>& l)
  {
    std::cout << "{";
    for( std::set<int>::iterator itr1 = l.begin(); itr1 != l.end(); itr1++ )
    {
      std::cout << *itr1 << ", " ;
    }
    std::cout << " }" << std::endl;
  }


  void Utils::print_pair_set(std::set<std::pair<int,int> >& l)
  {
    std::cout << "{";
    for( std::set<std::pair<int,int> >::iterator itr1 = l.begin(); itr1 != l.end(); itr1++ )
    {
      std::cout << "("<<itr1->first <<","<<itr1->second<<")"<< ", " ;
    }
    std::cout << " }" << std::endl;
  }

  void Utils::print_int_map(std::map<int,int>& l)
  {
    std::cout << "{";
    for( std::map<int,int>::iterator itr1 = l.begin(); itr1 != l.end(); itr1++ )
    {
      std::cout << "("<<itr1->first <<","<<itr1->second<<")"<< ", " ;
    }
    std::cout << " }" << std::endl;
  }

  bool Utils::does_intersect(std::set<int>::iterator begin_set_1, std::set<int>::iterator end_set_1,
			     std::set<int>::iterator begin_set_2, std::set<int>::iterator end_set_2)
  {
    while ( begin_set_1 != end_set_1 && begin_set_2 != end_set_2 )
    {
      if ( *begin_set_1 == *begin_set_2 )
	return true;
      else if ( *begin_set_1 < *begin_set_2 )
      {
	begin_set_1++;
      }
      else if ( *begin_set_2 < *begin_set_1 )
      {
	begin_set_2++;
      }
    }
    return false;


    /*
    std::set<int> set_1;
    std::set<int> set_2;
    set_1.insert(2);
    set_1.insert(3);
    set_1.insert(15);
    set_1.insert(22);
    set_1.insert(17);

    set_2.insert(21);
    set_2.insert(25);
    set_2.insert(14);
    set_2.insert(1);
    set_2.insert(17);

    std::set<int>::iterator begin_set_1 = set_1.begin();
    std::set<int>::iterator end_set_1 = set_1.end();
    std::set<int>::iterator begin_set_2 = set_2.begin();
    std::set<int>::iterator end_set_2 = set_2.end();

    while ( begin_set_1 != end_set_1 && begin_set_2 != end_set_2 )
    {
      if ( *begin_set_1 == *begin_set_2 )
      {
	std::cout << "Found intersection = " << *begin_set_1 << std::endl;
	break;
      }
      else if ( *begin_set_1 < *begin_set_2 )
      {
	begin_set_1++;
      }
      else if ( *begin_set_2 < *begin_set_1 )
      {
	begin_set_2++;
      }
    }

    std::exit(0);
    */

  }

  void Utils::remove_from_vector(int index, std::vector<int>& vlist)
  {
    std::vector<int>::iterator itr = vlist.begin();
     for ( ; itr != vlist.end(); itr++ )
     {
       if ( index == *itr )
       {
	 vlist.erase(itr);
	 return;
       }
     }
  }


  void Utils::cartesian_product(int set_lvl, std::vector<std::vector<int> >& sets, 
				std::vector<std::vector<int> >& prods, std::vector<int> curr_pr)
  {
    //base condition
    if (set_lvl == sets.size())
    {
      prods.push_back(curr_pr);
    }
    else
    {
      int ele_lvl = 0;
      while (ele_lvl < sets[set_lvl].size())
      {
	int n_set_lvl = set_lvl+1;
	curr_pr.push_back(sets[set_lvl][ele_lvl]);
	ele_lvl++;
	cartesian_product(n_set_lvl, sets, prods, curr_pr);
	curr_pr.pop_back();
      }
    }
  }

  void Utils::swap(std::vector<int>& vec, int from, int to)
  {
    int temp = vec[to];
    vec[to] = vec[from];
    vec[from] = temp;
  }

  std::vector<int> Utils::set_to_vector(const std::set<int>& from_set)
  {
    std::vector<int> to_vec;
    for (std::set<int>::iterator itr = from_set.begin(); itr != from_set.end(); ++itr)
    {
      to_vec.push_back(*itr);
    }

    to_vec.resize(from_set.size());
    return to_vec;
  }

  ///increasing order of keys
  void Utils::sort_set_ascending (std::set<int>& set_to_sort, std::map<int,int>& keys, std::vector<int>& sorted_set)
  {
    for (std::set<int>::iterator itr = set_to_sort.begin(); itr != set_to_sort.end(); ++itr)
    {
      int curr_ele = *itr;
      int key = keys[curr_ele];

      sorted_set.push_back(curr_ele);

      int inserted_pos = sorted_set.size()-1;

      int prev_pos = inserted_pos-1; 

      while (prev_pos >= 0)
      {
	int prev_key = keys[sorted_set[prev_pos]];

	if (prev_key > key)
	{
	  int temp = sorted_set[prev_pos];
	  sorted_set[prev_pos] = curr_ele;
	  sorted_set[inserted_pos] = temp;
	  inserted_pos = prev_pos;
	}

	prev_pos--;
      }
    }
  }

  //decreasing order of keys
  void Utils::sort_set_descending (std::set<int>& set_to_sort, std::map<int,int>& keys, std::vector<int>& sorted_set)
  {
    for (std::set<int>::iterator itr = set_to_sort.begin(); itr != set_to_sort.end(); ++itr)
    {
      int curr_ele = *itr;
      int key = keys[curr_ele];

      sorted_set.push_back(curr_ele);

      int inserted_pos = sorted_set.size()-1;

      int prev_pos = inserted_pos-1; 

      while (prev_pos >= 0) 
      {
	int prev_key = keys[sorted_set[prev_pos]];

	if (prev_key < key)
	{
	  int temp = sorted_set[prev_pos];
	  sorted_set[prev_pos] = curr_ele;
	  sorted_set[inserted_pos] = temp;
	  inserted_pos = prev_pos;
	}

	prev_pos--;
      }
    }
  }

  //////////////////////////////////////////////////////////
  //        String Extranction Methods            //////////////
  /////////////////////////////////////////////////////////

  void Utils::trim(std::string& str)
  {
    //str = "    THIS     ";
    //std::cout << "We are triming " << str << " ...." << std::endl;
    std::string::size_type start_pos = str.find_first_not_of(" ", 0);
    str = str.substr(start_pos, str.size());
    std::string::size_type last_pos = str.find_first_of(" ", 0);
    str = str.substr(0, last_pos);
    //std::cout << "We have trimed " << str << " ...." << std::endl;
  }

  void Utils::tokenize(std::string& str, std::vector<std::string>& tokens, const std::string& delim)
  { 
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delim, 0);
   
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delim, lastPos);

    while (std::string::npos != pos || std::string::npos != lastPos)
    {
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
      // Skip delimiters.  Note the "not_of"
      lastPos = str.find_first_not_of(delim, pos);
      // Find next "non-delimiter"
      pos = str.find_first_of(delim, lastPos);
    }

    for (int i = 0; i < tokens.size(); ++i )
    {
      Utils::trim(tokens[i]);
    }
  }

  void Utils::get_int_interval(std::string& str, std::pair<int,int>& interval)
  {
    //std::cout << "Interval ... " << str << std::endl;
    std::vector<std::string> tokens;
    tokens.reserve(2);
    Utils::tokenize(str, tokens, std::string(" [,]"));
    assert(tokens.size() <= 2);
    if (tokens.size() > 0)
      interval.first = atoi(tokens[0].c_str());
    if (tokens.size() > 1)
      interval.second = atoi(tokens[1].c_str());
  }

  //(,,,,)
  void Utils::extract_from_first_bracket(std::string& str, std::vector<std::string>& tokens)
  {
    Utils::tokenize(str, tokens, std::string(" (,)"));
    for (int i = 0; i < tokens.size(); ++i )
    {
      Utils::trim(tokens[i]);
    }
  }

  //{,,,,}
  void Utils::extract_from_second_bracket(std::string& str, std::vector<std::string>& tokens)
  {
    Utils::tokenize(str, tokens, std::string(" {,}"));
    //assert(tokens.size() > 0);
    for (int i = 0; i < tokens.size(); ++i )
    {
      Utils::trim(tokens[i]);
    }
  }


  //[,,,,]
  void Utils::extract_from_third_bracket(std::string& str, std::vector<std::string>& tokens)
  {
    Utils::tokenize(str, tokens, std::string(" [,]"));
    //assert(tokens.size() > 0);
    for (int i = 0; i < tokens.size(); ++i )
    {
      Utils::trim(tokens[i]);
    }
  }

  ///<,,,,,>
  void Utils::extract_from_angle_bracket(std::string& str, std::vector<std::string>& tokens)
  {
    Utils::tokenize(str, tokens, std::string(" <,>"));
    //assert(tokens.size() > 0);
    for (int i = 0; i < tokens.size(); ++i )
    {
      Utils::trim(tokens[i]);
    }
  }

  
  //******************************************************************************************************
  //*******************************  Random number generations   *****************************************
  //******************************************************************************************************
  
  //use this first function to seed the random number generator,
  //call this before any of the other functions
  void Utils::initrand()
  {
    srand((unsigned)(time(0)));
  }
  
  //generates a psuedo-random integer between 0 and 32767
  int Utils::randint()
  {
    return rand();
  }
  
  //generates a psuedo-random integer between 0 and max
  int Utils::randint(int max)
  {
    return int(max*rand()/(RAND_MAX+1.0));
  }
  
  //generates a psuedo-random integer between min and max
  int Utils::randint(int min, int max)
  {
    if (min>max)
    {
      return max+int((min-max+1)*rand()/(RAND_MAX+1.0));
    }
    else
    {
      return min+int((max-min+1)*rand()/(RAND_MAX+1.0));
    }
  }
  
  //generates a psuedo-random float between 0.0 and 0.999...
  float Utils::randfloat()
  {
    return rand()/(float(RAND_MAX)+1);
  }
  
  //generates a psuedo-random float between 0.0 and max
  float Utils::randfloat(float max)
  {
    return randfloat()*max;
  }
  
  //generates a psuedo-random float between min and max
  float Utils::randfloat(float min, float max)
  {
    if (min>max)
    {
      return randfloat()*(min-max)+max;    
    }
    else
    {
      return randfloat()*(max-min)+min;
    }    
  }
  
  //generates a psuedo-random double between 0.0 and 0.999...
  double Utils::randdouble()
  {
    return rand()/(double(RAND_MAX)+1);
  }
  
  //generates a psuedo-random double between 0.0 and max
  double Utils::randdouble(double max)
  {
    return randdouble()*max;
  }
  
  //generates a psuedo-random double between min and max
  double Utils::randdouble(double min, double max)
  {
    if (min>max)
    {
      return randdouble()*(min-max)+max;    
    }
    else
    {
      return randdouble()*(max-min)+min;
    }
  }
  
  char* Utils::randstr(char* str, unsigned long length)
  {
    //make we were passed a valid pointer
    if (!str)
    {
      return 0;
    }
    
    //make sure the string is supposed to contain something
    if (!length)
    {
      return 0;
    }
    
    //put random characters into the string, give both
    //upper and lower case numbers an equal chance at
    //being used
    unsigned long x;
    for (x=0;x<length-1;x++)
    {
      if (!randint(1))
      {
	str[x]=(char)randint(65,90);
      }
      else
      {
	str[x]=(char)randint(97,122);
      }
    }
    
    //null terminate the string
    str[x]=0;
    
    return str;
  }


}
