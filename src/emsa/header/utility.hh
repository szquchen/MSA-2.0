#ifndef UTILITY_H_
#define UTILITY_H_

#include <string>
namespace Utility{
  using std::string;
  using std::size_t;
  
  string basename(const string &fn){ 
    string bn; 
    size_t np=fn.find_last_of('.'); 
    if(np==string::npos) bn = fn; 
    else bn=fn.substr(0,np); 
    return bn; 
  } 
  string basename(const char *cfn){ 
    string fn=cfn;
    string bn; 
    size_t np=fn.find_last_of('.'); 
    if(np==string::npos) bn = fn; 
    else bn=fn.substr(0,np); 
    return bn; 
  } 
}

#endif
