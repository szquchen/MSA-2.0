#ifndef MOLECULE_H_
#define MOLECULE_H_

#include <algorithm>
#include <vector>
namespace molecule{
  using std::vector;
  using std::min;
  using std::max;
  
  vector<size_t> mole_bond(const vector<size_t> &mole){
    // given a molecule (an arrangement of atoms in a molecule), construct the
    // bond length variable vector
    size_t n = mole.size();
    size_t m = n*(n-1)/2;
  
    vector<size_t> bond(m);
  
    size_t k = 0;
    for(size_t i=0;i<n-1;i++){
      for(size_t j=i+1;j<n;j++){
	size_t a = min(mole[i], mole[j]);
	size_t b = max(mole[i], mole[j]);
	bond[k] = n*a - (a+1)*a/2 + b - a - 1;
	k++;
      }
    }
    return bond;
  }
}
#endif
