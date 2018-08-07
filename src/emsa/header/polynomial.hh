#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <iostream>
#include <set>
#include "monomial.hh"

namespace POLY{
  using namespace MONO;
  using std::ostream;
  using std::set;
  using std::vector;
  
  class polynomial{
  private:
    size_t        id_;
    monomial      seed_;
    set<monomial,monocmp> orbit_;
  public:
    polynomial();
    polynomial(const monomial & seed, const size_t id=0);
    polynomial(const set<monomial,monocmp> & orbit, const size_t id=0);
    
    void         id(const size_t index){ id_ = index;};
    size_t       id()const{return id_;};
    size_t       size()const{return orbit_.size();};
    size_t       degree()const{return seed_.degree();};
    bool         getid(const vector<polynomial>& basis);
    bool         operator<(const polynomial &rhs)const;
    bool         has(const monomial &mono)const;
    bool         contains(const polynomial &poly)const;
    bool         empty()const{return orbit_.empty();};
    bool         is_01()const;
    bool         is_basic()const;
    bool         is_basic(const vector<monomial> &)const;
    vector<size_t> fact(const monomial &ma, 
			const monomial &mb, 
			const set<monomial,monocmp> & bseed,
			const vector<polynomial>& vb,
			const vector<monomial>& prime)const;
    vector<size_t> decomp(const set<monomial,monocmp> & bseed,
			  const vector<polynomial>& vb,
			  const vector<monomial>& prime)const;
    monomial     seed()const{return seed_;};
    monomial     comfactor()const;
    polynomial & operator=(const polynomial &other);
    polynomial & operator-(const polynomial &other);
    vector<monomial>  operator*(const polynomial &rhs)const;
    vector<monomial>  monos()const; // return all the monos in the poly
    
    friend ostream & operator<<(ostream & os,  const polynomial & poly);
    friend void remove(vector<monomial> &, const polynomial &);
  };
  
  struct polycmp {
    bool operator() (const polynomial& lhs, const polynomial& rhs) const
    {return lhs < rhs;}
  };
  
  
  bool is_in(const monomial &, const vector<monomial> &);
  polynomial get_orbit(const monomial& seed, const vector<polynomial> &vb);
  
  // decomposition an sequence of monomials into disjoint polynomials that are
  // characterized by their seeds. the set of unique seeds is needed as an
  // argument
  vector<monomial> partition(vector<monomial> &, const set<monomial,monocmp>&,
			     const vector<polynomial>&);
  bool check_order(const vector<size_t>&pid);
}
#endif
