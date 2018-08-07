#ifndef MONOMIAL_H_
#define MONOMIAL_H_

#include <iostream>
#include <vector>

namespace MONO{
  using std::ostream;
  using std::istream;
  using std::vector;
  
  class monomial{
  private:
    size_t width_; 
    size_t max_degree_;
    // size_t coef_;
    size_t *value_;
    size_t degree_;
    bool overflow_;
    bool full(){return degree_>=max_degree_?true:false;};
  public:
    monomial();
    monomial(const monomial & mono);
    monomial(const size_t n, const size_t d);
    monomial(const size_t d, const vector<size_t> &v);
    monomial(const size_t n, const size_t d, const size_t *v);
    monomial & permute(const size_t *p);
    monomial & permute(const vector<size_t> &p);
    monomial & operator=(const monomial &other);
    monomial operator*(const monomial &rhs);
    vector<size_t> out()const;
    bool is_01()const;
    bool is_basic()const;
    bool is_basic(const vector<monomial> &)const;
    bool is_prime(const vector<monomial> &)const;
    bool is_zero()const;
    bool has_factor(const monomial &)const;
    void decomp(monomial &ma, monomial &mb, size_t)const;
    void pdecomp(const monomial &ma, monomial &mb)const;
    // void inc(){coef_++;};	// increase the coefficients
    // void dec(){coef_--;};	// increase the coefficients
    bool operator++(int);
    bool operator<(const monomial &rhs)const;
    bool operator==(const monomial &rhs)const;
    size_t degree()const{return degree_;};
    size_t entropy()const;
    size_t size()const{return width_;};
    void reset();
    ~monomial();
    
    bool top()const {return overflow_;};
    
    friend ostream & operator<<(ostream & os,  const monomial & mono);
    friend bool operator>>(istream & is, monomial & mono);
  };
  
  class permutation{
  private:
    size_t ng_;
    size_t np_;
    size_t *group_;
    size_t *perm_;
  public:
    permutation();
    permutation(const size_t & np);
    permutation(const size_t & ng, const size_t *group);
    ~permutation();
    permutation & setp(const size_t &np, const size_t *perm);
    permutation & reset();	// default is to reset all the perm
    permutation & reset(const size_t idx); // default is to reset all the perm
    bool operator++(int);
    vector<size_t> seq()const;
    friend ostream & operator<<(ostream & os,  const permutation & perm);
    friend monomial permute(const monomial & mono, const permutation & p);
  };

  // monomial permute(const monomial & mono, const size_t *p);
  struct monocmp {
    bool operator() (const monomial& lhs, const monomial& rhs) const
    {return lhs < rhs;}
  };
  
}
#endif
