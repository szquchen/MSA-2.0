#include <iomanip>
#include <numeric>
#include <algorithm>
#include "monomial.hh"

namespace MONO{
  monomial::~monomial(){
    delete [] value_;
  }
  monomial::monomial(){
    degree_ = 0;
    width_ = 1;
    // coef_ = 1;
    max_degree_ = 0;
    value_ = new size_t[1];
    overflow_ = false;
  }

  monomial::monomial(const monomial & mono){
    // if(width_ > 0) { delete [] value_;};
    degree_ = mono.degree_;
    width_  = mono.width_;
    // coef_   = mono.coef_;
    max_degree_ = mono.max_degree_;
    overflow_ = mono.overflow_;
    value_ = new size_t[width_];
    for(size_t i=0; i< width_; i++){
      value_[i] = mono.value_[i];
    }
  }
  
  monomial & monomial::operator=(const monomial &mono){
    if (this != & mono) {
      degree_ = mono.degree_;
      width_  = mono.width_;
      // coef_ = mono.coef_;
      max_degree_ = mono.max_degree_;
      overflow_ = mono.overflow_;
      delete [] value_;
      value_ = new size_t[width_];
      std::copy(mono.value_, mono.value_ + mono.width_, value_);
    }
    return *this;
  }
  vector<size_t> monomial::out()const{
    using std::copy;
    vector<size_t> var(width_);
    copy(value_,value_+width_,var.begin());
    return var;
  }
  
  monomial::monomial(const size_t n, const size_t d){
    degree_ = 0;
    width_ = n;
    // coef_ = 1;
    max_degree_ = d;
    overflow_ = false;
    value_ = new size_t[n];
  }
  monomial::monomial(const size_t d, const vector<size_t> &v){
    using std::max;
    using std::accumulate;
    size_t total_degree = accumulate(v.begin(),v.end(),0);
    
    degree_ = total_degree;
    width_ = v.size();
    // coef_ = 1;
    max_degree_ =  max(total_degree, d);
    if(full() && v[0] == max_degree_){
      overflow_ = true;
    }else{
      overflow_ = false;
    }
    value_ = new size_t[width_];
    for(size_t i=0;i<width_;i++){
      value_[i] = v[i];
    }
  }
  
  monomial::monomial(const size_t n, const size_t d, const size_t *v){
    using std::max;
    using std::accumulate;
    size_t total_degree = accumulate(v,v+n,0);
    
    degree_ = total_degree;
    width_ = n;
    // coef_ = 1;
    max_degree_ =  max(total_degree, d);
    if(full() && v[0] == max_degree_){
      overflow_ = true;
    }else{
      overflow_ = false;
    }
    value_ = new size_t[width_];
    for(size_t i=0;i<width_;i++){
      value_[i] = v[i];
    }
  }
  
  bool monomial::is_01()const{
    using std::max_element;
    if(*max_element(value_,value_+width_)==1) return true;
    else return false;
  }
  
  bool monomial::is_basic()const{
    using std::max_element;
    if(*max_element(value_,value_+width_)==1 && entropy() ==1) return true;
    else return false;
  }
  
  
  bool monomial::is_basic(const vector<monomial> & prime)const{
    using std::max_element;
    if(*max_element(value_,value_+width_)>1) return false;

    vector<monomial>::const_iterator vmit;
    for(vmit=prime.begin();vmit!=prime.end();vmit++)
      if(this->has_factor(*vmit)) return false;
    return true;
  }
  
  bool monomial::is_prime(const vector<monomial> & prime)const{
    using std::max_element;
    if(*max_element(value_,value_+width_)>1) return false;
    
    vector<monomial>::const_iterator vmit;
    for(vmit=prime.begin();vmit!=prime.end();vmit++)
      if(this->has_factor(*vmit)) return false;
    return true;
  }
  
  bool monomial::has_factor(const monomial &mo)const{
    for(size_t i=0;i<width_;i++){
      if(value_[i] < mo.value_[i]) return false;
    }
    return true;
  }
  
  bool monomial::is_zero()const{
    using std::max_element;
    if(*max_element(value_,value_+width_)==0) return true;
    else return false;
  }
  
  monomial monomial::operator*(const monomial &rhs){
    monomial prod(rhs); 
    prod.degree_ += degree_;
    for(size_t i=0;i<width_;i++){
      prod.value_[i] += value_[i];
    }
    return prod;
  }
  
  void monomial::decomp(monomial &ma, monomial &mb, size_t sw)const{
    ma = *this;
    mb = *this;
    bool up = true;
    // this is one approach
    switch(sw){
    case 3:
      // std::cout << "case 0" << std::endl;
      for(size_t i=0;i<width_;i++){
	if(this->value_[i] > 1) {
	  mb.value_[i] = 1;
	  ma.value_[i] = value_[i] - 1;
	} else {
	  mb.value_[i] = 0;
	  ma.value_[i] = value_[i];
	}
      }
      break;
    case 2:
      ma.reset();
      for(size_t i=0;i<width_;i++){
	if(this->value_[i] > 0){
	  ma.value_[i] = 1;
	  mb.value_[i] = value_[i]-1;
	  break;
	}
      }
      break;
    case 1:
      mb.reset();
      for(size_t i=0;i<width_;i++){
	if(this->value_[i] > 0){
	  ma.value_[i] = value_[i]-1;
	  mb.value_[i] = 1;
	  break;
	}
      }
      break;
    case 0:
      // this is called smart decomposition, the idea is to seperate the
      // non-zero elements "equally" into two monomials
      ma.reset();
      mb.reset();
      
      for(size_t i=0;i<width_;i++){
	if(value_[i] == 0) continue;
	else if(value_[i] == 1){
	  if(up) {
	    ma.value_[i] = 1;
	    up = false;
	  }else{
	    mb.value_[i] = 1;
	    up = true;
	  }
	}else{
	  // ma.value_[i] = value_[i]/2;
	  // mb.value_[i] = value_[i] - ma.value_[i];
	  ma.value_[i] = value_[i] - 1;
	  mb.value_[i] = 1;
	}
      }
      break;
    default:
      ma.reset();
      mb.reset();
      
      for(size_t i=0;i<width_;i++){
	if(value_[i] == 0) continue;
	else if(value_[i] == 1){
	  if(up) {
	    ma.value_[i] = 1;
	    up = false;
	  }else{
	    mb.value_[i] = 1;
	    up = true;
	  }
	}else{
	  // ma.value_[i] = value_[i]/2;
	  // mb.value_[i] = value_[i] - ma.value_[i];
	  ma.value_[i] = value_[i] - 1;
	  mb.value_[i] = 1;
	}
      }
    }
  }
  
  void monomial::pdecomp(const monomial &ma, monomial &mb)const{
    // here we know that ma is already a factor of the current monomial
    mb = *this;
    for(size_t i=0;i<width_;i++) mb.value_[i] = value_[i] - ma.value_[i];
  }
  
  
  monomial & monomial::permute(const size_t *p){
    // here we suppose the permutation p is a feasible permutation
    size_t *tv = new size_t[width_];
    for(size_t i=0;i<width_;i++){
      tv[i] = value_[p[i]];
    }
    for(size_t i=0;i<width_;i++){
      value_[i] = tv[i];
    }
    delete [] tv;
    return *this;
  }

  monomial & monomial::permute(const vector<size_t> &p){
    // here we suppose the permutation p is a feasible permutation
    size_t *tv = new size_t[width_];
    for(size_t i=0;i<width_;i++){
      tv[i] = value_[p[i]];
    }
    for(size_t i=0;i<width_;i++){
      value_[i] = tv[i];
    }
    delete [] tv;
    return *this;
  }
  
  
  
  bool monomial::operator++(int){ 
    if(overflow_) return false;
    size_t i = 0;
    while(i < width_){
      if(full()){
	degree_ -= value_[i];
	value_[i++]=0;
      }else{
	value_[i]++;
	degree_++;
	break;
      }
    };
    if(i==width_) {
      overflow_ = true; 
      return false;
    }
    
    return true;
  }
  
  void monomial::reset(){
    degree_ = 0;
    for(size_t i=0;i<width_;i++) value_[i] = 0;
    overflow_ = false;
  }

  bool monomial::operator<(const monomial &rhs)const{
    using std::lexicographical_compare;
    
    if(this->degree() < rhs.degree()) return true;
    else if(this->degree() > rhs.degree()) return false;
    else{
      if(this->entropy() > rhs.entropy()) return true;
      else if(this->entropy() < rhs.entropy()) return false;
      else return lexicographical_compare(value_,value_+width_,
					  rhs.value_, rhs.value_+rhs.width_);
    }
  }
  
  bool monomial::operator==(const monomial &rhs)const{
    using std::equal;
    return equal(value_,value_+width_, rhs.value_);
  }
  
  size_t monomial::entropy()const{
    size_t e = 0;
    for(size_t i=0;i<width_;i++) if (value_[i] > 0) e++;
    return e;
  };

  // ============================================================
  permutation::~permutation(){
    delete [] group_;
    delete [] perm_;
  }
  permutation::permutation(){
    ng_ = 1;
    np_ = 1;
    group_ = new size_t[1];
    perm_ = new size_t[1];
  }
  permutation::permutation(const size_t &np){
    ng_ = 1;
    np_ = np;
    group_ = new size_t[1];
    perm_  = new size_t[np_];
    for(size_t i =0;i<np;i++){perm_[i] = i;};
  }
  
  permutation::permutation(const size_t & ng, const size_t *group){
    using std::accumulate;
    ng_ = ng;
    np_ = accumulate(group,group+ng,0);
    group_ = new size_t[ng_];
    perm_  = new size_t[np_];
    for(size_t i=0; i < ng_; i++){group_[i] = group[i];};
    for(size_t i=0; i < np_; i++){perm_[i] = i;};
  };
  
  permutation & permutation::setp(const size_t &np, const size_t *perm){
    if(np != np_) return *this;
    for(size_t i=0; i<np;i++){perm_[i] = perm[i];};
    return *this;
  };
  permutation & permutation::reset(){
    using std::sort;
    sort(perm_,perm_+ np_);
    return *this;
  };
  permutation & permutation::reset(const size_t idx){
    using std::sort;
    using std::accumulate;
    sort(perm_+accumulate(group_,group_+idx,0), 
	 perm_+accumulate(group_,group_+idx+1,0));
    return *this;
  };
  
  bool permutation::operator++(int){ 
    using std::next_permutation;
    using std::accumulate;
    size_t i = 0;
    while(i < ng_){
      size_t a = accumulate(group_,group_+i,0);
      size_t z = accumulate(group_,group_+i+1,0);
      if(next_permutation(perm_+a,perm_+z)){
	i = 0;
	break;
      }else{
	reset(i);
	i++;
      }
    }
    if(i == ng_) {return false;}
    else{return true;};
  };
  vector<size_t> permutation::seq()const{
    vector<size_t> perm(perm_,perm_+np_);
    return perm;
  }
  // ============================================================
  ostream & operator<<(ostream &os, const monomial & mono) {
    // using std::setw;
    // using std::right;
    // os << right;
    for(size_t i=0; i < mono.width_; i++){
      // os << setw(3) << mono.value_[i];
      //os << " " << mono.value_[i];
      os << mono.value_[i] << " ";
    }
    return os;
  }
  
  
  bool operator>>(istream & is, monomial & mono){
    // here we suppose that mono width has been set already
    if(is.eof()) return false;
    if(is >> mono.value_[0]){;
      mono.degree_ = mono.value_[0];
      for(size_t i=1;i<mono.width_;i++){
	is >> mono.value_[i];
	mono.degree_+= mono.value_[i];
      }
      mono.overflow_ = false;
      return true;
    }else return false;
    
  }

  ostream & operator<<(ostream &os, const permutation & perm) {
    // using std::setw;
    // using std::right;
    // os << right;
    for(size_t i=0; i < perm.np_; i++){
      // os << setw(3) << perm.perm_[i];
      os << " " << perm.perm_[i];
    }
    return os;
  }
  // this one is a friend of permutation
  monomial permute(const monomial & mono, const permutation & p){
    monomial np(mono);
    np.permute(p.perm_);
    return np;
  }
  // ============================================================
  monomial permute(const monomial & mono, const size_t *p){
    monomial np(mono);
    np.permute(p);
    return np;
  }
  
}
  
