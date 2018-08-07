#include <iomanip>
#include <numeric>
#include <algorithm>
#include "monomial.hh"
#include "polynomial.hh"

namespace POLY{
  polynomial::polynomial(){
    id_   = 0;
    seed_ = monomial();
    orbit_.insert(seed_);
  }
  
  polynomial::polynomial(const monomial & seed, const size_t id){
    id_   = id;
    seed_ = seed;
    orbit_.insert(seed_);
  }
  
  polynomial::polynomial(const set<monomial,monocmp> & orbit, const size_t id){
    id_    = id;
    seed_  = *(orbit.begin());
    orbit_ = orbit;
  }
  
  monomial polynomial::comfactor()const{
    // get the common factor among all the monomials in the polynomial
    // if the poly has only one mono, then the factor is the monomial
    // less one order
    set<monomial,monocmp>::iterator  it;
    size_t n = seed_.size();
    
    std::vector<size_t> vc(n);
    
    it = orbit_.begin();
    vc = it->out();
    if(size()==1 && seed_.entropy() == 1){
      for(size_t i=0;i<n;i++) {
  	if(vc[i] > 1){
  	  vc[i]-- ;
  	  break;
  	}else vc[i] = 0;
      }
    }else{
      it++;
      while(it!=orbit_.end()){
  	std::vector<size_t> mono = it->out();
  	for(size_t i=0;i<n;i++) vc[i] = std::min(vc[i], mono[i]);
  	it++;
      }
    }
    monomial comm(0,vc);
    return comm;
  }
  
  bool polynomial::operator<(const polynomial &rhs)const{
    // lhs < rhs iff the max(lhs) < max(rhs)
    // compare only the last monomial.
    set<monomial,monocmp>::iterator pa,pb,pt;
    pa = this->orbit_.end(); pa--;
    pb =   rhs.orbit_.end(); pb--;
    return *pa < *pb;
  }
  
  bool polynomial::contains(const polynomial &poly)const{
    monomial mono;
    mono = poly.seed();
    set<monomial,monocmp>::iterator pt;
    for(pt=orbit_.begin();pt!=orbit_.end();pt++){
      if(pt->has_factor(mono)) return true;
    }
    return false;
  }
  
  bool polynomial::has(const monomial &mono)const{
    set<monomial,monocmp>::iterator pit;
    for(pit = orbit_.begin();pit != orbit_.end();pit++){
      if(*pit == mono) return true;
    }
    return false;
  };
  
  bool polynomial::is_01()const{    
    return seed_.is_01();
  };
  bool polynomial::is_basic()const{    
    return seed_.is_basic();
  };
  bool polynomial::is_basic(const vector<monomial> & prime)const{    
    return seed_.is_basic(prime);
  };
  
  bool polynomial::getid(const vector<polynomial>& basis){
    vector<polynomial>::const_iterator pit;
    
    for(pit=basis.begin();pit!=basis.end();pit++) 
      if(seed_ == (*pit).seed_){
	id_ = (*pit).id_;
	return true;
      }
    return false;
  }
  
  vector<monomial>  polynomial::monos()const{
    // get all the monos from a poly
    vector<monomial> vmono;
    set<monomial,monocmp>::iterator mit ;
    for(mit=orbit_.begin();mit!=orbit_.end();mit++)
      vmono.push_back(*mit);
    return vmono;
  }
  polynomial & polynomial::operator=(const polynomial &poly){
    if (this != & poly) {
      id_    = poly.id_;
      seed_  = poly.seed_;
      orbit_ = poly.orbit_;
    }
    return *this;
  }
  
  polynomial & polynomial::operator-(const polynomial &poly){
    if (this != & poly) {
      set<monomial,monocmp>::iterator mit ;
      for(mit=poly.orbit_.begin();mit!=poly.orbit_.end();mit++){
	if(this->has(*mit)) orbit_.erase(*mit);
      }
    }
    return *this;
  };
  vector<monomial> polynomial::operator*(const polynomial &poly)const{
    vector<monomial> prod;
    
    set<monomial,monocmp>::iterator pa,pb;
    for(pa = orbit_.begin();pa!=orbit_.end();pa++){
      for(pb = poly.orbit_.begin();pb!=poly.orbit_.end();pb++){
	monomial ma = *pa;
	monomial mb = *pb;
	monomial mp = ma*mb;
	prod.push_back(mp);
      }
    }
    return prod;
  }
  
  vector<size_t> polynomial::fact(const monomial &ma, 
				  const monomial &mb, 
				  const set<monomial,monocmp> & bseed,
				  const vector<polynomial>& vb,
				  const vector<monomial>& prime)const{
    
    // we want to factor orb(m) = orb(m1)*orb(m2) - orb(m3)...  in this
    // function, we suppose we know m, m1 and m2. all we want to know is m3,
    // etc. the result to return is a vector of polynomials. the first one in
    // this vector is orb(ma), and the second one is orb(mb), then the left
    // ones as orb(m3) ...

    // actually, I do not need return a vector of polynomials, their ID is
    // enough

    vector<size_t> pid;
    pid.push_back(id_);		// save the current poly id
    polynomial pa,pb;
    
    pa = get_orbit(ma, vb); pid.push_back(pa.id());
    pb = get_orbit(mb, vb); pid.push_back(pb.id());
    
    vector<monomial> prod = pa*pb;
    remove(prod,*this);
    vector<monomial> pseeds = partition(prod, bseed, vb);
    vector<monomial>::iterator vmit;
    for(vmit=pseeds.begin();vmit!=pseeds.end();vmit++) {
      pa = get_orbit(*vmit,vb);
      pid.push_back(pa.id());
    }
    return pid;
  }
  
  vector<size_t> polynomial::decomp(const set<monomial,monocmp> & bseed,
				    const vector<polynomial>& vb,
				    const vector<monomial>& prime)const{
    // we want to factor orb(m) = orb(m1)*orb(m2) - orb(m3)...
    // here we only know m, and we want to try several difference m1 and m2 until
    // we get some successful decomposition. if all the possible factorizations
    // were tried, and all failed, then return failed.
    // ZHEN
    // this function is used to decompose a high order poly into the product
    // of two lower order poly less some other poly
    // the return values would be:
    // 1. the first poly: mono + orbit
    // 2. the second poly: mono + orbit
    // 3. a vector of left poly
    // 4. a stat indicator to show whether the decomposition is successful or not
    //    note that of all the left polys are rank lower than the current poly,
    //    then, the decomposition is successful, otherwize, not
    
    vector<size_t> pid;
    monomial ma, mb;
    
    // check if the current poly contains a prime factor
    vector<monomial>::const_iterator vmit;
    for(vmit=prime.begin();vmit!=prime.end();vmit++)
      if(seed_.has_factor(*vmit)){
	ma = *vmit;
	seed_.pdecomp(ma,mb);
	pid = fact(ma,mb,bseed, vb, prime);
	return pid;
      }
    
    
    vector< vector<size_t> > sfact;
    vector<polynomial>::const_iterator it;
    it=vb.begin(); it++;
    for(;it->degree() <= degree()/2+1;it++){	
      if(this->contains(*it)){ 
	set<monomial,monocmp>::const_iterator sit;
	for(sit = it->orbit_.begin();sit!=it->orbit_.end();sit++){
	  if(seed_.has_factor(*sit)) {ma = *sit;break;}
	}
	seed_.pdecomp(ma,mb);
	pid = fact(ma,mb,bseed, vb, prime);
	sfact.push_back(pid);
	//if(check_order(pid)) return pid;
      }
    }
    
    pid = sfact[0];
    bool cstat = check_order(pid);
    for(size_t i=1; i<sfact.size(); i++)
      if(cstat){
	if(sfact[i].size() < pid.size() && check_order(sfact[i]))
	  pid = sfact[i];
      }else{
	if(check_order(sfact[i])){
	  pid = sfact[i];
	  cstat = true;
	}else
	  if(sfact[i].size() < pid.size()) pid = sfact[i];
      }
    
    /*
    for(;it->id()< id_;it++){	
      if(it->is_01()){
    	if(this->contains(*it)){ // *it may be a factor: m1
	  // *it = orb(m1);
    	  set<monomial,monocmp>::const_iterator sit;
    	  for(sit = it->orbit_.begin();sit!=it->orbit_.end();sit++){
    	    if(seed_.has_factor(*sit)) ma = *sit;
    	  }
    	  seed_.pdecomp(ma,mb);
    	  pid = fact(ma,mb,bseed, vb, prime);
    	  if(check_order(pid)) return pid;
    	}
      }
    }
    */
    
    // seed_.decomp(ma,mb,0);
    // pid = fact(ma,mb,bseed, vb, prime);
    
    // for(size_t s=0;s<5;s++){
    //   seed_.decomp(ma,mb,s);
    //   pid = fact(ma,mb,bseed, vb, prime);
    //   if(check_order(pid)) return pid;
    // }
    return pid;
  }
  // ============================================================
  void remove(vector<monomial> & vp, const polynomial &ip){
    set<monomial,monocmp>::const_iterator pa;
    vector<monomial>::iterator vit;
    
    for(pa = ip.orbit_.begin();pa!=ip.orbit_.end();pa++)  // poly
      for(vit=vp.begin();vit!=vp.end();vit++)		   // mono seq
	if(*vit == *pa) {vp.erase(vit); break;}
  };
  
  ostream & operator<<(ostream &os, const polynomial & poly) {
    using std::setw;
    using std::endl;
    
    size_t i=1;
    set<monomial,monocmp>::iterator mit;
    for(mit=poly.orbit_.begin();mit!=poly.orbit_.end();mit++,i++){
      os << setw(4) << poly.id_      << " "
	 << setw(2) << poly.degree() << " "
         << setw(2) << poly.size()   << " "
         << setw(2) << i << " : "
         << *mit    << endl;
    }
    
    // os << setw(4) << poly.id_ << "|"
    //    << "   "
    //    << poly.seed_; 
    // set<monomial,monocmp>::iterator mit = poly.orbit_.begin() ;
    // mit++;
    // for(;mit!=poly.orbit_.end();mit++){
    //   os << "+" << *mit  ;
    // }
    
    return os;
  }
  
  vector<monomial>  partition(vector<monomial> & vp, 
			      const set<monomial,monocmp>& seeds,
			      const vector<polynomial>& vb){
    vector<monomial> pseeds;
    set<monomial,monocmp>::const_iterator sit;
    for(sit=seeds.begin();sit!=seeds.end();sit++){
      if(is_in(*sit, vp)){
	pseeds.push_back(*sit);
	polynomial poly = get_orbit(*sit,vb);
	remove(vp,poly);
      };
    }
    
    vector<monomial> dseeds(pseeds);
    vector<monomial>::iterator it;
    
    while(vp.size() > 0){
      for(it=dseeds.begin();it!=dseeds.end();it++){
    	if(is_in(*it, vp)){
    	  pseeds.push_back(*it);
    	  polynomial poly = get_orbit(*it,vb);
    	  remove(vp,poly);
    	};
      }
    }
    return pseeds;
  }

  polynomial get_orbit(const monomial& seed, const vector<polynomial> &vb){
    vector<polynomial>::const_iterator vit;
    for(vit=vb.begin();vit!=vb.end();vit++){
      if((*vit).has(seed)) return *vit;;
    }
    return *vit;
  }
  bool is_in(const monomial & m, const vector<monomial> & vm){
    vector<monomial>::const_iterator vit;
    for(vit = vm.begin();vit!=vm.end();vit++) if(m == *vit)  return true;
    return false;
  };
  bool check_order(const vector<size_t>&pid){
    using std::min_element;
    size_t n = pid.size();
    if(n==3){ 
      if(*min_element(pid.begin(),pid.begin()+3)==0)  return false;
      else return true;
    }
    else{
      for(size_t i=3;i<n;i++) if(pid[i] > pid[0]) return false;
      return true;
    }
  };
}
  
