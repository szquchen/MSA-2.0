#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <iterator>
#include "molecule_simple.hh"
#include "polynomial.hh"
#include "monomial.hh"

using namespace std;
using namespace molecule;
using namespace MONO;
using namespace POLY;

int main(int argc, char **argv){
  
  if (argc < 3 ) {
    cerr << "please specify the molecule type!"<< endl;
    cerr << "symtax: "<< endl;
    cerr << "./msa [max degree] [molecule configuration] "<< endl;
    return 1;
  }
  
  
  // the first arguments: the maximum degree
  size_t mxd = atoi(argv[1]); 
  
  // check the molecule
  
  // read in the molecule configuration
  size_t ng = argc - 2;
  size_t *mol = new size_t[ng];
  size_t natom = 0;
  ostringstream strBuff;
  strBuff << "MOL_";
  for(size_t i=0; i < ng;i++) {
    mol[i] = atoi(argv[i+2]);
    natom += mol[i];
    strBuff << mol[i] << "_";
  }
  strBuff << mxd;
  
  size_t nbond = natom * (natom - 1)/2;
  
  // print out the molecule configuration
  // cout << "MOLECULE:";
  // copy(mol,mol+ng, ostream_iterator<int>(cout, " ")); 
  // cout << endl;
  
  monomial mono(nbond, mxd);
  monomial pm;			// permuted monomial
  
  set<monomial,  monocmp> bseed; // basis seed, seed for each mono oribit
  set<monomial,  monocmp> orbit; // orbit of a monomial
  set<polynomial,polycmp> basis;
  
  pair<set<monomial,monocmp>::iterator, bool> stat;
  
  // construct the invariant basis 
  permutation pmol(ng,mol);
  do{
    orbit.clear();
    orbit.insert(mono);
    do{
      // permutation for the bond length variable
      vector<size_t> pbond = mole_bond(pmol.seq()); 
      pm = mono;
      pm.permute(pbond);
      orbit.insert(pm);
    }while(pmol++);
    
    stat = bseed.insert(*orbit.begin());
    if(stat.second){
      polynomial poly(orbit);
      basis.insert(poly);
    }
  }while(mono++);
  
  // copy the invariant basis into a vector
  vector<polynomial> vbasis(basis.size());
  set<polynomial, polycmp>::iterator spit;
  size_t i=0;
  for(spit=basis.begin();spit!=basis.end();spit++){
    vbasis[i] = *spit;
    vbasis[i].id(i);
    i++;
  }
  
  // clear the basis in the set structure
  basis.clear();
  
  // ========================================================================
  // now we want to loop over all the polynomials, and find out all those
  // single term polynomials that cannot be factorized.
  
  vector<polynomial>::iterator vpit;
  vector<monomial> prime;
  vector<monomial>::iterator vmit;
  vpit = vbasis.begin(); 
  if(vpit->seed().is_zero()) vpit++;	// skip the first poly if it is zero;
  
  for(;vpit!=vbasis.end();vpit++){
    if(vpit->size()==1){		// single term poly
      monomial seed;
      seed = vpit->seed();
      if(seed.is_prime(prime))	
	prime.push_back(seed);
    }
  }
  
  // ========================================
  // show the prime polynomials
  // ========================================
  /*
  for(vmit=prime.begin();vmit!=prime.end();vmit++) 
  cout << *vmit << endl;
  cout << "--" << endl;
  */
  // ========================================
  // save the entire originial invariant basis
  // ========================================
  string file_bas = strBuff.str() + ".BAS";
  ofstream fbas;
  fbas.open(file_bas.c_str());
  for(vpit=vbasis.begin();vpit!=vbasis.end();vpit++) fbas << *vpit;
  fbas.clear();
  fbas.close();
  
  // factor a high order term into a product and some minus
  string file_foc=strBuff.str() + ".FOC";
  string file_map=strBuff.str() + ".MAP";
  ofstream ffoc, fmap;
  ffoc.open(file_foc.c_str());
  fmap.open(file_map.c_str());
  
  // this is the vector to save the information for every basis, the first
  // element is used to save the basis character:
  // ====================
  // 0: prime  mono
  // 1: second mono
  // 2: prime  poly
  // 3: second poly
  // ========================================
  
  vector< vector<size_t> > fbasis; 
  vector<size_t> cur_basis,tmpv;
  vector<monomial> foc_amonos, monos;
  polynomial pa,pb,ph;
  vpit=vbasis.begin();
  //ffoc << "@"  << setw(3) << vpit->degree() << "|" << *vpit;
  ffoc << *vpit;
  monos = vpit->monos();
  foc_amonos.insert(foc_amonos.end(),monos.begin(),monos.end());
  
  // put the mono into the fbasis;
  cur_basis.push_back(2);	// character, prime mono
  cur_basis.push_back(0);	// 
  fbasis.push_back(cur_basis);
  cur_basis.clear();
  
  vpit++;
  
  for(;vpit!=vbasis.end();vpit++){ // loop over the basis set
    ph = *vpit;
    if(is_in(ph.seed(), prime)) {
      //ffoc << "@" << setw(3) << ph.degree() << "|" << ph;
      ffoc << ph;
      monos = ph.monos();
      // there are supposed to be prime_poly with charater 2, and we
      // save the mono ids
      cur_basis.push_back(2);
      for(size_t i=0;i<monos.size();i++)
	cur_basis.push_back(foc_amonos.size()+i);
      fbasis.push_back(cur_basis);
      cur_basis.clear();
      foc_amonos.insert(foc_amonos.end(),monos.begin(),monos.end());
    }
    else{
      vector<size_t> dd = ph.decomp(bseed, vbasis, prime);
      
      if(check_order(dd)){
	// fmap << " " << setw(3) << vpit->degree() << "|";
	// fmap << setw(4) << dd[0] << "|";
	// fmap << setw(4) << dd[0] << setw(4);
	for(size_t i=0;i<dd.size();i++) fmap << setw(4) << dd[i] << " ";
	fmap << endl;
	vector<size_t>::iterator it;
	it = dd.begin();
	it++;
	cur_basis.push_back(3);
	cur_basis.insert(cur_basis.end(),it,dd.end());
	fbasis.push_back(cur_basis);
	cur_basis.clear();
      }
      else {
	// ffoc << "#" << setw(3) << vpit->degree() << "|" << endl;
	ffoc << ph;
	monos = ph.monos();
	cur_basis.push_back(2);
	for(size_t i=0;i<monos.size();i++)
	  cur_basis.push_back(foc_amonos.size()+i);
	fbasis.push_back(cur_basis);
	cur_basis.clear();
	foc_amonos.insert(foc_amonos.end(),monos.begin(),monos.end());
      }
    }  
  }
  ffoc.close(); 
  fmap.close(); 
  
  // for(vmit=foc_amonos.begin();vmit !=foc_amonos.end();vmit++)
  //   cout << *vmit << endl;
  vector< vector<size_t> > mono_comp;
  // we need to save the composition of every mono with a status indicator
  // 0: prime mono
  // 1: second mono
  
  vmit=foc_amonos.begin();
  
  cur_basis.push_back(0);	// prime mono
  tmpv = vmit->out();		// the decomposition
  cur_basis.insert(cur_basis.end(),tmpv.begin(),tmpv.end());
  mono_comp.push_back(cur_basis);
  cur_basis.clear();
  vmit++;
  for(;vmit!=foc_amonos.end();vmit++){
    // cout << setw(3) << vmit - foc_amonos.begin() << '|'
    // 	 << *vmit << " : ";
    vector<monomial>::iterator it;
    it=foc_amonos.begin();
    it++;
    for(;it!=vmit;it++){
      if(vmit->has_factor(*it)) {
	size_t ida;
	ida = it-foc_amonos.begin();
	monomial mb;
	vmit->pdecomp(*it, mb);
	size_t idb = 0;
	vector<monomial>::iterator ita;
	ita=foc_amonos.begin();
	ita++;
	for(;ita!=vmit;ita++){
	  if(*ita == mb) {
	    idb = ita - foc_amonos.begin();
	    break;
	  }
	}
	if(idb==0) continue;
	// second mono
	cur_basis.push_back(1);	  // second mono
	cur_basis.push_back(ida); // second mono decomp
	cur_basis.push_back(idb); // second mono decomp
	mono_comp.push_back(cur_basis);
	cur_basis.clear();
	
	// cout << *it << "(" << ida << ")" << "+ " 
	//      << mb  << "(" << idb << ")"
	//      << endl;
	
	break;
      }
    }
    if(it==vmit) {
      // vmit is a prime mono
      cur_basis.push_back(0);	// prime mono
      tmpv = vmit->out();	// the decomposition
      cur_basis.insert(cur_basis.end(),tmpv.begin(),tmpv.end());
      mono_comp.push_back(cur_basis);
      cur_basis.clear();
      // cout << endl;
    }
  }
  vector< vector<size_t> >::iterator vvit;
  // show the mono_set configuration
  string file_mono = strBuff.str() + ".MONO";
  ofstream fmono;
  fmono.open(file_mono.c_str());
  for(vvit = mono_comp.begin();vvit!=mono_comp.end();vvit++){
    copy(vvit->begin(), vvit->end(), ostream_iterator<size_t>(fmono, " ")); 
    fmono << endl;
  }
  fmono.close();
  
  // show the fbasis configuration
  string file_poly = strBuff.str() + ".POLY";
  ofstream fpoly;
  fpoly.open(file_poly.c_str());
  for(vvit = fbasis.begin();vvit!=fbasis.end();vvit++){
    copy(vvit->begin(), vvit->end(), ostream_iterator<size_t>(fpoly, " ")); 
    fpoly << endl;
  }
  fpoly.close();

  return 0;
}


