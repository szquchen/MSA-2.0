#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "poly_basis.hh"

namespace poly_basis{
  using std::vector;
  using std::string;

  size_t read_coef(const char * file_coef, vector<double> & vc){
    using std::ifstream;
    
    vc.clear();
    
    ifstream fin;
    fin.open(file_coef);
    if(!fin.good()) {
      std::cerr << "FILE READING ERROR: " << file_coef << std::endl;
      exit(1);
    }
    double coef;
    while(fin >> coef) vc.push_back(coef);
    fin.close();
    return vc.size();
    
  }
  void bdv_trans(double & x){
    x = pow(1.0-exp(-0.2*x),2);
  }
  size_t read_mono(const char *file_mono, 
		 vector< vector<size_t> > & mono_comp){
    using std::ifstream;
    using std::istringstream;
    
    mono_comp.clear();
    
    ifstream fin;
    string linBuff;
    
    fin.open(file_mono);
    if(!fin.good()) {
      std::cerr << "FILE READING ERROR: " << file_mono << std::endl;
      exit(1);
    }
    while(getline(fin,linBuff)){
      istringstream issBuff;
      issBuff.str(linBuff);
      vector<size_t> cmp;
      size_t elt;
      while(issBuff >> elt) cmp.push_back(elt);
      mono_comp.push_back(cmp);
    }
    fin.close(); fin.clear();
    return mono_comp.size();
  }
  
  size_t read_poly(const char *file_poly, 
		 vector< vector<size_t> > & poly_comp){
    using std::ifstream;
    using std::istringstream;
    
    poly_comp.clear();
    
    ifstream fin;
    string linBuff;
    
    fin.open(file_poly);
    if(!fin.good()) {
      std::cerr << "FILE READING ERROR: " << file_poly << std::endl;
      exit(1);
    }
    while(getline(fin,linBuff)){
      istringstream issBuff;
      issBuff.str(linBuff);
      vector<size_t> cmp;
      size_t elt;
      while(issBuff >> elt) cmp.push_back(elt);
      poly_comp.push_back(cmp);
    }
    fin.close(); fin.clear();
    return poly_comp.size();
  }
  
  
  void mono_evaluate(vector<double> & vmono, 
		     const vector< vector<size_t> > & mono_comp,
		     const vector<double> & vx){
    vector< vector<size_t> >::const_iterator vvit;
    size_t k=0;
    for(vvit = mono_comp.begin();vvit!=mono_comp.end();vvit++){
      vector<size_t> comp = *vvit;
      vector<size_t>::iterator it;
      double value=1.0;
      it = comp.begin();
      size_t stat = *it;
      it++;
      if(stat == 0){		// prime mono
	for(;it!=comp.end();it++) {
	  if(*it == 0) continue;
	  else if(*it==1) value *= vx[it-comp.begin()-1];
	  else value *= pow(vx[it-comp.begin()], *it);
	}
	//vmono.push_back(value);
	vmono[k++] = value;
      }else{			// second mono
	value = vmono[comp[1]]*vmono[comp[2]];
	// vmono.push_back(value);
	vmono[k++] = value;
      }
    }
  }
  void poly_evaluate(vector<double> & vpoly, 
		     const vector< vector<size_t> > & poly_comp,
		     const vector<double> & vmono){
    vector< vector<size_t> >::const_iterator vvit;
    size_t k=0;
    for(vvit = poly_comp.begin();vvit!=poly_comp.end();vvit++){
      vector<size_t> comp = *vvit;
      vector<size_t>::iterator it;
      double value=0;
      it = comp.begin();
      size_t stat = *it;
      it++;
      if(stat == 2){		// prime poly
	for(;it!=comp.end();it++) 
	  value += vmono[*it];
	// vpoly.push_back(value);
	vpoly[k++] = value;
      }else{			// second poly
	value = vpoly[comp[1]]*vpoly[comp[2]];
	it++; it++;
	for(;it!=comp.end();it++) 
	  value -= vpoly[*it];
	// vpoly.push_back(value);
	vpoly[k++] = value;
      }
    }
  }
}
