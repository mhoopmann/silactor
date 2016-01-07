#pragma once

#include "CKronik.h"

using namespace std;

typedef struct uniquePep {
	string seq;
	int charge;
} uniquePep;

typedef struct MRT{
  double mz;
  float rt1;
  float rt2;
} MRT;

typedef struct pepSILAC{
	char peptide[64];
	double ratio;
	double std;
	double rValue;
	float hIntensity;
	float lIntensity;
	float rTime;
	int charge;
} pepSILAC;

class silac {
public:
  int charge;
  int label;          //index to label
  int pepIndex1;      //To Singlets
  int pepIndex2;
  vector<float>* area1;
  vector<float>* area2;
  double monoMass;    //Light chain mono mass used for tracking
  float rTime;        //Light chain rTime
  double da;          //daltons difference in mass between heavy and light
  double ppm;         //ppm error of mass difference
  float rval;
  float pval;
  float corr;
  double mz;
  float slope;
  float intercept;

    //Constructors & Destructor
	silac(){
		area1 = new vector<float>;
		area2 = new vector<float>;
	}
	silac(const silac& s){
		area1 = new vector<float>;
		area2 = new vector<float>;
		for(unsigned int i=0;i<s.area1->size();i++) area1->push_back(s.area1->at(i));
    for(unsigned int i=0;i<s.area2->size();i++) area2->push_back(s.area2->at(i));
    charge=s.charge;
    label=s.label;
    pepIndex1=s.pepIndex1;
    pepIndex2=s.pepIndex2;
    monoMass=s.monoMass;
    rTime=s.rTime;
    da=s.da;
    ppm=s.ppm;
    rval=s.rval;
    pval=s.pval;
    corr=s.corr;
    slope=s.slope;
    intercept=s.intercept;
    mz=s.mz;
	}
	~silac(){
		if(area1) delete area1;
		if(area2) delete area2;
	}

  //Copy operator
	silac& operator=(const silac& s){
		if(this!=&s){
			delete area1;
	    delete area2;
			area1 = new vector<float>;
			area2 = new vector<float>;
			for(unsigned int i=0;i<s.area1->size();i++) area1->push_back(s.area1->at(i));
			for(unsigned int i=0;i<s.area2->size();i++) area2->push_back(s.area2->at(i));
      charge=s.charge;
      label=s.label;
      pepIndex1=s.pepIndex1;
      pepIndex2=s.pepIndex2;
      monoMass=s.monoMass;
      rTime=s.rTime;
      da=s.da;
      ppm=s.ppm;
      rval=s.rval;
      pval=s.pval;
      corr=s.corr;
      slope=s.slope;
      intercept=s.intercept;
      mz=s.mz;
		}
		return *this;
	}
};

//These are all the Kronik peptides. May or may not be
//involved with pairs
class singlet {
public:
  char paired;
  int charge;
  double monoMass;
  string peptide;
  string protein;
  double mz;
  float rTime; 
  vector<float>* area;

  //Constructors & Destructor
	singlet(){
		paired=0;
		area = new vector<float>;
	}
	singlet(const singlet& s){
		area = new vector<float>;
		for(unsigned int i=0;i<s.area->size();i++) area->push_back(s.area->at(i));
		paired=s.paired;
    charge=s.charge;
    monoMass=s.monoMass;
    peptide=s.peptide;
    protein=s.protein;
    rTime=s.rTime;
    mz=s.mz;
	}
	~singlet(){
		if(area) delete area;
	}

  //Copy operator
	singlet& operator=(const singlet& s){
		if(this!=&s){
			delete area;
			area = new vector<float>;
			for(unsigned int i=0;i<s.area->size();i++) area->push_back(s.area->at(i));
		  paired=s.paired;
      charge=s.charge;
      monoMass=s.monoMass;
      peptide=s.peptide;
      protein=s.protein;
      rTime=s.rTime;
      mz=s.mz;
		}
		return *this;
	}

};

class silacProtein {
public:
  string protein;
  float ratio;
  float std;
 
  //Constructors & Destructor
	silacProtein(){}
	silacProtein(const silacProtein& s){
    protein=s.protein;
    ratio=s.ratio;
    std=s.std;
	}
	~silacProtein(){}

  //Copy operator
	silacProtein& operator=(const silacProtein& s){
		if(this!=&s){
			protein=s.protein;
      ratio=s.ratio;
      std=s.std;
		}
		return *this;
	}

};

class CSILACtor{
public:
  //Constructors and destructors
  CSILACtor();
  //CSILACtor(const CKronik2& c);
  CSILACtor(const CSILACtor& c);
  ~CSILACtor();

  //Operator overrides
  CSILACtor& operator=(const CSILACtor& c);
  silac& operator[ ](const unsigned int i);

  //Data accessors
  void clear();
  unsigned int size();

  //Data modifiers
  void addLabel(double d);
	void setCorrThreshold(double d);
  //void setKronik(const CKronik2& c);

  //Analysis Functions
  void addReplicate(CSILACtor& c);
  void analyze(CKronik& kro, float minRT, float maxRT, vector<uniquePep>& v);
  void proteinSummary(char* out, char* target, bool bTrace=false);

  vector<singlet>* peps;
  vector<silacProtein>* proteins;
  double xVal;
	bool bVerbose;

  void reduceAreas();

protected:
private:
  //Data Members
  //CKronik2 kro;
  
  vector<silac>* pairs;
  vector<double>* labels;

  int replicates;
	double corrThreshold;

  //Analysis Functions
  bool checkSILAC(silac& s, bool update);
  void countTryptic(char* seq, int& k, int& r);
  int countMods(char* seq);
  bool matchSeq(char* seq1,char* seq2);
  void scorePair(int index, float& ratio, float& std, float& symRatio, float& symStd, float& fracRatio, float& fracStd);
	void sortPairs();
	void sortSinglets();

	//For searching
	int binarySearch(vector<silac>* p, double mass);
	int binarySearch(vector<singlet>* p, double mass);

	//For sorting
	static int comparePairs(const void *p1, const void *p2);
	static int compareSinglets(const void *p1, const void *p2);

};