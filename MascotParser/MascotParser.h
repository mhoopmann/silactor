#ifndef _MASCOTPARSER_H
#define _MASCOTPARSER_H

#include <vector>
#include <iostream>

using namespace std;

typedef struct Mascot {
  char file[256];
	int scanNum;
	int charge;
	char gene[64];
	double intensity;
	int count;
	char sequence_Long[96];
	double zeroMass;
	char sequence[96];
  float rTime;
  float mascotScore;
} Mascot;

typedef struct Mod {
  char aa[20];
  double modMass;
} Mod;

class MascotParser{
private:
	vector<Mascot> m;
  vector<Mod> fixMod;
  vector<Mod> varMod;
  bool bDistiller;
  bool bRTime;

	static int compareScanNum(const void *p1,const void *p2);
protected:
public:

	//Constructors & Destructors
	MascotParser();
	MascotParser(const MascotParser&);
	~MascotParser();

  //Operator Functions
  MascotParser& operator=(const MascotParser& p);
  Mascot& operator[ ](const unsigned int i);

  void add(Mascot& d);
	Mascot& at(unsigned int i);
	void clear();
  void erase(unsigned int i);
	void excel(char *fn);
  bool hasRTime();
  bool isDistiller();
	void load(char *fn);
	void readText(char *fn);
  void reduce();
	void save(char *fn);
	int size();
	void sortScanNum();

  double calcMonoMass(char *seq);

};

#endif

