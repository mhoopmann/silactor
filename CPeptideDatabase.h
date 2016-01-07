#pragma once

#include "MascotParser.h"
#include "MSReader.h"
#include "Spectrum.h"
#include <vector>
#include <cmath>

using namespace std;
using namespace MSToolkit;

typedef struct sDBEntry{
  string strName;
  string strSeq;
} sDBEntry;

typedef struct MascotLite {
  int charge;
  double monoMass;
  char peptide[64];
  char protein[32];
  float rTime;        //Light chain rTime
  int scanNum;
  int fileID;
} MascotLite;

typedef struct DBpep {
  char peptide[64];
  char protein[32];
  float firstRT;
  float lastRT;
  int charge;
  double monoMass;
} DBpep;

class CPeptideDatabase{
public:

  //Constructors and destructors
  CPeptideDatabase();
  ~CPeptideDatabase();

  //Operator overrides
  DBpep& operator[ ](const unsigned int i);

  //Data accessors and functions
  void clear();
  void clearMascot();
  void clearPeps();
  size_t size();

  //Analysis functions
  void addPeptide(char* pep, char* prot, float firstRT, float lastRT, int charge, double monoMass);
  bool readMascot(char* mascotFile, char* dataFile);
	bool readPepList(char* fn);
	bool readTextPepList(char* fn);
	void loadDB(char* fname);
  bool buildDB(int fileCount, int pepCount, float rtDrift);
  void exportDB(const char* fname);

protected:
private:

  int fileID;

  vector<MascotLite> mascot;
  vector<DBpep> peps;
	vector<sDBEntry> vDB;

  //Sorting Functions
  void sortPeptide();
  static int comparePeptide(const void *p1, const void *p2);

};
