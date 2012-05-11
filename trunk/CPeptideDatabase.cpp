#include "CPeptideDatabase.h"

//------------------------------------------
//      Constructors & Destructors
//------------------------------------------
CPeptideDatabase::CPeptideDatabase(){
  fileID=0;
}

CPeptideDatabase::~CPeptideDatabase(){
}

//------------------------------------------
//         Operator Overloads
//------------------------------------------
DBpep& CPeptideDatabase::operator[ ](const unsigned int i){
  return peps[i];
}
//------------------------------------------
//      Data accessors and functions
//------------------------------------------
void CPeptideDatabase::clear(){
  //clears all data
  mascot.clear();
  peps.clear();
}

void CPeptideDatabase::clearMascot(){
  //clears mascot data to release memory
  mascot.clear();
}

void CPeptideDatabase::clearPeps(){
  //clears peptide database but leaves Mascot data intact.
  peps.clear();
}

unsigned int CPeptideDatabase::size(){
  return peps.size();
}

//------------------------------------------
//           Analysis Functions
//------------------------------------------
void CPeptideDatabase::addPeptide(char* pep, char* prot, float firstRT, float lastRT, int charge, double monoMass){
  int i;
  DBpep d;
  bool bMatch=false;

  for(i=0;i<peps.size();i++){
    if(strcmp(peps[i].peptide,pep)==0 && strcmp(peps[i].protein,prot)==0 && peps[i].charge==charge){
      if(peps[i].firstRT>firstRT) peps[i].firstRT=firstRT;
      if(peps[i].lastRT<lastRT) peps[i].lastRT=lastRT;
      bMatch=true;
      break;
    }
  }
  if(!bMatch){
    strcpy(d.peptide,pep);
    strcpy(d.protein,prot);
    d.charge=charge;
    d.firstRT=firstRT;
    d.lastRT=lastRT;
    d.monoMass=monoMass;
    peps.push_back(d);
    //cout << "Added " << pep << " " << d.peptide << " -- " << prot << " " << d.protein << " " << d.charge << endl;
  }
}

bool CPeptideDatabase::readMascot(char* mascotFile, char* dataFile){
  MSReader r;
  MascotParser p;
  Spectrum s;
  MascotLite m;

  int i;
  int percent=0;

  p.clear();
  p.readText(mascotFile);

  if(p.isDistiller()){
    cout << mascotFile << "\t" << p.size() << " peptides from Distiller, no need for RT lookup." << endl;
  } else if(p.hasRTime()){
    cout << mascotFile << "\t" << p.size() << " peptides have RT already, no need for lookup." << endl;
  } else {
    cout << mascotFile << "\tLooking up retention times for " << p.size() << " Mascot peptides..." << endl;
    cout << percent;
    r.readFile(dataFile,s);
  }

  for(i=0;i<p.size();i++){
    if(p.isDistiller()||p.hasRTime()){ 
      m.rTime=p[i].rTime;
    } else {
      if( (int)i*100/p.size() > percent){
        percent=(int)i*100/p.size();
        cout << "\b\b\b" << percent;
      }

      r.readFile(NULL,s,p[i].scanNum);
      if(s.getScanNumber()==0){
        cout << "Scan not found. Please check that correct data file is used. Stopping Mascot parsing." << endl;
        return false;
      }
      m.rTime=s.getRTime();
    }
    m.charge=p[i].charge;
    m.fileID=fileID;
    m.monoMass=p[i].zeroMass;
    m.scanNum=p[i].scanNum;
    strcpy(m.peptide,p[i].sequence_Long);
    strcpy(m.protein,p[i].gene);
    mascot.push_back(m);
  }
  cout << "  Done!" << endl;

  fileID++;

  return true;
}

bool CPeptideDatabase::readPepList(char* fn){
	FILE* f=fopen(fn,"rt");

	char pep[64];
	float rt1,rt2;
	int charge;
	double mass;
	char basePep[64];
	char name[32];
	int i,j;

	bool bFound;

	while(!feof(f)){
		fscanf(f,"%s\t%f\t%f\t%d\t%lf\n",pep,&rt1,&rt2,&charge,&mass);

		//match base peptide to database;
		strcpy(basePep,"");
		for(i=0;i<strlen(pep);i++){
			if(isalpha(pep[i])) strncat(basePep,&pep[i],1);
		}

		//iterate through all proteins for match
		bFound=false;
		for(i=0;i<vDB.size();i++){
			if(strstr(&vDB[i].strSeq[0],basePep)!=NULL){
				bFound=true;
				strncpy(name,&vDB[i].strName[0],31);
				name[31]='\0';
				addPeptide(pep,name,rt1,rt2,charge,mass);
			}
		}
		if(!bFound) cout << basePep << " not found." << endl;

	}

	fclose(f);
	return true;

}

bool CPeptideDatabase::readTextPepList(char* fn){

	int i;
  FILE* f;
	char str[2048];
	char* tok;
	DBpep d;

	int iPepIndex=-1;
	int iProtIndex=-1;
	int iFirstRT=-1;
	int iLastRT=-1;
	int iCharge=-1;
	int iMass=-1;


	f=fopen(fn,"rt");
	if(f==NULL){
		cout << "ERROR: Cannot open " << fn << endl;
		return false;
	}

	//read file header
	i=0;
	fgets(str,2048,f);
	tok=strtok(str,"\t\n");
	while(tok!=NULL){
		if(strcmp(tok,"peptide")==0) iPepIndex=i;
		if(strcmp(tok,"charge")==0) iCharge=i;
		if(strcmp(tok,"protein")==0) iProtIndex=i;
		if(strcmp(tok,"firstRT")==0) iFirstRT=i;
		if(strcmp(tok,"lastRT")==0) iLastRT=i;
		if(strcmp(tok,"mass")==0) iMass=i;
		i++;
		tok=strtok(NULL,"\t\n");
	}

	//Check that all columns are present:
	if(iPepIndex<0) {
		cout << "ERROR: Dabase missing 'peptide' column." << endl;
		return false;
	}
	if(iProtIndex<0) {
		cout << "ERROR: Dabase missing 'protein' column." << endl;
		return false;
	}
	if(iCharge<0) {
		cout << "ERROR: Dabase missing 'charge' column." << endl;
		return false;
	}
	if(iFirstRT<0) {
		cout << "ERROR: Dabase missing 'firstRT' column." << endl;
		return false;
	}
	if(iLastRT<0) {
		cout << "ERROR: Dabase missing 'lastRT' column." << endl;
		return false;
	}
	if(iMass<0) {
		cout << "ERROR: Dabase missing 'mass' column." << endl;
		return false;
	}

	//read in peptides
	peps.clear();
	while(!feof(f)){
		fgets(str,2048,f);
		i=0;
		tok=strtok(str,"\t\n");
		while(tok!=NULL){
			if(i==iPepIndex) strcpy(d.peptide,tok);
			if(i==iProtIndex) strcpy(d.protein,tok);
			if(i==iCharge) d.charge=atoi(tok);
			if(i==iFirstRT) d.firstRT=(float)atof(tok);
			if(i==iLastRT) d.lastRT=(float)atof(tok);
			if(i==iMass) d.monoMass=atof(tok);
			i++;
			tok=strtok(NULL,"\t\n");
		}
		peps.push_back(d);
	}
	fclose(f);

	return true;
}

bool CPeptideDatabase::buildDB(int fileCount, int pepCount, float rtDrift){
  int i,j,k,n;
  float lowRT;
  float highRT;
  bool bMatch;
  vector<int> fc;
  DBpep d;

  sortPeptide();
  i=0;
  /*
  while(i<mascot.size()){
    cout << mascot[i].peptide << "\t" << mascot[i].charge << endl;
    i++;
  }
  exit(0);
  */

  while(i<mascot.size()){
    j=i+1;
    while(j<mascot.size()){
      if(strcmp(mascot[i].peptide,mascot[j].peptide) || mascot[i].charge!=mascot[j].charge) break;
      j++;
    }

    //check pepCount
    if(j-i<pepCount){
      i=j;
      continue;
    }

    //check retention time
    lowRT=mascot[i].rTime;
    highRT=mascot[i].rTime;
    for(k=i+1;k<j;k++){
      if(mascot[k].rTime<lowRT) lowRT=mascot[k].rTime;
      if(mascot[k].rTime>highRT) highRT=mascot[k].rTime;
    }
    if(highRT-lowRT>rtDrift){
      i=j;
      continue;
    }

    //check file count
    for(k=i;k<j;k++){
      bMatch=false;
      for(n=0;n<fc.size();n++){
        if(fc[n]==mascot[k].fileID) bMatch=true;
      }
      if(!bMatch) fc.push_back(mascot[k].fileID);
    }
    if(fc.size()<fileCount){
      i=j;
      continue;
    }

    //Add peptide
    d.charge=mascot[i].charge;
    d.firstRT=lowRT;
    d.lastRT=highRT;
    d.monoMass=mascot[i].monoMass;
    strcpy(d.peptide,mascot[i].peptide);
    strcpy(d.protein,mascot[i].protein);
    peps.push_back(d);

    i=j;

  }

  cout << peps.size() << " consistent peptides from a total of " << mascot.size() << " entries." << endl;
	
	double avg=0.0;
	double std=0.0;
	double dif=0.0;
	for(i=0;i<peps.size();i++)	avg+=(peps[i].lastRT-peps[i].firstRT);
	avg/=peps.size();
	for(i=0;i<peps.size();i++) dif+=pow(( (peps[i].lastRT-peps[i].firstRT) -avg),2);
  std=sqrt(dif/peps.size());
	cout << "Average retention time length: " << avg << " +/- " << std << endl;

	//clear mascot entries to save memory
	mascot.clear();

  return true;
}

void CPeptideDatabase::exportDB(char* fname){
  FILE *f=fopen(fname,"wt");
  int i;

  for(i=0;i<peps.size();i++){
    fprintf(f,"%s\t%s\t%.6lf\t%d\t%.4f\t%.4f\n",peps[i].protein,peps[i].peptide,peps[i].monoMass,peps[i].charge,peps[i].firstRT,peps[i].lastRT);
  }
  fclose(f);
}

//------------------------------------------
//           Sorting Functions
//------------------------------------------
void CPeptideDatabase::sortPeptide(){
	qsort(&mascot[0],mascot.size(),sizeof(MascotLite),comparePeptide);
}
int CPeptideDatabase::comparePeptide(const void *p1, const void *p2){
  const MascotLite d1 = *(MascotLite *)p1;
  const MascotLite d2 = *(MascotLite *)p2;
  if(strcmp(d1.peptide,d2.peptide)==0){
    if(d1.charge==d2.charge) return 0;
    else if(d1.charge<d2.charge) return -1;
    else return 1;
  } else {
    return strcmp(d1.peptide,d2.peptide);
  }
}



void CPeptideDatabase::loadDB(char* fname){

  //char szBuf[8192];
  FILE *fptr;
  sDBEntry dbe;
  int iTmpCh;

  if ((fptr=fopen(fname, "rb")) == NULL) {
		cout << "File open failed: " << fname  << endl;
    exit(-1);
  }

  //fgets(szBuf, 8192, fptr);
  iTmpCh=getc(fptr);

  while(!feof(fptr))  {
    dbe.strName="";
    dbe.strSeq="";
    if (iTmpCh!='>')  {
      cout << "Error" << endl;
      break;
    }
    while( (iTmpCh=getc(fptr))!='\n' && iTmpCh!='\r' && iTmpCh!=EOF ) {
      if(dbe.strName.size()<128) dbe.strName+=iTmpCh;
    }
    while( (iTmpCh=getc(fptr))!='>' && iTmpCh!=EOF ) {
      if (isalnum(iTmpCh))  dbe.strSeq+=toupper(iTmpCh);
    }
		//if(strstr(&dbe.strName[0],"RAND")!=NULL) vDB.push_back(dbe);
		vDB.push_back(dbe);
  }
  fclose(fptr);

}