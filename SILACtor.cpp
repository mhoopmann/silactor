#include "CKronik2.h"
#include "MascotParser.h"
#include "CSILACtor.h"
#include "CPeptideDatabase.h"
#include <vector>

using namespace std;

typedef struct files {
  string mascotFile;
  string dataFile;
  string hkFile;
  string outFile;
  double xVal;
  int ID;
} files;

typedef struct replicates {
  string outFile;
  string listFile;
  int ID;
  vector<int> IDs;
} replicates;

typedef struct singleScans {
  CSILACtor ss;
  int ID;
} singleScans;

/*
typedef struct sAnalysis {
  char cType;
  string outFile;
  int ID;
  vector<
  */

int binarySearch(CKronik2& p, double mass);
bool checkSILAC(silac& s,int index);
void matchMascot(CKronik2& p, MascotParser& mp);
void matchMascot2(CKronik2& p, MascotParser& mp);
void matchMascot3(CKronik2& p, CPeptideDatabase& mp, vector<uniquePep>& v);
bool parser(char* in);
void SILACtor(CKronik2& p, MascotParser& mp, int index);
void timePoint(int index);
//vector<int> pepID;
//vector<silac> hits;
vector<int> masterID;

vector<double> silacMass;
vector<files> myDB;
vector<files> mySILAC;
vector<replicates> myReps;
vector<replicates> myTime;

float firstRT;
float lastRT;
float proteinStd;

int dbMinFileCount;
int dbMinEventCount;
float dbRTTolerance;
double dbMatchTolerance;

double corrThreshold;

bool bDBx;
bool bDB;
bool bDBA;
bool bVerbose;

vector<singleScans> singles;

int main(int argc, char* argv[]){

	cout << "SILACtor version 1.05, Build date: May 9 2012" << endl;
	cout << "Copyright 2010-2012, Michael R. Hoopmann, James E. Bruce" << endl;
	cout << "University of Washington, Department of Genome Sciences" << endl;
	cout << "http://brucelab.gs.washington.edu/SILACtor.php" << endl;
	cout << "\nSee documentation for usage details." << endl;

	if(argc!=2){
		exit(1);
	}		

  //argv[1] = Hardklor output

  CKronik2 p;
  MascotParser mp;
  singleScans s;
  CPeptideDatabase d;
  
  int i,j,k;

  //Initialize variables
  firstRT=0.0;
  lastRT=300.0;
  dbMinFileCount=1;
  dbMinEventCount=1;
  dbRTTolerance=10.0f;
	dbMatchTolerance=10.0;
	corrThreshold=0.9;

  //Parse settings
	bDBx=false;
	bDB=false;
	bDBA=false;
	bVerbose=false;
  parser(argv[1]);

  //Set silac masses
  for(i=0;i<silacMass.size();i++) s.ss.addLabel(silacMass[i]);

  //Build peptide sequence database
	if(bDBx){
		d.loadDB(&myDB[0].dataFile[0]);
		d.readPepList(&myDB[0].mascotFile[0]);
	} else {
		if(bDB && bDBA){
			cout << "ERROR: Configuration file uses both Mascot and AMT libraries (DB and DBA tags). Please use only one." << endl;
			exit(10);
		} else if(bDBA){
			d.readTextPepList(&myDB[0].dataFile[0]);
		} else {
			for(i=0;i<myDB.size();i++){
				d.readMascot(&myDB[i].mascotFile[0],&myDB[i].dataFile[0]);
			}
			d.buildDB(dbMinFileCount,dbMinEventCount,dbRTTolerance);
		}
	}

	/* This is for using a specific peptide database instead
  FILE* f=fopen("pepDatabase.txt","rt");
  char str1[64];
  char str2[64];
  double d1;
  float f1;
  float f2;
  if(f!=NULL){
    while(!feof(f)){
      fscanf(f,"%s\t%s\t%lf\t%d\t%f\t%f\n",str1,str2,&d1,&i,&f1,&f2);
      if(feof(f)) continue;
      d.addPeptide(str2,str1,f1,f2,i,d1);
    }
    fclose(f);
  }

	cout << d.size() << endl;
	*/

  d.exportDB("pepDatabase.txt");
  
	vector<uniquePep> vu;
  //Analyze all data files
  for(i=0;i<mySILAC.size();i++){
    cout << "\nAnalyzing " << &mySILAC[i].hkFile[0] << endl;

    p.clear();
    p.setGapTol(1);
    p.setMatchTol(2);
    p.setPPMTol(8.0);
    p.processHK(&mySILAC[i].hkFile[0],"krout.kro");
    matchMascot3(p,d,vu);

    //Run SILACtor
    s.ss.clear();
		s.ss.bVerbose=bVerbose;
    s.ID=mySILAC[i].ID;
    cout << "Performing SILAC analysis..." << endl;
		s.ss.setCorrThreshold(corrThreshold);
    s.ss.analyze(p,firstRT,lastRT,vu);
    s.ss.proteinSummary(&mySILAC[i].outFile[0],"",false);
    s.ss.xVal=mySILAC[i].xVal;
    singles.push_back(s);
    cout << "Done!" << endl;
  }

	//clear last Kronik to save memory
	p.clear();

  //Combine SILAC results into replicates
  for(i=0;i<myReps.size();i++){
    s.ID=myReps[i].ID;
    for(j=0;j<singles.size();j++){
      if(singles[j].ID==myReps[i].IDs[0]) {
        s.ss=singles[j].ss;
        break;
      }
    }

    cout << "\nCombining replicates.." << endl;
    for(j=1;j<myReps[i].IDs.size();j++){
      for(k=0;k<singles.size();k++){
        if(singles[k].ID==myReps[i].IDs[j]) {
          s.ss.addReplicate(singles[k].ss);
          break;
        }
      }
    }
    s.ss.proteinSummary(&myReps[i].outFile[0],&myReps[i].listFile[0],false);
    singles.push_back(s);
    cout << "Done!" << endl;
  }

  //Do any timepoint analyses
  for(i=0;i<myTime.size();i++) timePoint(i);

  return 0;
}

void matchMascot2(CKronik2& p, MascotParser& mp){
  int i,j;
  double ppm;
  int match=0;
  int multi=0;
  int difPep=0;
  int difProt=0;
  int err=0;
  bool bMatch=false;
  p.sortMonoMass();

  double avg=0.0;
  double dif=0.0;
  double std=0.0;

  vector<double> vm;

  for(i=0;i<mp.size();i++){
    bMatch=false;
    for(j=0;j<p.size();j++){
      if(mp[i].charge!=p[j].charge) continue;
      ppm = (p[j].monoMass-mp[i].zeroMass)/mp[i].zeroMass * 1000000;

      //Change this to use RT - must gather from RAW files...
      if(fabs(ppm)<10.0 && (mp[i].scanNum > p[j].lowScan-80) && (mp[i].scanNum < p[j].highScan+40)){

        vm.push_back(ppm);

        if(!bMatch){
          bMatch=true;
          match++;
        } else {
          multi++;
        }
        if(!strcmp(p[j].sequence,"")){
          strcpy(p[j].sequence,mp[i].sequence_Long);
          strcpy(p[j].gene,mp[i].gene);  
        } else if(!strcmp(p[j].sequence,mp[i].sequence_Long)){
          if(!strcmp(p[j].gene,mp[i].gene)) difProt++;
          else err++;
        } else {
          difPep++;
        }
      }
    }
    //if(!bMatch) cout << mp[i].scanNum << " " << mp[i].sequence_Long << " " << mp[i].zeroMass << endl;
  }

  for(i=0;i<vm.size();i++) avg+=vm[i];
  avg/=vm.size();
  for(i=0;i<vm.size();i++) dif+=pow((vm[i]-avg),2);
  std=sqrt(dif/vm.size());
  cout << avg << " +/- " << std << endl;

  cout << match << " of " << mp.size() << " (" << (double)(match)/(double)mp.size()*100.0 << "%) were matched to a precursor peptide ion." << endl;
  cout << difPep+difProt << " matches were to an already identified peptide." << endl;
  cout << difProt << " matches were to different proteins." << endl;
  cout << multi << " sequences matched to multiple precursors." << endl;
  cout << err << " redundantly identified peptides from Mascot." << endl;

}

void matchMascot3(CKronik2& p, CPeptideDatabase& mp, vector<uniquePep>& v){
  int i,j;
  double ppm;
  int match=0;
  int multi=0;
  int difPep=0;
  int difProt=0;
  int err=0;
  bool bMatch=false;
  p.sortMonoMass();

  double avg=0.0;
  double dif=0.0;
  double std=0.0;

  vector<double> vm;
	v.clear();
	uniquePep up;

	unsigned int k;
	int pos,lower,upper;

	//Sort peptides by monoMass for binary searching
	p.sortMonoMass();

  for(i=0;i<mp.size();i++){
    bMatch=false;

		//Find boundaries to search by
		//1) finding closest mass
		//2) extending boundaries until ppm cutoff
		pos=binarySearch(p,mp[i].monoMass);
		if(pos==-1) continue;
		lower=pos;
		upper=pos;
		for(j=pos-1;j>-1;j--){
			ppm = (p[j].monoMass-mp[i].monoMass)/mp[i].monoMass * 1000000;
			if(ppm>-dbMatchTolerance) lower=j;
			else break;
		}
		for(j=pos+1;j<p.size();j++){
			ppm = (p[j].monoMass-mp[i].monoMass)/mp[i].monoMass * 1000000;
			if(ppm<dbMatchTolerance) upper=j;
			else break;
		}

    for(j=lower;j<=upper;j++){
      if(mp[i].charge!=p[j].charge) continue;     

			//already checked ppm
      //Change this to use RT - must gather from RAW files...
      if( /*fabs(ppm)<10.0 && */ p[j].rTime>mp[i].firstRT-1 && p[j].rTime<mp[i].lastRT+1){

				ppm = (p[j].monoMass-mp[i].monoMass)/mp[i].monoMass * 1000000;
        vm.push_back(ppm);

				//cout << mp[i].peptide << " " << mp[i].monoMass << " " << p[j].monoMass << " " << ppm << endl;

        if(!bMatch){
          match++;
					up.seq=mp[i].peptide;
					up.charge=mp[i].charge;
					for(k=0;k<v.size();k++){
						if(v[k].seq.compare(up.seq)==0 && v[k].charge==up.charge){
							bMatch=true;
							break;
						}
					}
					if(!bMatch) v.push_back(up);
					bMatch=true;
        } else {
          multi++;
        }
        if(!strcmp(p[j].sequence,"")){
          strcpy(p[j].sequence,mp[i].peptide);
          strcpy(p[j].gene,mp[i].protein);  
        } else if(!strcmp(p[j].sequence,mp[i].peptide)){
          if(!strcmp(p[j].gene,mp[i].protein)) difProt++;
          else err++;
        } else {
          //cout << p[j].sequence << " " << mp[i].peptide << endl;
          difPep++;
        }
      }
    }
		//if(!bMatch) cout << mp[i].firstRT << " " << mp[i].peptide << " " << mp[i].monoMass << endl;
  }

  for(i=0;i<vm.size();i++) avg+=vm[i];
  avg/=vm.size();
  for(i=0;i<vm.size();i++) dif+=pow((vm[i]-avg),2);
  std=sqrt(dif/vm.size());
  cout << "Mass accuracy of peptide precursors (when matched to sequence): " << avg << " +/- " << std << endl;

  cout << match << " of " << mp.size() << " (" << (double)(match)/(double)mp.size()*100.0 << "%) were matched to a precursor peptide ion." << endl;
  cout << difPep+difProt << " matches were to an already identified peptide." << endl;
  cout << difProt << " matches were to different proteins." << endl;
  cout << multi << " sequences matched to multiple precursors." << endl;
  cout << err << " redundantly identified peptides from Mascot." << endl;

}

int binarySearch(CKronik2& p, double mass){
	int lower,mid,upper;
	int sz=p.size();
	double lowMass=mass-(mass/1000000*dbMatchTolerance);
	double highMass=mass+(mass/1000000*dbMatchTolerance);

	if(sz<1) return -1;

	//Find a mass within ppm tolerance. Use binary search to for speed.
	//Stop searching when first mass is found (although a second may exist).
	mid=sz/2;
	lower=0;
	upper=sz;
	
	while(p[mid].monoMass<lowMass || p[mid].monoMass>highMass){
		if(lower>=upper) break;
		if(mass<p[mid].monoMass){
			upper=mid-1;
			mid=(lower+upper)/2;
		} else {
			lower=mid+1;
			mid=(lower+upper)/2;
		}
		if(mid==sz) {
			mid--;
			break;
		}
	}

	//Check that mass is correct
	if(p[mid].monoMass<lowMass || p[mid].monoMass>highMass) return -1;
	return mid;
}



bool parser(char* in) {
  FILE* f;
  char str[512];
  char* tok;
  string s,s2;
  int i,j;
  bool bGood=true;
  bool bMatch;
  files myFile;
  replicates rep;

  f=fopen(in,"rt");
  if(f==NULL) {
    cout << "Cannot open " << in << " to parse." << endl;
    return false;
  }

  while(!feof(f)){
    if(fgets(str,512,f)==NULL) continue;

    tok=strtok(str,"\t\n");
    if(tok==NULL) continue;
    if(tok[0]=='#') continue;

    //Read in database files
    if(!strcmp(tok,"DB") || !strcmp(tok,"DBx")){
			if(!strcmp(tok,"DBx")) bDBx=true;
      bGood=true;
      tok=strtok(NULL,"\t\n");
      if(tok==NULL) bGood=false;
      else myFile.mascotFile=tok;
      tok=strtok(NULL,"\t\n");
      if(tok==NULL) bGood=false;
      else myFile.dataFile=tok;

      if(bGood) myDB.push_back(myFile);
      else cout << "Bad DB line: " << str << endl;
			bDB=true;

		} else if(!strcmp(tok,"DBA")) {
			//Read in AMT library file
			tok=strtok(NULL,"\t\n");
      if(tok==NULL) cout << "Bad DBA line" << endl;
			else myFile.dataFile=tok;
			myFile.mascotFile="";

			myDB.push_back(myFile);
			bDBA=true;

      //Read in silac masses
    } else if(!strcmp(tok,"S")){
      tok=strtok(NULL,"\t\n");
      if(tok==NULL) cout << "Bad S line" << endl;
      else silacMass.push_back(atof(tok));

      //Read in files for SILAC analysis
    } else if(!strcmp(tok,"F")){
      bGood=true;
      tok=strtok(NULL,"\t\n");
      if(tok==NULL) bGood=false;
      else {
        myFile.ID=atoi(tok);
        for(i=0;i<masterID.size();i++){
          if(myFile.ID==masterID[i]) bGood=false;
        }
      }
      tok=strtok(NULL,"\t\n");
      if(tok==NULL) bGood=false;
      else myFile.hkFile=tok;
      tok=strtok(NULL,"\t\n");
      if(tok==NULL) bGood=false;
      else myFile.outFile=tok;
      tok=strtok(NULL,"\t\n");
      if(tok==NULL) bGood=false;
      else myFile.xVal=atof(tok);

      if(bGood) {
        mySILAC.push_back(myFile);
        masterID.push_back(myFile.ID);
      } else cout << "Bad F line: " << str << endl;

      //Read in parameters
    } else if(!strcmp(tok,"P")){
      tok=strtok(NULL,"\t\n");
      if(tok==NULL){
        cout << "Bad P line" << endl;
        continue;
      } else {
        if(!strcmp(tok,"RTime")){
          tok=strtok(NULL,"\t\n");
          if(tok==NULL) {
            cout << "Bad RTime value." << endl;
            continue;
          }
          firstRT=(float)atof(tok);
          tok=strtok(NULL,"\t\n");
          if(tok==NULL) {
            cout << "Bad RTime value." << endl;
            continue;
          }
          lastRT=(float)atof(tok);
        } else if(!strcmp(tok,"ProteinStDev")){
          tok=strtok(NULL,"\t\n");
          if(tok==NULL) {
            cout << "Bad ProteinStDev value." << endl;
            continue;
          }
          proteinStd=(float)atof(tok);
				} else if(!strcmp(tok,"DBMatchTolerance")){
          tok=strtok(NULL,"\t\n");
          if(tok==NULL) {
            cout << "Bad DBMatchTolerance value." << endl;
            continue;
          }
          dbMatchTolerance=atof(tok);
        } else if(!strcmp(tok,"DBMinFileCount")){
          tok=strtok(NULL,"\t\n");
          if(tok==NULL) {
            cout << "Bad DBMinFileCount value." << endl;
            continue;
          }
          dbMinFileCount=atoi(tok);
        } else if(!strcmp(tok,"DBMinEventCount")){
          tok=strtok(NULL,"\t\n");
          if(tok==NULL) {
            cout << "Bad DBMinEventCount value." << endl;
            continue;
          }
          dbMinEventCount=atoi(tok);
        } else if(!strcmp(tok,"DBRTTolerance")){
          tok=strtok(NULL,"\t\n");
          if(tok==NULL) {
            cout << "Bad DBRTTolerance value." << endl;
            continue;
          }
          dbRTTolerance=(float)atof(tok);
				} else if(!strcmp(tok,"CorrThreshold")){
          tok=strtok(NULL,"\t\n");
          if(tok==NULL) {
            cout << "Bad CorrThreshold value." << endl;
            continue;
          }
          corrThreshold=atof(tok);
				} else if(!strcmp(tok,"Verbose")){
					bVerbose=true;
        } else {
          cout << "Unknown parameter: " << tok << endl;
        }
      }
      
      //Parse replicate analyses
    } else if(!strcmp(tok,"R")){
      bGood=true;
      rep.IDs.clear();
      tok=strtok(NULL,"\t\n");
      if(tok==NULL) bGood=false;
      else {
        rep.ID=atoi(tok);
        for(i=0;i<masterID.size();i++){
          if(rep.ID==masterID[i]) bGood=false;
        }
      }
      tok=strtok(NULL,"\t\n");
      if(tok==NULL) bGood=false;
      else rep.outFile=tok;
      tok=strtok(NULL,"\t\n");
      if(tok==NULL) bGood=false;
      else rep.listFile=tok;

      tok=strtok(NULL,"\t\n");
      while(tok!=NULL){
        i=atoi(tok);
        bMatch=false;
        for(j=0;j<mySILAC.size();j++){
          if(i==mySILAC[j].ID) {
            rep.IDs.push_back(i);
            bMatch=true;
          }
        }
        if(!bMatch) bGood=false;
        tok=strtok(NULL,"\t\n");
      }

      if(bGood && rep.IDs.size()>1) {
        myReps.push_back(rep);
        masterID.push_back(rep.ID);
      } else cout << "Bad R line: " << str << endl;

      //Parse timepoint analyses
    } else if(!strcmp(tok,"T")){
      bGood=true;
      rep.IDs.clear();
      tok=strtok(NULL,"\t\n");
      if(tok==NULL) bGood=false;
      else {
        rep.ID=atoi(tok);
        for(i=0;i<masterID.size();i++){
          if(rep.ID==masterID[i]) bGood=false;
        }
      }
      tok=strtok(NULL,"\t\n");
      if(tok==NULL) bGood=false;
      else rep.outFile=tok;

      tok=strtok(NULL,"\t\n");
      while(tok!=NULL){
        i=atoi(tok);
        bMatch=false;
        for(j=0;j<mySILAC.size();j++){
          if(i==mySILAC[j].ID) {
            rep.IDs.push_back(i);
            bMatch=true;
          }
        }
        for(j=0;j<myReps.size();j++){
          if(i==myReps[j].ID) {
            rep.IDs.push_back(i);
            bMatch=true;
          }
        }
        if(!bMatch) bGood=false;
        tok=strtok(NULL,"\t\n");
      }

      if(bGood && rep.IDs.size()>1) {
        myTime.push_back(rep);
        masterID.push_back(rep.ID);
      } else cout << "Bad T line: " << str << endl;

    }
  }

  fclose(f);
  return true;

}

void timePoint(int index){
  int i,j,k,n;
  FILE* f;
  vector<string> prots;
  string protein;
  vector<int> sIndex;
  bool bMatch;

  vector<int> analIndex;
  vector<int> protIndex;
  int matchCount;

  vector<double> yVal;
  vector<double> xVal;

	vector<double> timeOut;
	vector<double> riaOut;
	vector<double> stdOut;

  CSILACtor s,s2;

  double sumx=0.0;
  double sumy=0.0;
  double sumxx=0.0;
  double sumxy=0.0;
  double avgx=0.0;
  double avgy=0.0;
  double b=0.0;
  double A=0.0;
  double R2=0.0;
  double sqx=0.0;
  double sqy=0.0;
  double sqxy=0.0;
  double sqye=0.0;

  cout << "Performing timepoint...";

  f=fopen(&myTime[index].outFile[0],"wt");
  if(f==NULL) {
    cout << "Cannot open " << &myTime[index].outFile[0] << endl;
    return;
  }
	fprintf(f,"Protein\tKloss\tR^2\tTimePoints\tRIA\tStDev\n");

  FILE *f2=fopen("timepoints2.txt","wt");
  fprintf(f2,"Protein\tKloss\tR^2\n");

	//Set time points for output
	for(i=0;i<singles.size();i++){
		timeOut.push_back(singles[i].ss.xVal);
	}
	for(i=0;i<timeOut.size()-1;i++){
		for(j=i+1;j<timeOut.size();j++){
			if(timeOut[j]<timeOut[i]){
				b=timeOut[j];
				timeOut[j]=timeOut[i];
				timeOut[i]=b;
			}
		}
	}
	for(i=0;i<timeOut.size()-1;i++){
		for(j=i+1;j<timeOut.size();j++){
			if(timeOut[j]==timeOut[i]){
				timeOut.erase(timeOut.begin()+j);
				j--;
			}
		}
	}

  //Get indexes from singles
  for(i=0;i<myTime[index].IDs.size();i++){
    for(j=0;j<singles.size();j++){
      if(singles[j].ID==myTime[index].IDs[i]) sIndex.push_back(j);
    }
  }

  //Check all requested silac analyses
  for(i=0;i<sIndex.size();i++){

    if(i==0) {
      s=singles[sIndex[i]].ss;
      s.reduceAreas();
    } else {
      s2=singles[sIndex[i]].ss;
      s2.reduceAreas();
      s.addReplicate(s2);
    }

    //Check all proteins
    for(j=0;j<singles[sIndex[i]].ss.proteins->size();j++){

      //Get protein name
      protein=singles[sIndex[i]].ss.proteins->at(j).protein;

      //Check if protein already analyzed
      bMatch=false;
      for(k=0;k<prots.size();k++){
        if(prots[k].compare(protein)==0) {
          bMatch=true;
          break;
        }
      }
      if(bMatch) continue;

      analIndex.clear();
      protIndex.clear();
      matchCount=0;

      k=0;
      while(k<i){
        analIndex.push_back(-1);
        protIndex.push_back(-1);
        k++;
      }

      //Output protein score
      analIndex.push_back(sIndex[i]);
      protIndex.push_back(j);
      matchCount++;

      //Check remaining replicates
      for(k=i+1;k<sIndex.size();k++){

        //Check proteins for match
        bMatch=false;
        for(n=0;n<singles[sIndex[k]].ss.proteins->size();n++){

          if(protein.compare(singles[sIndex[k]].ss.proteins->at(n).protein)==0){
            bMatch=true;
            analIndex.push_back(sIndex[k]);
            protIndex.push_back(n);
            matchCount++;
            break;
          }
        }
        if(!bMatch) {
          analIndex.push_back(-1);
          protIndex.push_back(-1);
        }

      }

      prots.push_back(protein);

      if(matchCount>2){
        xVal.clear();
        yVal.clear();
				riaOut.clear();
				stdOut.clear();
        //fprintf(f,"\n%s:\n",&protein[0]);
        for(n=0;n<analIndex.size();n++){
          if(analIndex[n]==-1){
						riaOut.push_back(0.0);
						stdOut.push_back(0.0);
          } else {
            //fprintf(f,"%.0lf\t%.4f\t%.4f\n",singles[analIndex[n]].ss.xVal,singles[analIndex[n]].ss.proteins[protIndex[n]].ratio,singles[analIndex[n]].ss.proteins[protIndex[n]].std);
            xVal.push_back(singles[analIndex[n]].ss.xVal);
            yVal.push_back(log(singles[analIndex[n]].ss.proteins->at(protIndex[n]).ratio));
						riaOut.push_back(singles[analIndex[n]].ss.proteins->at(protIndex[n]).ratio);
						stdOut.push_back(singles[analIndex[n]].ss.proteins->at(protIndex[n]).std);
          }
        }

        //Compute R^2 and function
        //Get ln for each value
        sumx=0.0;
        sumy=0.0;
        sumxx=0.0;
        sumxy=0.0;
        sqx=0.0;
        sqy=0.0;
        sqxy=0.0;
        sqye=0.0;
        for(n=0;n<xVal.size();n++){
          sumx+=xVal[n];
          sumy+=yVal[n];
          sumxy+=(xVal[n]*yVal[n]);
          sumxx+=(xVal[n]*xVal[n]);
        }
        n=xVal.size();
        avgx=sumx/n;
        avgy=sumy/n;

        b = sumxy / sumxx;
        
        //Comput R2
        for(k=0;k<n;k++){
          sqy  += (yVal[k]-avgy) * (yVal[k]-avgy);
          sqye += (yVal[k]-b*xVal[k]) * (yVal[k]-b*xVal[k]);
        }
        R2 = 1.0 - sqye/sqy;

        //fprintf(f,"y=e^%.4lfx\t%.4lf\t%.4lf\n",b,R2,-b);
        fprintf(f2,"%s\t%.8lf\t%.4lf\n",&protein[0],-b,R2);

				fprintf(f,"%s\t%.4lf\t%.4lf\t",&protein[0],-b,R2);
				for(k=0;k<timeOut.size()-1;k++) fprintf(f,"%.1lf,",timeOut[k]);
				fprintf(f,"%.1lf\t",timeOut[timeOut.size()-1]);
				for(k=0;k<riaOut.size()-1;k++) fprintf(f,"%.4lf,",riaOut[k]);
				fprintf(f,"%.4lf\t",riaOut[riaOut.size()-1]);
				for(k=0;k<stdOut.size()-1;k++) fprintf(f,"%.4lf,",stdOut[k]);
				fprintf(f,"%.4lf\n",stdOut[stdOut.size()-1]);

      }

    }
  }

  fclose(f);

  cout << "Done!" << endl;

  s.proteinSummary("dump.txt","",true);


}



/*

void SILACtor(CKronik2& p, MascotParser& mp, int index){

  int i,j;
  int count=0;
  double dif;
  double ppm;

  double avgInc=0.0;
  int avgCount=0;
  silac s;

  CKronik2 tp;
  //double rval;
  //double pval;
  double slope;
  double intercept;
  //double ratio;

  for(i=0;i<p.size();i++){
    for(j=i-100;j<p.size();j++){

      s.ratio[0]=s.ratio[1]=s.ratio[2]=0;
      s.pval[0]=s.pval[1]=s.pval[2]=0;
      s.rval[0]=s.rval[1]=s.rval[2]=0;
      s.intensity[0]=s.intensity[1]=s.intensity[2]=0;
      
      if(i==j) continue;
      if(j<0) continue;

      if(p[i].charge!=p[j].charge) continue;
      if(fabs(p[i].rTime - p[j].rTime) > 0.50) continue;

      dif = p[j].monoMass-6.0201324-p[i].monoMass;
      
      ppm = dif/p[i].monoMass*1000000;
      if(fabs(ppm)<5.0) {
        tp.clear();
        tp.add(p[i]);
        tp.add(p[j]);

        tp.pearson(0,1,true,true,s.rval[index],s.pval[index],slope,intercept,s.ratio[index]);

        s.charge=p[i].charge;
        s.monoMass=p[i].monoMass;
        s.ppm=ppm;
        s.da=p[j].monoMass-p[i].monoMass;
        s.rTime=p[i].rTime;
        s.intensity[index]=p[i].intensity;
        s.ratio[index]=p[i].intensity / p[j].intensity; //(p.at(i).intensity+p.at(j).intensity);
        if(pepID[i]>-1){
          strcpy(s.peptide,mp[pepID[i]].sequence_Long);
          strcpy(s.protein,mp[pepID[i]].gene);
        } else if(pepID[j]>-1){
          strcpy(s.peptide,mp[pepID[j]].sequence_Long);
          strcpy(s.protein,mp[pepID[j]].gene);
        } else {
          strcpy(s.peptide,"");
          strcpy(s.protein,"");
        }
        count++;
        if(!checkSILAC(s,index))hits.push_back(s);        

      }

      dif = p[j].monoMass-8.0142036-p[i].monoMass;
      ppm = dif/p[i].monoMass*1000000;
      if(fabs(ppm)<5.0) {

        tp.clear();
        tp.add(p[i]);
        tp.add(p[j]);

        tp.pearson(0,1,true,true,s.rval[index],s.pval[index],slope,intercept,s.ratio[index]);

        s.charge=p[i].charge;
        s.monoMass=p[i].monoMass;
        s.ppm=ppm;
        s.da=p[j].monoMass-p[i].monoMass;
        s.rTime=p[i].rTime;
        s.intensity[index]=p[i].intensity;
        //s.ratio[index]= p[i].intensity / p[j].intensity; //(p.at(i).intensity+p.at(j).intensity);
        if(pepID[i]>-1){
          strcpy(s.peptide,mp[pepID[i]].sequence_Long);
          strcpy(s.protein,mp[pepID[i]].gene);
        } else if(pepID[j]>-1){
          strcpy(s.peptide,mp[pepID[j]].sequence_Long);
          strcpy(s.protein,mp[pepID[j]].gene);
        } else {
          strcpy(s.peptide,"");
          strcpy(s.protein,"");
        }
        count++;
        if(!checkSILAC(s,index))hits.push_back(s);

        if(p[i].datapoints>9){
          avgInc+=( p[j].intensity / (p[i].intensity+p[j].intensity));
          avgCount++;
        } 
      }

      if(ppm>100.0) break;
    }
  }

  cout << "Total pairs: " << count << endl;
  cout << "Average incorporation (10+ datapoints): " << avgInc/avgCount << endl;

}
*/