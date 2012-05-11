#include "CSILACtor.h"

//-------------------------------------
//   Constructors and destructors
//-------------------------------------
CSILACtor::CSILACtor(){
  replicates=0;
  xVal=0;
	corrThreshold=0.9;
	bVerbose=false;
	peps = new vector<singlet>;
  proteins = new vector<silacProtein>;
	pairs = new vector<silac>;
  labels = new vector<double>;
};
CSILACtor::CSILACtor(const CSILACtor& c){
  unsigned int i;
	peps = new vector<singlet>;
  proteins = new vector<silacProtein>;
	pairs = new vector<silac>;
  labels = new vector<double>;
  for(i=0;i<c.peps->size();i++) peps->push_back(c.peps->at(i));
  for(i=0;i<c.pairs->size();i++) pairs->push_back(c.pairs->at(i));
  for(i=0;i<c.labels->size();i++) labels->push_back(c.labels->at(i));
  for(i=0;i<c.proteins->size();i++) proteins->push_back(c.proteins->at(i));
  replicates=c.replicates;
	corrThreshold=c.corrThreshold;
  xVal=c.xVal;
	bVerbose=c.bVerbose;
}
CSILACtor::~CSILACtor(){
	if(peps) delete peps;
  if(proteins) delete proteins;
	if(pairs) delete pairs;
  if(labels) delete labels;
};

//-------------------------------------
//   Operator overrides
//-------------------------------------
CSILACtor& CSILACtor::operator=(const CSILACtor& c){
  if(this!=&c){
    unsigned int i;
    delete peps;
    delete proteins;
    delete pairs;
		delete labels;
		peps = new vector<singlet>;
		proteins = new vector<silacProtein>;
		pairs = new vector<silac>;
		labels = new vector<double>;
    for(i=0;i<c.peps->size();i++) peps->push_back(c.peps->at(i));
    for(i=0;i<c.pairs->size();i++) pairs->push_back(c.pairs->at(i));
    for(i=0;i<c.labels->size();i++) labels->push_back(c.labels->at(i));
    for(i=0;i<c.proteins->size();i++) proteins->push_back(c.proteins->at(i));
    replicates=c.replicates;
		corrThreshold=c.corrThreshold;
    xVal=c.xVal;
		bVerbose=c.bVerbose;
  }
  return *this;
}
silac& CSILACtor::operator[ ](const unsigned int i){
  return pairs->at(i);
}

//-------------------------------------
//   Data accessors
//-------------------------------------
void CSILACtor::clear(){
  delete peps;
  delete pairs;
  delete proteins;
	peps = new vector<singlet>;
	proteins = new vector<silacProtein>;
	pairs = new vector<silac>;
}
unsigned int CSILACtor::size(){
  return pairs->size();
}

//-------------------------------------
//   Data modifiers
//-------------------------------------
void CSILACtor::addLabel(double d){
  labels->push_back(d);
}
void CSILACtor::setCorrThreshold(double d){
	corrThreshold=d;
}


//-------------------------------------
//   Analysis Functions
//-------------------------------------
void CSILACtor::analyze(CKronik2& kro, float minRT, float maxRT, vector<uniquePep>& v){

  unsigned int i;
	int j,k;
  int count=0;
  double dif;
  double ppm;
	double maxLabel;

	double rval;
	double pval;

  int lys,arg;

  double avgInc=0.0;
  int avgCount=0;
  silac s;
  singlet p;

  double slope;
  double intercept;

  float area1,area2;

  //Filter Kronik2 results
  kro.filterRT(minRT,maxRT);
  cout << kro.size() << " from " << minRT << " to " << maxRT << endl;
  kro.removeContaminants(20.0);
  cout << kro.size() << " after removing contaminants." << endl;

  //Sort for speed
  kro.sortMonoMass();

	//Find maximum mass difference needed;
	maxLabel=0.0;
	for(i=0;i<labels->size();i++){
		if(labels->at(i)>maxLabel) maxLabel=labels->at(i);
	}
	maxLabel+=0.5; //Add some extra for buffer

  //All Kronik2 results are stored as single peptides
  for(i=0;i<kro.size();i++){
    p.area->clear();
    p.charge=kro[i].charge;
    p.area->push_back(kro[i].sumIntensity);
    p.monoMass=kro[i].monoMass;
    p.paired=0;
    p.peptide=kro[i].sequence;
    p.protein=kro[i].gene;
    p.rTime=kro[i].rTime;
    p.mz=kro[i].basePeak;
    peps->push_back(p);
  }

  //Iterate through all Kronik2 IDs to find pairs
	bool bMatch=false;
  for(i=0;i<kro.size();i++){
		bMatch=false;
    for(j=i+1;j<kro.size();j++){

			if( (kro[j].monoMass-kro[i].monoMass) > maxLabel) break;

      s.pval=1.0;
      s.rval=0.0;
      
      if(i==j) continue;
      if(j<0) continue;

      if(kro[i].charge!=kro[j].charge) continue;
      if(fabs(kro[i].rTime - kro[j].rTime) > 0.50) continue;

      for(k=0;k<labels->size();k++){
        dif = kro[j].monoMass-labels->at(k)-kro[i].monoMass;
      
        ppm = dif/kro[i].monoMass*1000000;
        if(fabs(ppm)<5.0) {
   
          kro.pearson(i,j,true,true,rval,pval,slope,intercept,area1,area2);
					s.rval=(float)rval;
					s.pval=(float)pval;

					//if(strcmp(kro[i].sequence,"")!=0 || strcmp(kro[j].sequence,"")!=0) {
					//	cout << kro[i].basePeak << "\t" << kro[i].charge << "\t" << kro[i].rTime << "\t" << kro[i].sequence << " matches " << kro[j].sequence << "\t" << ppm << "\t" << s.rval << "\t" << i << " - " << j << endl;
					//	bMatch=true;
					//}

					if(s.rval<(float)corrThreshold) {
						if(bVerbose){
							if(strcmp(kro[i].sequence,"")!=0 || strcmp(kro[j].sequence,"")!=0) {
								cout << "FAIL Corr: " << s.rval << "\t";
								cout << kro[i].monoMass << " " << kro[i].basePeak << " " << kro[i].charge << "" << kro[i].lowScan << "-" << kro[i].highScan << " " << kro[i].sequence;
								cout << " matches ";
								cout << kro[j].monoMass << " " << kro[j].basePeak << " " << kro[j].charge << "" << kro[j].lowScan << "-" << kro[j].highScan << " " << kro[j].sequence;
								cout << "\tPPM: " << ppm << endl;
							}
						}
						continue;
					}
          
          s.area1->clear();
          s.area2->clear();
          s.area1->push_back(area1);
          s.area2->push_back(area2);
          s.label=k;
          s.charge=kro[i].charge;
          s.monoMass=kro[i].monoMass;
          s.ppm=ppm;
          s.da=kro[j].monoMass-kro[i].monoMass;
          s.rTime=kro[i].rTime;
          s.pepIndex1=i;
          s.pepIndex2=j;
          s.mz=kro[i].basePeak;
          
          //Match Sequence
          
          if(strcmp(kro[i].sequence,"")!=0 && strcmp(kro[j].sequence,"")!=0){
						if(!matchSeq(kro[i].sequence,kro[j].sequence)) {
							//cout << "FAIL match: " << kro[i].sequence << " " << kro[j].sequence << endl;
							continue;
						}
          }

          //Do not count pair if light half has a modification (if it has a sequence)
          //or if the modifiable amino acids do not equal mod mass
          if(strcmp(kro[i].sequence,"")!=0){           
						if(countMods(kro[i].sequence)>0) {
							//cout << "FAIL match 2" << endl;
							continue;
						}
						if(strcmp(kro[j].sequence,"")!=0) {
							countTryptic(kro[j].sequence,lys,arg);
							ppm=lys*8.0142036+arg*6.0201324;
							if(fabs(ppm-labels->at(k))>1.0) {
								//cout << "FAIL match 3" << endl;
								continue;
							}
						}
          } 

          //Do not count pair if heavy half's modifications do not equal labels (if it has a sequence)
          if(strcmp(kro[j].sequence,"")!=0){  
            countTryptic(kro[j].sequence,lys,arg);
            ppm=lys*8.0142036+arg*6.0201324;
						if(fabs(ppm-labels->at(k))>1.0) {
							//cout << "FAIL match 4" << endl;
							continue;
						}
          }
          
          count++;

					
          peps->at(i).paired=1;
          peps->at(j).paired=1;
					//if(strcmp(kro[i].sequence,"")!=0 || strcmp(kro[j].sequence,"")!=0) cout << &peps->at(i).peptide[0] << " " << peps->at(i).paired << "\t\t" << &peps->at(j).peptide[0] << " " << peps->at(j).paired << endl;
					if(!checkSILAC(s,true)) {
						//if(strcmp(kro[i].sequence,"")!=0 || strcmp(kro[j].sequence,"")!=0) cout << "New pair" << endl;
						pairs->push_back(s);        
					//} else {
						//if(strcmp(kro[i].sequence,"")!=0 || strcmp(kro[j].sequence,"")!=0) cout << "Old pair" << endl;
					}
        }
      }
    }
		//if(!bMatch && strcmp(kro[i].sequence,"")!=0) cout << kro[i].sequence << "\t" << kro[i].charge << "\t" << kro[i].basePeak << "\t" << kro[i].monoMass << "\t" << kro[i].bestScan << " has no matches" << endl;
  }

	vector<int> vIndex;
	for(k=0;k<v.size();k++) vIndex.push_back(0);

  cout << "Total pairs: " << count << " from " << kro.size() << " peptide isotope distributions." << endl;
  j=0;
	cout << endl;
	for(i=0;i<kro.size();i++) {
		if(peps->at(i).paired==1) j++;
		if(strcmp(&peps->at(i).peptide[0],"")!=0 && peps->at(i).paired==1) {
			for(k=0;k<v.size();k++){
				if(v[k].seq.compare(peps->at(i).peptide)==0 && v[k].charge==peps->at(i).charge) {
					vIndex[k]=1;
					break;
				}
			}
		}
		if(bVerbose){
			if(strcmp(&peps->at(i).peptide[0],"")!=0 && peps->at(i).paired==0) cout << "No match: " << &peps->at(i).peptide[0] << " mass: " << peps->at(i).monoMass << " charge: " << peps->at(i).charge << " rtime: " << peps->at(i).rTime << endl;
		}
	}
  cout << "\n" << j << " PIDs are paired. " << (double)j/kro.size()*100.0 << "%" << endl;
	j=0;
	for(k=0;k<vIndex.size();k++) j+=vIndex[k];
	cout << j << " of " << vIndex.size() << " matched peptide sequences are paired." << endl;

  //FILE* fp=fopen("pairs.txt","at");
  //for(i=0;i<pairs->size();i++) fprintf(fp,"%.1f\n",kro[pairs->at(i).pepIndex1].intensity);
  //fprintf(fp,"\n\n");
  //fclose(fp);
  
  //cout << "Average incorporation (10+ datapoints): " << avgInc/avgCount << endl;

  double avg=0.0;
  double std=0.0;
  dif=0.0;
  for(i=0;i<pairs->size();i++){
    avg+=pairs->at(i).ppm;
    //cout << pairs->at(i).monoMass << "\t" << pairs->at(i).ratio << " (" << pairs->at(i).rval << ")\t" << pairs->at(i).peptide << "\t" << pairs->at(i).protein << "\t" << kro[pairs->at(i).pepIndex1].bestScan << "\t" << kro[pairs->at(i).pepIndex1].basePeak << endl;
  }
  avg/=pairs->size();
  for(i=0;i<pairs->size();i++) dif+=pow((pairs->at(i).ppm-avg),2);
  std=sqrt(dif/pairs->size());
  cout << "\nGlobal mass accuracy of pairs (ppm): " << avg << " +/- " << std << endl;

}

bool CSILACtor::checkSILAC(silac& s, bool update){
  int i;
  double ppm;

  //Check if peptide already found
  for(i=0;i<pairs->size();i++){
    ppm=(s.monoMass-pairs->at(i).monoMass)/pairs->at(i).monoMass*1000000;
    if(ppm>20.0) break;
    if(fabs(ppm)<10.0 && fabs(s.rTime-pairs->at(i).rTime)<5.0 && s.charge==pairs->at(i).charge && pairs->at(i).label==s.label){

      if(pairs->at(i).area1->at(0)>0 && s.area1->at(0)>pairs->at(i).area1->at(0) && update){
        pairs->at(i).rTime=s.rTime;
        pairs->at(i).ppm=s.ppm;
        pairs->at(i).da=s.da;
        pairs->at(i).monoMass=s.monoMass;
        pairs->at(i).area1->at(0)=s.area1->at(0);
        pairs->at(i).area2->at(0)=s.area2->at(0);
        pairs->at(i).pepIndex1=s.pepIndex1;
        pairs->at(i).pepIndex2=s.pepIndex2;
        pairs->at(i).corr=s.corr;
        pairs->at(i).label=s.label;
        pairs->at(i).slope=s.slope;
        pairs->at(i).intercept=s.intercept;
        pairs->at(i).pval=s.pval;
        pairs->at(i).rval=s.rval;
        pairs->at(i).mz=s.mz;
      }
      return true;
    }
  }

  return false;
}

void CSILACtor::countTryptic(char* seq, int& k, int& r){
  int i;
  k=0;
  r=0;
  if(strlen(seq)<2) return;
  for(i=1;i<strlen(seq)-1;i++){
    if(seq[i]=='K' || seq[i]=='k') k++;
    if(seq[i]=='R' || seq[i]=='r') r++;
  }
}

int CSILACtor::countMods(char* seq){
  int i,j;
  j=0;
  for(i=0;i<strlen(seq);i++){
		if(seq[i]=='#' || seq[i]=='%') j++;
    //if(!isalpha(seq[i]) && seq[i]!='.' && seq[i]!='-') j++;
  }
  return j;
}

bool CSILACtor::matchSeq(char* seq1,char* seq2){
  int i=0;
  int j=0;

  while(i<strlen(seq1) && j<strlen(seq2)){
    if(!isalpha(seq1[i])){
      i++;
      continue;
    }
    if(!isalpha(seq2[j])){
      j++;
      continue;
    }
    if(seq1[i]!=seq2[j]) return false;
    i++;
    j++;
  }
  return true;
}

void CSILACtor::addReplicate(CSILACtor& c){
  int i,j;
  int n=pairs->size();
  bool bMatch;

	int pos,lower,upper;

  double ppm;

  //increment the replicates. Indicates depth of each area value
  replicates++;

  //if(xVal!=c.xVal) cout << "Warning!!! Combining technical replicates with different independent values!" << endl;

  //set up temporary singlet
  singlet s;
  s.area->push_back(0.0);
  for(i=0;i<replicates;i++) s.area->push_back(0.0);

  //Add singlets first
  //Increment existing singles area values
  for(i=0;i<peps->size();i++) peps->at(i).area->push_back(0.0);

	//Sort by monoMass for binary search
	c.sortSinglets();
	sortSinglets();

  //Match new singlets with old singlets
  vector<int> newIndex;
  for(i=0;i<c.peps->size();i++){

    //Check against existing singlets for match
    bMatch=false;

		//Find closest match within tolerance, then step both ways to find boundaries
		pos=binarySearch(peps,c.peps->at(i).monoMass);
		if(pos>-1) {
			lower=pos;
			upper=pos;
			for(j=pos-1;j>-1;j--){
				ppm = (peps->at(j).monoMass-c.peps->at(i).monoMass)/c.peps->at(i).monoMass * 1000000;
				if(ppm>-10.0) lower=j;
				else break;
			}
			for(j=pos+1;j<peps->size();j++){
				ppm = (peps->at(j).monoMass-c.peps->at(i).monoMass)/c.peps->at(i).monoMass * 1000000;
				if(ppm<10.0) upper=j;
				else break;
			}
		} else {
			lower=0;
			upper=-1;
		}

    for(j=lower;j<=upper;j++){
      
      //ppm = (c.peps->at(i).monoMass-peps->at(j).monoMass)/peps->at(j).monoMass*1000000;

      //if we have a match
			//we are also already at ppm tolerance
      if(/*fabs(ppm)<10.0 &&*/ fabs(c.peps->at(i).rTime-peps->at(j).rTime)<5.0 && c.peps->at(i).charge==peps->at(j).charge){

        bMatch=true;
        
        //Mark new index for this peptide
        newIndex.push_back(j);
        
        //set appropriate area
        peps->at(j).area->at(replicates)=c.peps->at(i).area->at(0);

        break;
      }
    }

    //if we do not have a match
    if(!bMatch){

      //set up temporary singlet
      s.area->at(replicates)=c.peps->at(i).area->at(0);
      s.charge=c.peps->at(i).charge;
      s.monoMass=c.peps->at(i).monoMass;
      s.paired=c.peps->at(i).paired;
      s.peptide=c.peps->at(i).peptide;
      s.protein=c.peps->at(i).protein;
      s.rTime=c.peps->at(i).rTime;
      s.mz=c.peps->at(i).mz;

      //add and mark new index
      peps->push_back(s);
      newIndex.push_back(peps->size()-1);
    }
  }

  //Now add pairs
  cout << "Adding: " << c.size() << " to: " << pairs->size() << endl;

	//Sort for binary search
	c.sortPairs();
	sortPairs();

  //Go through all pairs
  for(i=0;i<c.size();i++){

		bMatch=false;

		//Find closest match within tolerance, then step both ways to find boundaries
		pos=binarySearch(pairs,c[i].monoMass);
		if(pos>-1){
			lower=pos;
			upper=pos;
			for(j=pos-1;j>-1;j--){
				ppm = (pairs->at(j).monoMass-c[i].monoMass)/c[i].monoMass * 1000000;
				if(ppm>-10.0) lower=j;
				else break;
			}
			for(j=pos+1;j<pairs->size();j++){
				ppm = (pairs->at(j).monoMass-c[i].monoMass)/c[i].monoMass * 1000000;
				if(ppm<10.0) upper=j;
				else break;
			}
		} else {
			lower=0;
			upper=-1;
		}
    
    for(j=lower;j<=upper;j++){
      //ppm=(c[i].monoMass-pairs->at(j).monoMass)/pairs->at(j).monoMass*1000000;

      //if there is a match
      if(/*fabs(ppm)<10.0 && */ fabs(c[i].rTime-pairs->at(j).rTime)<5.0 && c[i].charge==pairs->at(j).charge && pairs->at(j).label==c[i].label){

        //only increment the areas
        pairs->at(j).area1->push_back(c[i].area1->at(0));
        pairs->at(j).area2->push_back(c[i].area2->at(0));

        bMatch=true;
        break;
      }
      if(bMatch) break;
    }

    //if there was no match
    if(!bMatch) {
      
      //Add the pair with the new peptide indexes
      c[i].pepIndex1=newIndex[c[i].pepIndex1];
      c[i].pepIndex2=newIndex[c[i].pepIndex2];
      pairs->push_back(c[i]);

      /*
      for(j=0;j<peps->size();j++) {
        if(c.peps[c[i].pepIndex1].protein.compare(peps->at(j).protein)==0 && c.peps[c[i].pepIndex1].peptide.compare(peps->at(j).peptide)==0 ){
          c[i].pepIndex1=j;
          break;
        }
      }

      for(j=0;j<peps->size();j++) {
        if(c.peps[c[i].pepIndex2].protein.compare(peps->at(j).protein)==0 && c.peps[c[i].pepIndex2].peptide.compare(peps->at(j).peptide)==0 ){
          c[i].pepIndex2=j;
          break;
        }
      }
      pairs->push_back(c[i]);
      */
    
    }

  }

  cout << "After: " << pairs->size() << endl;

	//clear out replicates (hopefully we don't need them) to conserve memory:
	c.clear();

}

void CSILACtor::proteinSummary(char* out, char* target, bool bTrace){
  string protein;
  string peptide;
  string peptide2;
  string s;

  silacProtein p;
	pepSILAC ps;
	vector<pepSILAC> vps;

  proteins->clear();

  float avg;
  float std;
  float symRatio;
  float symStd;
  float fracRatio;
  float fracStd;

  float dif;

  vector<float> symRatios;
  vector<float> symStds;
  vector<float> fracRatios;
  vector<float> fracStds;
  vector<float> avgs;
  vector<float> stds;
  vector<string> pepper;
  vector<float> knownPepRatios;

  vector<int> index;

  MRT m;
  vector<MRT> list;

  bool bMatch;

  int i,j,k;

  FILE* f;
  FILE* f2=fopen("trace.txt","wt");   

  //Open file for output
  f=fopen(out,"wt");
  if(f==NULL){
    cout << "Cannot open " << out << " to export results. Stopping analysis." << endl;
    return;
  }

  //Output header line
  //fprintf(f,"Protein\tPeptide\tCharge\tRTime\tMass\tAbundance\tTrue Ratio\t+/-\tSymmetric Ratio\t+/-\tFractional Ratio\t+/-\tR-value\tDifference\n");
	fprintf(f,"Protein\tRIA\tStDev\tPeptide\tCharge\tLight Abundance (avg)\tHeavy Abundance (avg)\tPeptide RIA\tPeptide StDev\tR-value\tOutlier\n");

  //Check all peptide pairs
  for(i=0;i<pairs->size()-1;i++){
   
    //check for protein id
    if(peps->at(pairs->at(i).pepIndex1).protein.compare("")==0 && peps->at(pairs->at(i).pepIndex2).protein.compare("")==0) {
     
      //if there are no IDs for either peptide, make targeted list IF ratio is interesting...
      //scorePair(i,avg,std,symRatio,symStd,fracRatio,fracStd);
      //if(fabs

      continue;
    }

    //store protein ID
    if(peps->at(pairs->at(i).pepIndex1).protein.compare("")!=0) {
      protein=peps->at(pairs->at(i).pepIndex1).protein;
      peptide=peps->at(pairs->at(i).pepIndex1).peptide;
    } else {
      protein=peps->at(pairs->at(i).pepIndex2).protein;
      peptide=peps->at(pairs->at(i).pepIndex2).peptide;
    }

    //check if protein already analyzed
    bMatch=false;
    for(j=0;j<proteins->size();j++){
      if(protein.compare(proteins->at(j).protein)==0) {
        bMatch=true;
        break;
      }
    }
    if(bMatch) continue;

    //add protein to checklist
    p.protein=protein;

    symRatios.clear();
    fracRatios.clear();
    symStds.clear();
    fracStds.clear();
    avgs.clear();
    stds.clear();
    index.clear();
    pepper.clear();

    //Get average
    scorePair(i,avg,std,symRatio,symStd,fracRatio,fracStd);

    //fprintf(f,"%s\t%s\t%d\t%.2f\t%.4lf\t%.0f\t%.2f\t%.4f\t%.2f\t%.4f\t%.2f\t%.4f\t%lf\t%lf\n",&protein[0],&peptide[0],pairs->at(i).charge,pairs->at(i).rTime,pairs->at(i).monoMass,pairs->at(i).area1[0]+pairs->at(i).area2[0],avg,std,symRatio,symStd,fracRatio,fracStd,pairs->at(i).rval,pairs->at(i).da);
    if(symRatio>0) symRatios.push_back(symRatio-1.0f);
    else symRatios.push_back(symRatio+1.0f);
    symStds.push_back(symStd);
    fracStds.push_back(fracStd);
    fracRatios.push_back(fracRatio);
    avgs.push_back(avg);
    stds.push_back(std);
    index.push_back(i);
    pepper.push_back(peptide);

    if(bTrace){
      fprintf(f2,"\n%s\t%s\n",&protein[0],&peptide[0]);
      for(j=0;j<=replicates;j++){
        fprintf(f2,"%d\t%.0f\t%.0f\n",j,peps->at(pairs->at(i).pepIndex1).area[j],peps->at(pairs->at(i).pepIndex2).area[j]);
      }
    }

    for(j=i+1;j<pairs->size();j++){
      //check for protein id
      if(protein.compare(peps->at(pairs->at(j).pepIndex1).protein)!=0 && protein.compare(peps->at(pairs->at(j).pepIndex2).protein)!=0) continue;
      
      if(peps->at(pairs->at(j).pepIndex1).protein.compare("")!=0) peptide2=peps->at(pairs->at(j).pepIndex1).peptide;
      else peptide2=peps->at(pairs->at(j).pepIndex2).peptide;

      scorePair(j,avg,std,symRatio,symStd,fracRatio,fracStd);
      //fprintf(f,"%s\t%s\t%d\t%.2f\t%.4lf\t%.0f\t%.2f\t%.4f\t%.2f\t%.4f\t%.2f\t%.4f\t%lf\t%lf\n",&protein[0],&peptide2[0],pairs->at(j).charge,pairs->at(j).rTime,pairs->at(j).monoMass,pairs->at(j).area1[0]+pairs->at(j).area2[0],avg,std,symRatio,symStd,fracRatio,fracStd,pairs->at(j).rval,pairs->at(j).da);
      if(symRatio>0) symRatios.push_back(symRatio-1.0f);
      else symRatios.push_back(symRatio+1.0f);
      fracRatios.push_back(fracRatio);
      symStds.push_back(symStd);
      fracStds.push_back(fracStd);
      avgs.push_back(avg);
      stds.push_back(std);
      index.push_back(j);
      pepper.push_back(peptide2);

      if(bTrace){
        fprintf(f2,"\n%s\t%s\n",&protein[0],&peptide2[0]);
        for(k=0;k<=replicates;k++){
          fprintf(f2,"%d\t%.0f\t%.0f\n",k,peps->at(pairs->at(j).pepIndex1).area[k],peps->at(pairs->at(j).pepIndex2).area[k]);
        }
      }

    }

    avg=0.0;
    for(j=0;j<(int)fracRatios.size();j++) avg+=fracRatios[j];
    avg/=fracRatios.size();

    dif=0.0;
    for(j=0;j<(int)fracRatios.size();j++) dif+=pow(fracRatios[j]-avg,2);
    std=sqrt(dif/fracRatios.size());

    //erase outliers
		vps.clear();
    while(true){
      dif=0.0;
      j=-1;
      for(k=0;k<fracRatios.size();k++){

        if(fabs(fracRatios[k]-avg)>2*std){
          if(j==-1){
            dif=fabs(fracRatios[k]-avg);
            j=k;
          } else if(fabs(fracRatios[k]-avg)>dif) {
            dif=fabs(fracRatios[k]-avg);
            j=k;
          }
        }
      }
      if(j==-1) break;

			ps.charge=pairs->at(index[j]).charge;
			ps.lIntensity=pairs->at(index[j]).area1->at(0);
			ps.hIntensity=pairs->at(index[j]).area2->at(0);
			strncpy(ps.peptide,&pepper[j][0],63);
			ps.ratio=fracRatios[j];
			ps.rTime=pairs->at(index[j]).rTime;
			ps.rValue=pairs->at(index[j]).rval;
			ps.std=fracStds[j];
			vps.push_back(ps);

      fracRatios.erase(fracRatios.begin()+j);
      symRatios.erase(symRatios.begin()+j);
      fracStds.erase(fracStds.begin()+j);
      symStds.erase(symStds.begin()+j);
      avgs.erase(avgs.begin()+j);
      stds.erase(stds.begin()+j);
      index.erase(index.begin()+j);
      pepper.erase(pepper.begin()+j);

      avg=0.0;
      for(j=0;j<fracRatios.size();j++) avg+=fracRatios[j];
      avg/=fracRatios.size();

      dif=0.0;
      for(j=0;j<fracRatios.size();j++) dif+=pow(fracRatios[j]-avg,2);
      std=sqrt(dif/fracRatios.size());
    }
    for(j=0;j<fracRatios.size();j++) knownPepRatios.push_back(fracRatios[j]);

    //output results
    for(j=0;j<index.size();j++){

			//report average heavy and light intensity
			ps.lIntensity=0.0;
			ps.hIntensity=0.0;
			for(k=0;k<pairs->at(index[j]).area1->size();k++){
	      ps.lIntensity=pairs->at(index[j]).area1->at(k);
				ps.hIntensity=pairs->at(index[j]).area2->at(k);
			}
			if(k>0) ps.lIntensity/=k;
			if(k>0) ps.hIntensity/=k;

      //fprintf(f,"%s\t%s\t%d\t%.2f\t%.4lf\t%.0f\t%.2f\t%.4f\t%.2f\t%.4f\t%.2f\t%.4f\t%lf\t%lf\n",&protein[0],&pepper[j][0],pairs[index[j]].charge,pairs[index[j]].rTime,pairs[index[j]].monoMass,pairs[index[j]].area1[0]+pairs[index[j]].area2[0],avgs[j],stds[j],symRatio,symStds[j],fracRatios[j],fracStds[j],pairs[index[j]].rval,pairs[index[j]].da);
			fprintf(f,"%s\t%.4f\t%.4f\t%s\t%d\t%.0f\t%.0f\t%.4f\t%.4f\t%.4lf\t%s\n",&protein[0],avg,std,&pepper[j][0],pairs->at(index[j]).charge,ps.lIntensity,ps.hIntensity,fracRatios[j],fracStds[j],pairs->at(index[j]).rval,"");
    }
		for(j=0;j<vps.size();j++){
			fprintf(f,"%s\t%.4f\t%.4f\t%s\t%d\t%.0f\t%.0f\t%.4f\t%.4f\t%.4lf\t%s\n",&protein[0],avg,std,vps[j].peptide,vps[j].charge,vps[j].lIntensity,vps[j].hIntensity,vps[j].ratio,vps[j].std,vps[j].rValue,"*");
		}
  
    //fprintf(f,"Average: %lf +/- %lf\n\n",avg,std);

    p.ratio=avg;
    p.std=std;

    proteins->push_back(p);

  }

  fclose(f);
  fclose(f2);

  //build targeted mass list
  double avg2=0.0;
	//FILE *xfxf=fopen("overallDist.txt","wt");
  for(j=0;j<pairs->size();j++) {
    scorePair(j,avg,std,symRatio,symStd,fracRatio,fracStd);
    avg2+=fracRatio;
		//fprintf(xfxf,"%.2lf\n",fracRatio);
  }
  avg2/=pairs->size();
	//fclose(xfxf);

  double dif2=0.0;
  for(j=0;j<pairs->size();j++) {
    scorePair(j,avg,std,symRatio,symStd,fracRatio,fracStd);
    dif2+=pow(fracRatio-avg2,2);
  }
  double std2=sqrt(dif2/pairs->size());

  if(strlen(target)<2) return;

	cout << "Target list is below: " << avg2-2*std2 << " and above: " << avg2+2*std2 << endl;
	cout << "Avg: " << avg2 << " Std: " << std2 << endl;
	cout << target << endl;

  for(i=0;i<pairs->size()-1;i++){
   
    //check for protein id
    if(peps->at(pairs->at(i).pepIndex1).protein.compare("")!=0 || peps->at(pairs->at(i).pepIndex2).protein.compare("")!=0) continue;

    scorePair(i,avg,std,symRatio,symStd,fracRatio,fracStd);

    //if we're a standard deviation away...
    if( (fracRatio>(avg2+2*std2) || fracRatio<(avg2-2*std2)) ){
      m.mz=pairs->at(i).mz;
      m.rt1=pairs->at(i).rTime-2.0f;
      m.rt2=pairs->at(i).rTime+2.0f;
      list.push_back(m);
			//cout << pairs->at(i).mz << "\t" << pairs->at(i).da << "\t" << pairs->at(i).charge << endl;
			m.mz=pairs->at(i).mz+pairs->at(i).da/pairs->at(i).charge;
      m.rt1=pairs->at(i).rTime-2.0f;
      m.rt2=pairs->at(i).rTime+2.0f;
      list.push_back(m);
    }
  }

  //remove redundancies
  int count=0;
  while(list.size()!=count){
    count=list.size();
    for(i=0;i<list.size()-1;i++){
      for(j=i+1;j<list.size();j++){
        if(fabs(list[i].mz-list[j].mz)<0.02 &&
           ((list[j].rt1 > list[i].rt1-0.02 && list[j].rt1 < list[i].rt2+0.02) ||
           (list[j].rt2 > list[i].rt1-0.02 && list[j].rt2 < list[i].rt2+0.02)) ){
          if(list[j].rt1<list[i].rt1) list[i].rt1=list[j].rt1;
          if(list[j].rt2>list[i].rt2) list[i].rt2=list[j].rt2;
          list.erase(list.begin()+j);
          j--;
        }
      }
    }
  }

  //Export list
  FILE* f3=fopen(target,"wt");
  if(f3==NULL){
    cout << "Cannot export mass and time list to " << target << endl;
    return;
  }
  for(i=0;i<list.size();i++){
    fprintf(f3,"%.4lf\t%.4f\t%.4f\n",list[i].mz,list[i].rt1,list[i].rt2);
  }
  fclose(f3);

}

void CSILACtor::scorePair(int index, float& ratio, float& std, float& symRatio, float& symStd, float& fracRatio, float& fracStd){
  
  unsigned int j;
  float dif;

  //Get ratio
  ratio=0.0;
  for(j=0;j<pairs->at(index).area1->size();j++){
    ratio += (pairs->at(index).area1->at(j)/pairs->at(index).area2->at(j));
  }
  ratio/=pairs->at(index).area1->size();

  //Get standard deviation
  dif=0.0;
  for(j=0;j<pairs->at(index).area1->size();j++){
    dif += pow(( (pairs->at(index).area1->at(j)/pairs->at(index).area2->at(j))-ratio),2);
  }
  std=sqrt(dif/(float)pairs->at(index).area1->size());

  //Get symmetrical ratio
  symRatio=0.0;
  for(j=0;j<pairs->at(index).area1->size();j++){
    if(pairs->at(index).area1->at(j)>pairs->at(index).area2->at(j)) symRatio += (pairs->at(index).area1->at(j)/pairs->at(index).area2->at(j)-1.0f);
    else symRatio += (-pairs->at(index).area2->at(j)/pairs->at(index).area1->at(j)+1.0f);
  }
  symRatio/=pairs->at(index).area1->size();

  //Get symmetrical standard deviation
  dif=0.0;
  for(j=0;j<pairs->at(index).area1->size();j++){
    if(pairs->at(index).area1->at(j)>pairs->at(index).area2->at(j)) dif += pow(( (pairs->at(index).area1->at(j)/pairs->at(index).area2->at(j)-1.0f)-symRatio),2);
    else dif += pow(( (-pairs->at(index).area2->at(j)/pairs->at(index).area1->at(j)+1.0f)-symRatio),2);
  }
  symStd=sqrt(dif/(float)pairs->at(index).area1->size());

  if(symRatio>0) symRatio+=1.0;
  else symRatio-=1.0;

  //Get fractional ratio
  fracRatio=0.0;
  for(j=0;j<pairs->at(index).area1->size();j++){
    fracRatio += (pairs->at(index).area1->at(j)/(pairs->at(index).area1->at(j)+pairs->at(index).area2->at(j)));
  }
  fracRatio/=pairs->at(index).area1->size();

  //Get standard deviation
  dif=0.0;
  for(j=0;j<pairs->at(index).area1->size();j++){
    dif += pow(( (pairs->at(index).area1->at(j)/(pairs->at(index).area1->at(j)+pairs->at(index).area2->at(j)))-fracRatio),2);
  }
  fracStd=sqrt(dif/(float)pairs->at(index).area1->size());
}

void CSILACtor::reduceAreas(){
  int i,j,k;
  float avg;
  cout << "Reducing" << endl;
  for(i=0;i<peps->size();i++){

    avg=0.0f;
    k=0;
    for(j=0;j<peps->at(i).area->size();j++) {
      if(peps->at(i).area->at(j)>0){
        avg+=peps->at(i).area->at(j);
        k++;
      }
    }
    avg/=k;

    peps->at(i).area->clear();
    peps->at(i).area->push_back(avg);
  }
  replicates=0;
  cout << "ok" << endl;
}

void CSILACtor::sortPairs(){
	qsort(&pairs->at(0),pairs->size(),sizeof(silac),comparePairs);
}

void CSILACtor::sortSinglets(){
	qsort(&peps->at(0),peps->size(),sizeof(singlet),compareSinglets);
}

int CSILACtor::comparePairs(const void *p1, const void *p2){
  const silac d1 = *(silac *)p1;
  const silac d2 = *(silac *)p2;
	if(d1.monoMass<d2.monoMass) return -1;
	else if(d1.monoMass>d2.monoMass) return 1;
  else return 0;
}

int CSILACtor::compareSinglets(const void *p1, const void *p2){
  const singlet d1 = *(singlet *)p1;
  const singlet d2 = *(singlet *)p2;
	if(d1.monoMass<d2.monoMass) return -1;
	else if(d1.monoMass>d2.monoMass) return 1;
  else return 0;
}

int CSILACtor::binarySearch(vector<silac>* p, double mass){
	int lower,mid,upper;
	int sz=p->size();
	double lowMass=mass-(mass/1000000*10.0);
	double highMass=mass+(mass/1000000*10.0);

	if(sz<1) return -1;

	//Find a mass within ppm tolerance. Use binary search to for speed.
	//Stop searching when first mass is found (although a second may exist).
	mid=sz/2;
	lower=0;
	upper=sz;
	
	while(p->at(mid).monoMass<lowMass || p->at(mid).monoMass>highMass){
		if(lower>=upper) break;
		if(mass<p->at(mid).monoMass){
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
	if(p->at(mid).monoMass<lowMass || p->at(mid).monoMass>highMass) return -1;
	return mid;
}

int CSILACtor::binarySearch(vector<singlet>* p, double mass){
	int lower,mid,upper;
	int sz=p->size();
	double lowMass=mass-(mass/1000000*10.0);
	double highMass=mass+(mass/1000000*10.0);

	if(sz<1) return -1;

	//Find a mass within ppm tolerance. Use binary search to for speed.
	//Stop searching when first mass is found (although a second may exist).
	mid=sz/2;
	lower=0;
	upper=sz;
	
	while(p->at(mid).monoMass<lowMass || p->at(mid).monoMass>highMass){
		if(lower>=upper) break;
		if(mass<p->at(mid).monoMass){
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
	if(p->at(mid).monoMass<lowMass || p->at(mid).monoMass>highMass) return -1;
	return mid;
}

