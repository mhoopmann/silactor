#include "MascotParser.h"

using namespace std;

MascotParser::MascotParser(){
  bDistiller=false;
  bRTime=false;
}

MascotParser::~MascotParser(){
}

MascotParser::MascotParser(const MascotParser& p){
  m.clear();
  fixMod.clear();
  varMod.clear();
	for(unsigned int i=0;i<p.m.size();i++) m.push_back(p.m[i]);
  for(unsigned int i=0;i<p.fixMod.size();i++) fixMod.push_back(p.fixMod[i]);
  for(unsigned int i=0;i<p.varMod.size();i++) varMod.push_back(p.varMod[i]);
  bDistiller=p.bDistiller;
  bRTime=p.bRTime;
}

MascotParser& MascotParser::operator=(const MascotParser& p){
  if(this!=&p){
    m.clear();
    fixMod.clear();
    varMod.clear();
	  for(unsigned int i=0;i<p.m.size();i++) m.push_back(p.m[i]);
    for(unsigned int i=0;i<p.fixMod.size();i++) fixMod.push_back(p.fixMod[i]);
    for(unsigned int i=0;i<p.varMod.size();i++) varMod.push_back(p.varMod[i]);
    bDistiller=p.bDistiller;
    bRTime=p.bRTime;
  } 
  return *this;
}

Mascot& MascotParser::operator[ ](const unsigned int i){
  return m[i];
}

void MascotParser::add(Mascot& p){
  m.push_back(p);
}

Mascot& MascotParser::at(unsigned int i){
	return m[i];
}

void MascotParser::clear(){
	m.clear();
  varMod.clear();
  fixMod.clear();
}

void MascotParser::erase(unsigned int i){
  m.erase(m.begin()+i);
}

void MascotParser::excel(char *fn){
	//FILE *f;
	//int i;

  printf("Excel output temporarily disabled. See Mike for help.\n");

  /*
	f=fopen(fn,"wt");
	if(f==NULL){
		cout << "Cannot open " << fn << endl;
		return;
	}

	fprintf(f,"Scan Number\tExact Monoisotopic Mass\tSEQUEST Monoisotopic Mass\tObserved Monoisotopic Mass\tCharge State\txCorr\tdeltaCN\tSpR\tSpScore\tIon%\tIntensity\t#\tSequence\tGene\tFile\n");
	for(i=0;i<DTA->size();i++){
		fprintf(f,"%d\t%lf\t%lf\t%lf\t%d\t%lf\t%lf\t%d\t%lf\t%lf\t%lf\t%d\t%s\t%s\t%s\n",DTA->at(i).scanNum,DTA->at(i).CalcMH_Exact,DTA->at(i).CalcMH_SQ,
			                                                                               DTA->at(i).ObsMH,DTA->at(i).charge,DTA->at(i).xCorr,
																																										 DTA->at(i).deltaCN,DTA->at(i).SpR,DTA->at(i).SpScore,
																																										 DTA->at(i).ionPercent,DTA->at(i).intensity,DTA->at(i).count,
																																										 DTA->at(i).sequence_Long,DTA->at(i).gene,DTA->at(i).file);
	}

	fclose(f);
  */
}

bool MascotParser::isDistiller(){
  return bDistiller;
}

bool MascotParser::hasRTime(){
  return bRTime;
}

void MascotParser::load(char *fn){
  Mascot d;
	FILE *f;
	size_t s;

	f=fopen(fn,"rb");
	if(f==NULL){
		cout << "Cannot open " << fn << endl;
		return;
	}

  clear();
	while(!feof(f)){
		s = fread(&d,sizeof(Mascot),1,f);
		m.push_back(d);
	}

	fclose(f);

}

void MascotParser::readText(char *fn){
  FILE *f;
  Mascot md;
  Mod mo;

  char str[1024];
  //char title[256];
  char *tok;
  char peptide[128];

  int i,j;
  char mData[100][256];  //wasteful, but quick and easy

  //Indexes for Mascot results
  int iCharge;
  int iScore;
  int iExpect;
  int iProtein;
  int iResBefore;
  int iSeq;
  int iResAfter;
  int iScanTitle;
  int iMods;
  bool bParen;

  int iState=0; //0=none, 1=fixed mods, 2=variable mods

  //staticC=true;

  f=fopen(fn,"rt");
  clear();
  
  //Read header, loading mods if needed
  while(!feof(f)){
    fgets(str,1024,f);

    tok=strtok(str,",\n");
    if(tok==NULL) continue;
    if(strcmp(tok,"\"Fixed modifications\"")==0 || strcmp(tok,"Fixed modifications")==0) iState=1;
    else if(strcmp(tok,"\"Variable modifications\"")==0 || strcmp(tok,"Variable modifications")==0) iState=2;
    else if(strcmp(tok,"prot_hit_num")==0) break;
    else if(strcmp(tok,"_DISTILLER_RAWFILE")==0) bDistiller=true;
  
    switch(iState){
      case 1:
        if(atoi(tok)>0) {
          tok=strtok(NULL,",\n");
          strcpy(peptide,tok);
          tok=strtok(NULL,",\n");
          mo.modMass=atof(tok);
          tok=strtok(peptide," \n");
          tok=strtok(NULL,"()\n");
          strcpy(mo.aa,tok);
          fixMod.push_back(mo);
        }
        break;

      case 2:
        if(atoi(tok)>0) {
          tok=strtok(NULL,",\n");
          tok=strtok(NULL,",\n");
          mo.modMass=atof(tok);
          if(varMod.size()==0) strcpy(mo.aa,"@");
          else if(varMod.size()==1) strcpy(mo.aa,"#");
          else if(varMod.size()==2) strcpy(mo.aa,"$");
          else if(varMod.size()==3) strcpy(mo.aa,"%");
          else if(varMod.size()==4) strcpy(mo.aa,"^");
          else if(varMod.size()==5) strcpy(mo.aa,"&");
          else cout << "WARNING: Too many variable mods." << endl;
          varMod.push_back(mo);
        }
        break;
      default:
        break;
    }

  }

  //get indexes for columns
  i=0;
  while(tok!=NULL){
    if(!strcmp(tok,"pep_score")) iScore=i;
    else if(!strcmp(tok,"prot_acc")) iProtein=i;
    else if(!strcmp(tok,"pep_expect")) iExpect=i;
    else if(!strcmp(tok,"pep_res_before")) iResBefore=i;
    else if(!strcmp(tok,"pep_seq")) iSeq=i;
    else if(!strcmp(tok,"pep_res_after")) iResAfter=i;
    else if(!strcmp(tok,"pep_scan_title")) iScanTitle=i;
    else if(!strcmp(tok,"pep_exp_z")) iCharge=i;
    else if(!strcmp(tok,"pep_var_mod_pos")) iMods=i;
    tok=strtok(NULL,",\n");
    i++;
  }


  //read PSMs until end of file
  while(!feof(f)){
    if(fgets(str,1024,f)==NULL) continue;

    //Tokenize entire string
    //This isn't easy because commas appear as parts of fields.
    bParen=false;
    i=0;
    strcpy(mData[i],"");
    for(j=0;j<(int)strlen(str);j++){
      if(str[j]==','){
        if(!bParen){
          i++;
          strcpy(mData[i],"");
        } else {
          strncat(mData[i],&str[j],1);
        }
      } else if(str[j]=='\"'){
        if(bParen) bParen=false;
        else bParen=true;
      } else {
        strncat(mData[i],&str[j],1);
      }
    }

    //populate the datafields

    //If there are variable mods, process them here
    if(strlen(mData[iMods])>0){
      j=0;
      for(i=0;i<(int)strlen(mData[iMods]);i++){
        if(i==0 || i==strlen(mData[iMods])-1){
          //do nothing; no sequence to add
        }else if(i==1 || i==strlen(mData[iMods]-2)){
          //do nothing for the dots
        } else { 
          md.sequence[j++]=mData[iSeq][i-2];
        }
        switch(mData[iMods][i]){
          case '1': md.sequence[j++]=varMod[0].aa[0]; break;
          case '2': md.sequence[j++]=varMod[1].aa[0]; break;
          case '3': md.sequence[j++]=varMod[2].aa[0]; break;
          case '4': md.sequence[j++]=varMod[3].aa[0]; break;
          case '5': md.sequence[j++]=varMod[4].aa[0]; break;
          case '6': md.sequence[j++]=varMod[5].aa[0]; break;
          default:  break;
        }
      }
      md.sequence[j]='\0';
    } else {
      strcpy(md.sequence,mData[iSeq]);
    }

    strcpy(peptide,mData[iResBefore]);    //The prior AA
    strcat(peptide,".");                  //Add a period
    strcat(peptide,md.sequence);          //The actual peptide sequence
    strcat(peptide,".");                  //Add a period
    strcat(peptide,mData[iResAfter]);     //The following AA
    strcpy(md.sequence_Long,peptide);

    md.mascotScore=(float)atof(mData[iScore]);
    md.charge=atoi(mData[iCharge]);
    strcpy(md.gene,mData[iProtein]);

    //get the scan number
    if(bDistiller){
      tok=strtok(mData[iScanTitle],": ()rt=\n");
      tok=strtok(NULL,": ()rt=\n");
      if(strcmp(tok,"Sum")==0){
        tok=strtok(NULL,": ()rt=\n");
        tok=strtok(NULL,": ()rt=\n");
        tok=strtok(NULL,": ()rt=\n");
        tok=strtok(NULL,": ()rt=\n");
        tok=strtok(NULL,": ()rt=\n");
      }
      tok=strtok(NULL,": ()rt=\n");
      md.scanNum=atoi(tok);
      tok=strtok(NULL,": ()rt=\n");
      md.rTime=(float)atof(tok);
      //cout << md.scanNum << " "<< md.rTime << endl;
    } else {
      if(strstr(mData[iScanTitle],"RT=")!=NULL){
        tok=strtok(mData[iScanTitle],".\n");
        tok=strtok(NULL,".\n");
        md.scanNum=atoi(tok);

        tok=strtok(NULL,"=\n");
        tok=strtok(NULL,"=\n");
        //cout << tok << endl;
        md.rTime=(float)(atof(tok)/60.0);

        bRTime=true;

      } else {
        tok=strtok(mData[iScanTitle],".\n");
        tok=strtok(NULL,".\n");
        md.scanNum=atoi(tok);
      }
    }

    md.zeroMass = calcMonoMass(md.sequence);
    //cout << md.sequence << " " << md.zeroMass << endl;
    m.push_back(md);

  }

  fclose(f);

}


void MascotParser::save(char *fn){

}

int MascotParser::size(){
	return m.size();
}


double MascotParser::calcMonoMass(char *seq){

	double mass=0.0;

	int H=0;
	int C=0;
	int N=0;
	int O=0;
	int S=0;
  int m[6];

	unsigned int i;
  unsigned int j;
  for(j=0;j<6;j++) m[j]=0;
	for(i=0;i<strlen(seq);i++){
		switch(seq[i]){
		case 'A':
			C+= 3;
			H+= 5;   
			N+= 1;   
			O+= 1;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'R':
			C+= 6;
			H+= 12;
			N+= 4;
			O+= 1;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'N':
	    C+= 4;
		  H+= 6;
			N+= 2;
			O+= 2;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;    
		case 'D':
			C+= 4;
			H+= 5;
			N+= 1;
			O+= 3;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'C':
			C+= 3;
			H+= 5;
			N+= 1;
			O+= 1;
			S+= 1;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'Q':
			C+= 5;
			H+= 8;
			N+= 2;
			O+= 2;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'E':
			C+= 5;
			H+= 7;
			N+= 1;
			O+= 3;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'G':
			C+= 2;
			H+= 3;
			N+= 1;
			O+= 1;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'H':
			C+= 6;
			H+= 7;
			N+= 3;
			O+= 1;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'I':
		case 'L':
			C+= 6;
			H+= 11;
			N+= 1;
			O+= 1;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'K':
			C+= 6;
			H+= 12;
			N+= 2;
			O+= 1;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
    case 'M':
			C+= 5;
			H+= 9;
			N+= 1;
			O+= 1;
			S+= 1;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'F':
			C+= 9;
			H+= 9;
	    N+= 1;
		  O+= 1;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'P':
	    C+= 5;
		  H+= 7;
			N+= 1;
			O+= 1;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'S':
			C+= 3;
	    H+= 5;
		  N+= 1;
			O+= 2;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'T':
	    C+= 4;
		  H+= 7;
			N+= 1;
			O+= 2;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'W':
	    C+= 11;
		  H+= 10;
			N+= 2;
			O+= 1;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'Y':
	    C+= 9;
		  H+= 9;
			N+= 1;
			O+= 2;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
		case 'V':
			C+= 5;
			H+= 9;
			N+= 1;
		  O+= 1;
      for(j=0;j<fixMod.size();j++) {
        if(strchr(fixMod[j].aa,seq[i])!=NULL) mass+= fixMod[j].modMass;
      }
			break;
    case '@':
      m[0]++;
      break;
    case '#':
      m[1]++;
      break;
    case '$':
      m[2]++;
      break;
    case '%':
      m[3]++;
      break;
    case '^':
      m[4]++;
      break;
    case '&':
      m[5]++;
      break;
    default:
			break;
		}
	}

	H+=2;
	O++;

	mass+=1.0078246*H;
	mass+=12.0000000*C;
	mass+=14.0030732*N;
	mass+=15.9949141*O;
	mass+=31.972070*S;
  for(j=0;j<varMod.size();j++) mass+=varMod[j].modMass*m[j];

	return mass;
}

//removes redundant peptides
void MascotParser::reduce(){
  unsigned int i;
  unsigned int j;
  
  cout << "Size before: " << m.size() << endl;

  for(i=0;i<m.size()-1;i++){
    for(j=i+1;j<m.size();j++){
      if(!strcmp(m[i].sequence,m[j].sequence)){
        
        //Accept j if it is charge 2 and i is not
        if(m[i].charge!=2 && m[j].charge==2){
          m[i]=m[j];
          m.erase(m.begin()+j);
          j--;
          continue;
        }

        //Accept i if it is charge 2 and j is not
        if(m[i].charge==2 && m[j].charge!=2){
          m.erase(m.begin()+j);
          j--;
          continue;
        }

        //Keep best score
        if(m[i].mascotScore > m[j].mascotScore){
          m.erase(m.begin()+j);
          j--;
        } else {
          m[i]=m[j];
          m.erase(m.begin()+j);
          j--;
        }
      }
    }
  }

  cout << "Size after: " << m.size() << endl;

}

void MascotParser::sortScanNum(){
  qsort(&m[0],m.size(),sizeof(Mascot),compareScanNum);
}


/* For the qsort */
int MascotParser::compareScanNum(const void *p1, const void *p2){
  const Mascot d1 = *(Mascot *)p1;
  const Mascot d2 = *(Mascot *)p2;
  if(d1.scanNum<d2.scanNum) return -1;
  else if(d1.scanNum>d2.scanNum) return 1;
  else return 0;
}
