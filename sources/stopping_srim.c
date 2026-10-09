/* stopping_srim.c -- pouvoirs d'arret electroniques SRIM-2013 / SRIM-2008 / TRIM-98 pour iradina
 * (fichier unique : lecture des tables + calcul + construction des tables par materiau).
 *
 *   stopping_model = 3 : SRIM 2012 (SCOEF03)
 *   stopping_model = 4 : SRIM 2008 (SCOEF03)
 *   stopping_model = 5 : TRIM 1998 (SCOEF.95A, SCOEF.95B)
 *
 * Reconstruit par decompilation de SR.exe (STOP2012 / STOP2006 / TRIM-98).
 *
 * Fichiers de donnees, dans DirectoryData :
 *   models 3 et 4 : SCOEF03.DAT (texte) ou SCOEF03.da_ (format compresse SZDD d'origine).
 *   model 5       : SCOEF.95A, SCOEF.95B, et SCOEFGAS.95 (seulement si un materiau est gazeux)
 * (les noms se changent dans les macros ci-dessous)
 *
 * Options par materiau (SetSRIMMaterialOptions, a appeler avant PrepareStoppingTablesSRIM ;
 * defaut : solide, correction de compose = 1) :
 *   gas      = 1 : cible gazeuse (coefficients "gaz" de Bethe en SRIM, SCOEFGAS.95 en TRIM-98)
 *   compound    : facteur multiplicatif de correction de compose (modeles 3 et 4 ; <=0 : 1)
 * Difference SRIM-2013 / SRIM-2008 (modele 3 / 4) : correction supplementaire de He dans H ou C
 * quand gas=1 et compound != 1 (absente en 2008). Sans gaz ni compose les deux sont identiques.
 * Le pouvoir d'arret d'un materiau = somme de Bragg ponderee par ElementsConc.
 *
 * Unites internes de calcul : energie en keV/u ; S en eV/(1e15 atomes/cm2).
 * Table par (materiau, projectile) : dE/dx [eV/nm] sur une grille log de 1 eV/u a 1e10 eV/u
 * (70 points par decade), interpolation lineaire ; sous 1 eV/u : ~sqrt(E) ; au-dessus : ~1/E.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "target.h"
#include "utils.h"
#include "stopping_srim.h"

#define SRIM_FILE_A  "SCOEF03.DAT"
#define SRIM_FILE_B  "SCOEF03.da_"
#define TRIM_FILE_A  "SCOEF.95A"
#define TRIM_FILE_B  "SCOEF.95B"
#define TRIM_FILE_G  "SCOEFGAS.95"

#define SRIM_NPD  70                       /* points par decade */
#define SRIM_NDEC 10                       /* de 1 eV/u a 1e10 eV/u */
#define SRIM_NG   (SRIM_NPD*SRIM_NDEC+1)
#define SRIM_EMAX_EVU 1.0e10

static float *SrimTab[MAX_NO_MATERIALS][MAX_ELEMENT_NO];

/* donnees chargees */
static double S12[94][100];     /* SCOEF03 : S12[Z][j] = colonne numerique j (1..93) */
static char   S12have[94];
static double S98[94][55];      /* TRIM-98 : S98[Z][col], Z=1..93 (93 = grille), col=1..54 */
static double G98[25][12];      /* SCOEFGAS.95 : G98[ligne][col], ligne 1..24, col 1..11 */

/* options par materiau */
static int    MatGas[MAX_NO_MATERIALS];
static double MatCompound[MAX_NO_MATERIALS];   /* 0 = non defini -> 1 */

void SetSRIMMaterialOptions(int material,int gas,double compound)
{
  if(material<0||material>=MAX_NO_MATERIALS){return;}
  MatGas[material]=(gas!=0);
  MatCompound[material]=compound;
}

/* ================================================================= SRIM 2013 / 2008 ===== */
#define E_MIN 2.0
#define AMU_KEV 931494.0
#define TWO_MEC2_EV 1022000.0
#define INV_ALPHA 137.036
#define K_UNIT 1961.16
#define H_E_MIN 2.0
#define H_E_MAX 30000.0
#define PHASE_RATIO 0.92
#define E_CUTOFF 1.0e7        /* keV/u : au-dela, S = 0 */

/* decompression Microsoft SZDD (.da_) ; buffer malloc (taille dans *n) ou NULL */
static unsigned char *szdd_expand(const unsigned char *raw, size_t rawlen, size_t *n)
{
  static const unsigned char magic[8]={'S','Z','D','D',0x88,0xf0,0x27,0x33};
  unsigned char *out,win[4096];
  size_t size,i=14,o=0;
  int pos=4096-16,ctl,bit,k;
  if(rawlen<14||memcmp(raw,magic,8)){return NULL;}
  size=(size_t)raw[10]|((size_t)raw[11]<<8)|((size_t)raw[12]<<16)|((size_t)raw[13]<<24);
  out=(unsigned char*)malloc(size+1);
  if(out==NULL){return NULL;}
  memset(win,' ',sizeof win);
  while(i<rawlen&&o<size){
    ctl=raw[i++];
    for(bit=0;bit<8;bit++){
      if(o>=size||i>=rawlen){break;}
      if((ctl>>bit)&1){
	unsigned char c=raw[i++];
	out[o++]=c; win[pos]=c; pos=(pos+1)&4095;
      } else {
	int lo,hi,off,len;
	if(i+1>=rawlen){break;}
	lo=raw[i]; hi=raw[i+1]; i+=2;
	off=lo|((hi&0xF0)<<4); len=(hi&0x0F)+3;
	for(k=0;k<len&&o<size;k++){
	  unsigned char c=win[(off+k)&4095];
	  out[o++]=c; win[pos]=c; pos=(pos+1)&4095;
	}
      }
    }
  }
  out[o]=0; *n=o;
  return out;
}

/* texte de SCOEF03 : chaque ligne = "symbole" "nom" valeurs... (separees par espaces ou virgules) */
static int parse_scoef(const char *text,size_t len)
{
  const char *p=text,*end=text+len,*eol;
  int nz=0;
  memset(S12,0,sizeof S12); memset(S12have,0,sizeof S12have);
  while(p<end){
    char *line,*q;
    double vals[128];
    int ntok=0,nv=0,Z,j;
    size_t ll;
    eol=(const char*)memchr(p,'\n',end-p);
    if(eol==NULL){eol=end;}
    ll=eol-p;
    line=(char*)malloc(ll+1);
    if(line==NULL){return -1;}
    memcpy(line,p,ll); line[ll]=0;
    p=(eol<end)?eol+1:end;
    q=line;
    while(*q){
      while(*q&&(isspace((unsigned char)*q)||*q==','||*q=='\r')){q++;}
      if(!*q){break;}
      if(*q=='"'){ /* symbole / nom */
	char *e2=strchr(q+1,'"');
	if(e2==NULL){e2=q+strlen(q);}
	ntok++;
	q=*e2?e2+1:e2;
      } else {
	char *e2=q;
	while(*e2&&!isspace((unsigned char)*e2)&&*e2!=','){e2++;}
	if(ntok<2){ntok++;}
	else if(nv<128){vals[nv++]=strtod(q,NULL);}
	q=e2;
      }
    }
    free(line);
    if(nv==0){continue;}
    Z=(int)vals[0];
    if(Z<1||Z>93){continue;}
    S12have[Z]=1;
    for(j=0;j<nv&&j<93;j++){S12[Z][j+1]=vals[j];}
    nz++;
  }
  return nz;
}

static int load_scoef03(const char *path)
{
  FILE *f=fopen(path,"rb");
  unsigned char *raw;
  long n;
  int rc;
  if(f==NULL){return -1;}
  fseek(f,0,SEEK_END); n=ftell(f); fseek(f,0,SEEK_SET);
  raw=(unsigned char*)malloc(n+1);
  if(raw==NULL){fclose(f);return -1;}
  if(fread(raw,1,n,f)!=(size_t)n){fclose(f);free(raw);return -1;}
  fclose(f);
  if(n>=4&&!memcmp(raw,"SZDD",4)){
    size_t m;
    unsigned char *t=szdd_expand(raw,n,&m);
    if(t==NULL){free(raw);return -1;}
    rc=parse_scoef((const char*)t,m);
    free(t);
  } else {
    rc=parse_scoef((const char*)raw,n);
  }
  free(raw);
  return rc>0?0:-1;
}

static double poly_ln(int z,int first,int last,double x)
{
  double s=0.0,p=1.0;
  int c;
  for(c=first;c<=last;c++){s+=S12[z][c]*p; p*=x;}
  return s;
}

/* coeur Bethe (0x43C280) */
static double bethe_core(int z1,int z2,double e,int gas)
{
  double beta2=1.0-1.0/pow(1.0+e/AMU_KEV,2);
  double ln_e=log(e/1000.0);
  double l0=log(TWO_MEC2_EV*beta2/(1.0-beta2))-beta2;
  double shell=poly_ln(z2,70,80,ln_e);
  double phase=gas?0.0:poly_ln(z2,81,88,ln_e);
  double i_ev=S12[z2][67];
  if(gas&&z1==1){i_ev=S12[z2][68];}
  else if(gas&&z1>1){i_ev=S12[z2][69];}
  double y2=(z1/INV_ALPHA)*(z1/INV_ALPHA)/beta2;
  double bloch=-y2*(1.20206-y2*(1.042-0.8549*y2+0.343*y2*y2));
  double v60=0.001*e;
  double va0=gas?1.35/pow(e,0.4):1.5/pow(e,0.4)+45000.0/(z2*pow(e,1.6));
  double barkas=z1*v60*va0/(v60+va0);
  double big_l=l0-phase-shell-log(i_ev)+barkas+bloch;
  return (double)z1*z1*z2/(K_UNIT*beta2)*big_l;
}

/* ajustement de type "proton" a 10 coefficients a[1..10] (0x43AB40...) */
static double proton_fit(const double *a,double e)
{
  double s_low=a[1]*pow(e,a[2])+a[3]*pow(e,a[4])+a[9]*pow(e,0.5);
  double s_high=(a[5]/pow(e,a[6]))*(log(a[7]/e+a[8]*e)+a[10]/pow(e,0.2));
  return s_low*s_high/(s_low+s_high);
}

/* H, He, Li (z1<=3) */
static double light_ion(int z1,int z2,double e,int gas)
{
  int c0,k;
  double e_max,e_core,a[11],ee,s,s_core;
  switch(z1){
  case 1:  c0=11; e_max=1000.0; break;
  case 2:  c0=21; e_max=2000.0; break;
  default: c0=31; e_max=2000.0; break;
  }
  e_core=1.2*e_max;
  if(e>e_core){return bethe_core(z1,z2,e,gas);}
  for(k=1;k<=10;k++){a[k]=S12[z2][c0+k-1];}
  ee=e<E_MIN?E_MIN:e;
  if(ee>e_max){ee=e_max;}
  s=proton_fit(a,ee);
  if(e>e_max){
    s_core=bethe_core(z1,z2,e_core,gas);
    return s+(s_core-s)*(e-e_max)/(e_core-e_max);
  }
  if(e<=E_MIN){return s*pow(e/E_MIN,0.5);}
  return s;
}

/* ions lourds (0x43D100) : echelle sur He avec charge effective */
static double heavy_ion(int z1,int z2,double e,int gas)
{
  double ve8=e<H_E_MIN?H_E_MIN:e,ln_e,pa,pb,s_he,out;
  double l1,l2,l3,l5,lx,g_e,g_m;
  if(ve8>H_E_MAX){ve8=H_E_MAX;}
  ln_e=log(ve8);
  pa=poly_ln(z1,41,53,ln_e);
  pb=poly_ln(z1,54,66,ln_e);
  s_he=light_ion(2,z2,ve8,gas);
  out=s_he*pa*((double)z1*z1)/4.0*pb;
  if(e>=H_E_MIN&&e<=H_E_MAX){return out;}
  if(e<H_E_MIN){return out*pow(e/H_E_MIN,0.5);}
  l1=S12[z1][89]; l2=S12[z1][90]; l3=S12[z1][91]; l5=S12[z1][93];
  lx=log10(e);
  g_e=l1+(1.0-l1)/(1.0+exp(-(lx-l2)/l3))+l5*exp(-fabs(l2+1.3-lx));
  lx=log10(H_E_MAX);
  g_m=l1+(1.0-l1)/(1.0+exp(-(lx-l2)/l3))+l5*exp(-fabs(l2+1.3-lx));
  return out*(g_e*bethe_core(2,z2,e,gas))/(g_m*bethe_core(2,z2,H_E_MAX,gas));
}

/* aiguillage (0x439530) : facteur de phase x compose, estompe a haute energie, puis S.
   gas : cible gazeuse ; compound : correction de compose ; version : 2012 ou 2006 */
static double srim12_S(int z1,int z2,double e,int gas,double compound,int version)
{
  static const int gas_el[14]={1,2,7,8,9,10,17,18,35,36,53,54,85,86};
  static const int noble[6]={2,10,18,36,54,86};
  double v30,a,v40,s;
  int i,natural=0,is_noble=0;
  if(e>E_CUTOFF){return 0.0;}
  if(compound<=0.0){compound=1.0;}
  for(i=0;i<6;i++){if(noble[i]==z2){is_noble=1;}}
  for(i=0;i<14;i++){if(gas_el[i]==z2){natural=1;}}
  v30=compound;                                      /* gas_requested = gas */
  if(!is_noble&&gas!=natural){v30*=natural?PHASE_RATIO:1.0/PHASE_RATIO;}
  a=1.48*(sqrt(e/25.0)-(z1==1?7.0:4.0));
  if(a>700.0){a=700.0;}
  v40=exp(a)+1.0;
  v30=1.0+(v30-1.0)/v40;
  if(version>=2012&&z1==2&&(z2==1||z2==6)&&compound!=1.0&&gas==1){
    if(e<=10.0){v30*=0.7;}
    else if(e<=60.0){v30*=0.15*log(e-10.0)/log(50.0)+0.7;}
    else if(e<=200.0){v30*=0.85+0.18*(e-60.0)/140.0;}
    else if(e<=400.0){v30*=1.03-0.03*(e-200.0)/200.0;}
  }
  s=(z1<=3)?light_ion(z1,z2,e,gas):heavy_ion(z1,z2,e,gas);
  return v30*s;                                      /* eV/(1e15 at/cm2) */
}

/* ================================================================= TRIM 1998 ============ */
/* A[Z][col] : col 1..8 donnees atomiques, 9..16 proton, 17..20 haute energie,
   21..39 "Qual"+18 lambda d'ecrantage, 40..54 "Qual"+14 corrections de Fermi ; ligne 93 = grille */

static int read_rows(const char *path,int nskip,int maxrows,int maxcols,int col0,double *dst,int stride)
{
  FILE *f=fopen(path,"r");
  char line[4096];
  int skipped=0,row=0,col;
  if(f==NULL){return -1;}
  while(fgets(line,sizeof line,f)){
    const char *p=line;
    char *q=line,*end;
    if(skipped<nskip){skipped++;continue;}
    while(*p&&isspace((unsigned char)*p)){p++;}
    if(!*p){continue;}
    if(row>=maxrows){break;}
    for(col=0;col<maxcols;col++){
      double v=strtod(q,&end);
      if(end==q){break;}
      dst[(row+1)*stride+col0+col]=v;
      q=end;
    }
    row++;
  }
  fclose(f);
  return row;
}

static int load_trim98(int need_gas)
{
  char p[1200];
  memset(S98,0,sizeof S98); memset(G98,0,sizeof G98);
  snprintf(p,sizeof p,"%s/%s",DirectoryData,TRIM_FILE_A);
  if(read_rows(p,2,93,16,1,&S98[0][0],55)<0){return -1;}
  snprintf(p,sizeof p,"%s/%s",DirectoryData,TRIM_FILE_B);
  if(read_rows(p,2,93,38,17,&S98[0][0],55)<0){return -1;}
  snprintf(p,sizeof p,"%s/%s",DirectoryData,TRIM_FILE_G);
  if(read_rows(p,2,24,11,1,&G98[0][0],12)<0&&need_gas){return -1;}
  return 0;
}

static double pw(double x,double y){return exp(y*log(x));}

/* S "proton" ; E en keV/u ; gidx = 0 (solide) ou 1..8 (gaz : H,He,N,O,Ne,Ar,Kr,Xe) */
static double proton98(int Z1,int Z2,double E,int gidx)
{
  const double *A=S98[Z2];
  double x,SL,SH,S,L;
  if(E<=10000.0){
    x=E>10.0?E:10.0;
    if(gidx){
      const double *G=G98[Z1>2?gidx+16:(Z1>1?gidx+8:gidx)];
      SL=G[6]*pw(x,G[7])+G[4]*pw(x,G[5]);
      SH=(G[8]/pw(x,G[9]))*log(G[10]/x+G[11]*x);
    } else {
      SL=A[11]*pw(x,A[12])+A[9]*pw(x,A[10]);
      SH=(A[13]/pw(x,A[14]))*log(A[15]/x+A[16]*x);
    }
    S=SH*SL/(SH+SL);
    if(E<=10.0){S*=pw(E/10.0,Z2<7?0.35:0.45);}
    return S;
  }
  L=log(E)/E;
  return A[17]+A[18]*L+A[19]*L*L+A[20]/L;
}

/* interpolation dans la ligne Z, abscisses = ligne 93 */
static double interp98(int Z,int lo,int hi,double x,int clamp_hi)
{
  const double *A93=S98[93],*r=S98[Z];
  int i=lo;
  while(i<=hi&&A93[i]<x){i++;}
  i--;
  if(i>clamp_hi){i=clamp_hi;}
  if(i<lo){i=lo;}
  return r[i]+(r[i+1]-r[i])*(x-A93[i])/(A93[i+1]-A93[i]);
}

static double fermi_corr98(int Z2,double E)
{
  const double *A93=S98[93],*r=S98[Z2];
  double x=E<9999.0?E:9999.0;
  int i=41;
  while(i<=53&&A93[i]<=x){i++;}
  i--;
  if(i<41){i=41;}
  if(i>53){i=53;}
  return r[i]+(r[i+1]-r[i])*(x-A93[i])/(A93[i+1]-A93[i]);
}

static double stop_He98(int Z2,double E,int gidx)
{
  double x=E>1.0?E:1.0,L=log(x),a,h,f,S;
  a=0.2865+0.1266*L-0.001429*L*L+0.02402*pow(L,3)-0.01135*pow(L,4)+0.001475*pow(L,5);
  if(a>30.0){a=30.0;}
  h=1.0-exp(-a);
  f=1.0+(0.007+0.00005*Z2)*exp(-(7.6-L)*(7.6-L));
  S=h*f*f*proton98(2,Z2,x,gidx)*4.0;
  if(E<=1.0){S*=sqrt(E);}
  return S;
}

static double stop_heavy98(int Z1,int Z2,double E,int gidx)
{
  double vF=S98[Z2][7];
  double z13=pow((double)Z1,1.0/3.0),z23=pow((double)Z1,0.6667);
  double vr=sqrt(E/25.0)/vF,vrel,yr,t,q,lam,u,zeta,lnE,gamma,yrmin,zg,Sp;
  double vmin,tt,v1,E1,ex;
  if(vr<1.0){vrel=(vr*vr*2.0/3.0+1.0-pow(vr,4)/15.0)*vF*3.0/4.0;}
  else{vrel=vF*vr*(1.0+1.0/(5.0*vr*vr));}
  yr=vrel/z23;
  if(yr<0.13){yr=0.13;}
  if(yr<1.0/z23){yr=1.0/z23;}
  t=-0.803*pow(yr,0.3)+1.3167*pow(yr,0.6)+0.38157*yr+0.008983*yr*yr;
  if(t>50.0){t=50.0;}
  q=1.0-exp(-t);
  if(q<0.0){q=0.0;}
  if(q>1.0){q=1.0;}
  lam=interp98(Z1,22,39,q,38)/z13;
  u=lam*4.0*vF/1.919;
  zeta=(1.0-q)*(1.0/(2.0*vF*vF))*log(u*u+1.0)+q;
  lnE=log(E);
  if(lnE<0.0){lnE=0.0;}
  gamma=(((0.0015*Z2+0.08)*(1.0/((double)Z1*Z1)))*exp(-(7.6-lnE)*(7.6-lnE))+1.0)*zeta;
  yrmin=1.0/z23>0.13?1.0/z23:0.13;
  if(yrmin<yr){
    zg=Z1*gamma;
    Sp=proton98(Z1,Z2,E,gidx);
    return zg*zg*Sp*fermi_corr98(Z2,E);
  }
  /* domaine basse vitesse */
  vmin=z23*0.13>1.0?z23*0.13:1.0;
  tt=vmin*vmin-0.8*vF*vF;
  if(tt<0.0){tt=0.0;}
  v1=(sqrt(tt)+vmin)*0.5;
  E1=v1*v1*25.0;
  Sp=proton98(Z1,Z2,E1,gidx)*fermi_corr98(Z2,E1);
  if(Z1==3){ex=0.55;}
  else if(Z2<7){ex=0.375;}
  else if(Z1>=18){ex=0.47;}
  else if(Z2==14||Z2==32){ex=0.375;}
  else{ex=0.47;}
  zg=Z1*gamma;
  return pow(E/E1,ex)*zg*zg*Sp;
}

/* gas = 0/1 ; gidx = rang de Z2 dans la liste des 8 gaz de SCOEFGAS.95, ou 0 */
static double trim98_S(int Z1,int Z2,double E,int gas)
{
  static const int gas_z[8]={1,2,7,8,10,18,36,54};
  int gidx=0,i;
  if(gas){for(i=0;i<8;i++){if(gas_z[i]==Z2){gidx=i+1;break;}}}
  if(Z1==1){return proton98(1,Z2,E,gidx);}
  if(Z1==2){return stop_He98(Z2,E,gidx);}
  return stop_heavy98(Z1,Z2,E,gidx);
}

static int trim98_is_gas(int z)
{
  return z==1||z==2||z==7||z==8||z==10||z==18||z==36||z==54;
}

/* ================================================================= interface iradina ==== */

static int srim_load_data(int model)
{
  char path[1200];
  if(model==5){
    int need_gas=0,m;
    for(m=0;m<NumberOfMaterials;m++){if(MatGas[m]){need_gas=1;}}
    if(load_trim98(need_gas)!=0){
      message_error(-2,"Cannot read TRIM-98 files (%s, %s%s) in %s\n",TRIM_FILE_A,TRIM_FILE_B,need_gas?", "TRIM_FILE_G:"",DirectoryData);
      return -2;
    }
    return 0;
  }
  snprintf(path,sizeof path,"%s/%s",DirectoryData,SRIM_FILE_A);
  if(load_scoef03(path)!=0){
    snprintf(path,sizeof path,"%s/%s",DirectoryData,SRIM_FILE_B);
    if(load_scoef03(path)!=0){
      message_error(-2,"Cannot read SRIM stopping data file %s or %s in %s\n",
		    SRIM_FILE_A,SRIM_FILE_B,DirectoryData);
      return -2;
    }
  }
  return 0;
}

/* S en eV/(1e15 atomes/cm2) pour l'ion Z1 sur l'element Z2, E en keV/u */
static double srim_element_S(int model,int z1,int z2,double e_kev_u,int gas,double compound)
{
  double s;
  if(z1<1||z1>ZMAXdpass||z2<1||z2>ZMAXdpass){return 0.0;}
  if(model==5){
    s=trim98_S(z1,z2,e_kev_u,gas);
  } else {
    if(!S12have[z1]||!S12have[z2]){return 0.0;}
    s=srim12_S(z1,z2,e_kev_u,gas,compound,(model==3)?2012:2006);
  }
  if(!isfinite(s)||s<0.0){s=0.0;}
  return s;
}

int PrepareStoppingTablesSRIM(int model)
{
  int i,j,k,ie,result,g;
  double sum,e_kev_u;

  if(model<3||model>5){return -2;}
  result=srim_load_data(model);
  if(result!=0){return result;}

  for(i=0;i<MAX_NO_MATERIALS;i++){for(j=0;j<MAX_ELEMENT_NO;j++){SrimTab[i][j]=NULL;}}

  for(i=0;i<NumberOfMaterials;i++){
    for(j=1;j<=ZMAXdpass;j++){ /* tous les projectiles possibles */
      if(existing_elements[j]!=1){continue;}
      SrimTab[i][j]=(float*)calloc(SRIM_NG,sizeof(float));
      if(SrimTab[i][j]==NULL){return -1;}
      if(ListOfMaterials[i].Is_Vacuum==1||ListOfMaterials[i].Density<=0.0f){continue;}
      for(ie=0;ie<SRIM_NG;ie++){
	e_kev_u=1.0e-3*pow(10.0,(double)ie/SRIM_NPD);
	sum=0.0;
	g=MatGas[i];
	for(k=0;k<ListOfMaterials[i].ElementCount;k++){ /* somme de Bragg */
	  if(model==5&&g&&!trim98_is_gas(ListOfMaterials[i].ElementsZ[k])){g=0;} /* TRIM-98 : un element non gazeux annule gaz (comme SR.exe) */
	  sum+=ListOfMaterials[i].ElementsConc[k]*srim_element_S(model,j,ListOfMaterials[i].ElementsZ[k],e_kev_u,g,MatCompound[i]);
	}
	/* eV/(1e15 at/cm2) -> eV/nm : *1e-15 * densite (at/cm3) * 1e-7 */
	SrimTab[i][j][ie]=(float)(sum*1e-15*1e-7*ListOfMaterials[i].Density);
      }
    }
  }
  return 0;
}

float ElectronicStoppingSRIM(int ionZ,float ionM,float E,int material)
{
  float *y;
  double Epu,x,t;
  int i;
  y=SrimTab[material][ionZ];
  if(y==NULL){return 0.0f;}
  Epu=(double)E/ionM;                /* eV/u */
  if(Epu<=0.0){return 0.0f;}
  x=log10(Epu)*SRIM_NPD;
  if(x<=0.0){return (float)(y[0]*sqrt(Epu));}                           /* E0 = 1 eV/u */
  if(x>=(double)(SRIM_NG-1)){return (float)(y[SRIM_NG-1]*(SRIM_EMAX_EVU/Epu));}
  i=(int)x;
  t=x-i;
  return (float)(y[i]+t*(y[i+1]-y[i]));
}
