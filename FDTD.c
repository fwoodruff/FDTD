//
//  main.c
//  scruffy4mm
//
//  Created by Frederick Benjamin Woodruff on 19/03/2017.
//  Copyright © 2017 Frederick Benjamin Woodruff. All rights reserved.
//
//
//  main.c
//  5mmWaveplate
//
//  Created by Frederick Benjamin Woodruff on 28/02/2017.
//  Copyright © 2017 Frederick Benjamin Woodruff. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <pthread.h>

#define DIMS 3
#define NTHREADS 5
#define strM 1023
#define speedC 299792458.0 // m/s
#define VAC_PMBY 4e-7*M_PI
#define VAC_PTVY 1./(VAC_PMBY*speedC*speedC)
#define FLOAT float
#define DOUBLE double // for extension to complex

#define max(a,b) \
({ __typeof__ (a) _a = (a); \
__typeof__ (b) _b = (b); \
_a > _b ? _a : _b; })
#define min(a,b) \
({ __typeof__ (a) _a = (a); \
__typeof__ (b) _b = (b); \
_a < _b ? _a : _b; })


char ExReferenceFilename[strM] = "a4mmDATAscrappy/ExPaperRef.csv"; // free space
char EyReferenceFilename[strM] = "a4mmDATAscrappy/EyPaperRef.csv";
char EzReferenceFilename[strM] = "a4mmDATAscrappy/EzPaperRef.csv";

char ExSampleFilename[strM] = "a4mmDATAscrappy/ExPaperSample.csv"; // pillar array
char EySampleFilename[strM] = "a4mmDATAscrappy/EyPaperSample.csv";
char EzSampleFilename[strM] = "a4mmDATAscrappy/EzPaperSample.csv";

char materialFile[strM] = "a4mmDATAscrappy/PaperMaterial.csv";

//char ExPtReferenceFilename[strM] = "ExPtFieldFSPaperf.csv";
//char EyPtReferenceFilename[strM] = "EyPtFieldFSPaperf.csv";
//char EzPtReferenceFilename[strM] = "EzPtFieldFSPaperf.csv";

//char ExPtSampleFilename[strM] = "ExPtFieldPAPaperf.csv";
//char EyPtSampleFilename[strM] = "EyPtFieldPAPaperf.csv";
//char EzPtSampleFilename[strM] = "EzPtFieldPAPaperf.csv";

char dExReferenceFileName[strM]= "a4mmDATAscrappy/dExPaper4mmReference.csv";
char dEyReferenceFileName[strM]= "a4mmDATAscrappy/dEyPaper4mmReference.csv";
char dEzReferenceFileName[strM]= "a4mmDATAscrappy/dEzPaper4mmReference.csv";

char dHxReferenceFileName[strM]= "a4mmDATAscrappy/dHxPaper4mmReference.csv";
char dHyReferenceFileName[strM]= "a4mmDATAscrappy/dHyPaper4mmReference.csv";
char dHzReferenceFileName[strM]= "a4mmDATAscrappy/dHzPaper4mmReference.csv";

char dExSampleFileName[strM]= "a4mmDATAscrappy/dExPaper4mmSample.csv";
char dEySampleFileName[strM]= "a4mmDATAscrappy/dEyPaper4mmSample.csv";
char dEzSampleFileName[strM]= "a4mmDATAscrappy/dEzPaper4mmSample.csv";

char dHxSampleFileName[strM]= "a4mmDATAscrappy/dHxPaper4mmSample.csv";
char dHySampleFileName[strM]= "a4mmDATAscrappy/dHyPaper4mmSample.csv";
char dHzSampleFileName[strM]= "a4mmDATAscrappy/dHzPaper4mmSample.csv";


// Domain
const DOUBLE CouS = 1.8; // 1/Courant stability factor
const FLOAT deltaS = 60e-6/2;
const DOUBLE domainX = 0.006; // metres
const DOUBLE domainY = 0.006; // metres
const DOUBLE domainZ = 0.006; // metres Should be 0.010, smaller value is an approximation
const DOUBLE timePeriod = 5e-11; // seconds



// Current Source
DOUBLE t0=1.2e-12; // seconds
DOUBLE FWHM = 1e-12; //  seconds (pulse's full width half maximum)
DOUBLE frq = 0.8e12; // Hz
DOUBLE phase = 0.4; // radians
const FLOAT sourceZ = domainZ; // m (z)
const FLOAT sourceY = domainY; // m (y)
const FLOAT sourceY0 = domainY/2; // m (y)
const FLOAT sourceX0 = 0;


// Waveplate
const DOUBLE nPaper = 1.41; // refractive index
const DOUBLE alphaPaper = 1460;
const DOUBLE nAir = 1.00027;
FLOAT thPaper = 120e-6; // thickness of a sheet

//float nVac = 2.4;
FLOAT thVac = 180e-6;

FLOAT width = 0.004; // x-direction width in metres
FLOAT xStart = 0.001+domainX-domainY; // metres



FILE* ExReferenceFile;
FILE* EyReferenceFile;
FILE* EzReferenceFile;

FILE* ExSampleFile;
FILE* EySampleFile;
FILE* EzSampleFile;

FILE* materialFieldOut;

//FILE* ExPtReferenceFile;
//FILE* EyPtReferenceFile;
//FILE* EzPtReferenceFile;

//FILE* ExPtSampleFile;
//FILE* EyPtSampleFile;
//FILE* EzPtSampleFile;

FILE* dExReferenceFile;
FILE* dEyReferenceFile;
FILE* dEzReferenceFile;

FILE* dExSampleFile;
FILE* dEySampleFile;
FILE* dEzSampleFile;

FILE* dHxReferenceFile;
FILE* dHyReferenceFile;
FILE* dHzReferenceFile;

FILE* dHxSampleFile;
FILE* dHySampleFile;
FILE* dHzSampleFile;

const FLOAT deltaT = deltaS/(speedC*CouS); // metres
const int GRIDX = domainX/deltaS;
const int GRIDY = domainY/deltaS;
const int GRIDZ = domainZ/deltaS;
const int GRIDYZ = GRIDY*GRIDZ;


const int MEM = GRIDX * GRIDY * GRIDZ;
int bounds[NTHREADS+1];



const int TSTEPS = timePeriod/deltaT;
DOUBLE sourceJz[TSTEPS];

const DOUBLE sigmaPaper = (alphaPaper*nPaper)/(speedC*VAC_PMBY);
const DOUBLE permPaper = nPaper*nPaper;
const DOUBLE permAir =  nAir*nAir;

// dielectric, conductivity, mu, mcond.
DOUBLE vac[4] = {permAir,0,1,0};
DOUBLE paper[4] = {permPaper,sigmaPaper,1,0};
FLOAT C_vac[4];
FLOAT C_paper[4];


FLOAT getC(char material,char type) {
    if(material=='p') {
        return C_paper[type];
    }
    if (material=='v') {
        return C_vac[type];
    }
    printf("error");
    return 0;
}

FLOAT Ex[MEM];
FLOAT Ey[MEM];
FLOAT Ez[MEM];
FLOAT Hx[MEM];
FLOAT Hy[MEM];
FLOAT Hz[MEM];
char materialType[MEM];


static FLOAT EyAtX0[GRIDY*GRIDZ];
static FLOAT HyAtXmax[GRIDY*GRIDZ];
static FLOAT EzAtX0[GRIDY*GRIDZ];
static FLOAT HzAtXmax[GRIDY*GRIDZ];
static FLOAT ExAtY0[GRIDX*GRIDZ];
static FLOAT HxAtYmax[GRIDX*GRIDZ];
static FLOAT EzAtY0[GRIDX*GRIDZ];
static FLOAT HzAtYmax[GRIDX*GRIDZ];
static FLOAT ExAtZ0[GRIDX*GRIDY];
static FLOAT HxAtZmax[GRIDX*GRIDY];
static FLOAT EyAtZ0[GRIDX*GRIDY];
static FLOAT HyAtZmax[GRIDX*GRIDY];

int IJKf(int,int,int);
DOUBLE current_timestamp();
void writeToFile(FILE*, FLOAT*);
void writePointToFile(FILE*, FLOAT*);
void getBounds();

int IJKf(int i, int j, int k) {
    return i*GRIDY*GRIDZ + j*GRIDZ + k;
}

double current_timestamp() {
    struct timeval te;
    gettimeofday(&te, NULL); // get current time
    double milliseconds = te.tv_sec*1000LL + te.tv_usec/1000; // caculate milliseconds
    // printf("milliseconds: %lld\n", milliseconds);
    return milliseconds;
}

void writeToFile(FILE* file, FLOAT *field) {
    for (int i=GRIDX-GRIDY+1; i<GRIDX-1; i++) {
        for (int j=1; j<GRIDY-1; j++) {
            int IJK = IJKf(i,j,GRIDZ/2);
            fprintf(file, "%.8f; ", field[IJK]);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}

void writeMaterialFile(FILE* file, char *field) {
    for (int i=GRIDX-GRIDY+1; i<GRIDX-1; i++) {
        for (int j=1; j<GRIDY-1; j++) {
            int IJK = IJKf(i,j,GRIDZ/2);
            fprintf(file, "%c;", field[IJK]);
        }
        fprintf(file, "\n");
    }
    fprintf(file, "\n");
}

void writePointToFile(FILE* file, FLOAT *field) {
    int IJK = IJKf(GRIDX-2,GRIDY/2,GRIDZ/2);
    fprintf(file, "%.8f; ", field[IJK]);
}

void takeWallToSingleValueAndWriteToFile(FILE* file, FLOAT *field) {
    double total=0;
    for (int j=0; j<GRIDY; j++) {
        for (int k=0; k<GRIDZ; k++) {
            int IJK = IJKf(GRIDX-2,j,k);
            total+=field[IJK];
        }
    }
    fprintf(file, "%.8f; ", total);
}



// Ca, Cb, Da, Db
void getConsts(FLOAT *Cvalues,DOUBLE *material) {
    Cvalues[0] =(1-material[1]*deltaT/(2*material[0]*VAC_PTVY))/(1+material[1]*deltaT/(2*material[0]*VAC_PTVY));
    Cvalues[1] = (deltaT/(material[0]*VAC_PTVY*deltaS))/(1+material[1]*deltaT/(2*material[0]*VAC_PTVY));
    Cvalues[2] = (1-(material[3]*deltaT/(2*material[2]*VAC_PMBY)))/(1+(material[3]*deltaT/(2*material[2]*VAC_PMBY)));
    Cvalues[3] = (deltaT/(material[2]*VAC_PMBY*deltaS))/(1+(material[3]*deltaT/(2*material[2]*VAC_PMBY)));
}



// calculated parameters

// ____________________________________________________
// MEMORY INIT.

int initABC() {
    for (int m = 0; m<GRIDX*GRIDZ; m++) {
        ExAtY0[m]=0;
        HxAtYmax[m]=0;
        EzAtY0[m]=0;
        HzAtYmax[m]=0;
    }
    for (int m = 0; m<GRIDY*GRIDZ; m++) {
        EyAtX0[m]=0;
        HyAtXmax[m]=0;
        EzAtX0[m]=0;
        HzAtXmax[m]=0;
    }
    for (int m = 0; m<GRIDX*GRIDY; m++) {
        ExAtZ0[m]=0;
        HxAtZmax[m]=0;
        EyAtZ0[m]=0;
        HyAtZmax[m]=0;
    }
    return 1;
}
void getBounds() {
    for (int tid=0; tid<NTHREADS+1; tid++) {
        bounds[tid]=tid*(GRIDX-2)/(NTHREADS)+1;
    }
}
int ABC(int X0Interface, int XmaxInterface, int Y0Interface, int YmaxInterface, int ZmaxInterface, int Z0Interface) {
    
    FLOAT abccoef = ((CouS )-1)/((CouS)+1);
    
    if(X0Interface) {
        //back
        for(int j=0;j<(GRIDY-1);j++) {
            for(int k=0;k<GRIDZ;k++) {
                Ey[0*(GRIDY*GRIDZ)+j*GRIDZ+k] = EyAtX0[j*GRIDZ+k]-abccoef*(Ey[1*GRIDY*GRIDZ+j*GRIDZ+k]-Ey[0*(GRIDY*GRIDZ)+j*GRIDZ+k]);
                EyAtX0[j*GRIDZ+k] = Ey[1*(GRIDY*GRIDZ)+j*GRIDZ+k];
            }
            for(int j=0;j<GRIDY;j++) {
                for(int k=0;k<(GRIDZ-1);k++) {
                    Ez[0*(GRIDY*GRIDZ)+j*GRIDZ+k] = EzAtX0[j*GRIDZ+k]-abccoef*(Ez[1*(GRIDY*GRIDZ)+j*GRIDZ+k]-Ez[0*(GRIDY*GRIDZ)+j*GRIDZ+k]);
                    EzAtX0[j*GRIDZ+k] = Ez[1*(GRIDY*GRIDZ)+j*GRIDZ+k];
                }
            }
        }
    }
    
    if(XmaxInterface) {
        //front
        for(int j=0;j<(GRIDY-1);j++) {
            for(int k=0;k<GRIDZ;k++) {
                Hy[(GRIDX-1)*(GRIDY*GRIDZ)+j*GRIDZ+k] = HyAtXmax[j*GRIDZ+k]-abccoef*(Hy[(GRIDX-2)*(GRIDY*GRIDZ)+j*GRIDZ+k]-Hy[(GRIDX-1)*(GRIDY*GRIDZ)+j*GRIDZ+k]);
                HyAtXmax[j*GRIDZ+k] = Hy[(GRIDX-2)*(GRIDY*GRIDZ)+j*GRIDZ+k];
            }
        }
        for(int j=0;j<GRIDY;j++) {
            for(int k=0;k<(GRIDZ-1);k++) {
                Hz[(GRIDX-1)*GRIDYZ+j*GRIDZ+k] = HzAtXmax[j*GRIDZ+k]-abccoef*(Hz[(GRIDX-2)*GRIDYZ+j*GRIDZ+k]-Hz[(GRIDX-1)*(GRIDY*GRIDZ)+j*GRIDZ+k]);
                HzAtXmax[j*GRIDZ+k] = Hz[(GRIDX-2)*(GRIDY*GRIDZ)+j*GRIDZ+k];
            }
        }
    }
    
    if(Y0Interface) {
        //left
        for(int i=0;i<(GRIDX-1);i++) {
            for(int k=0;k<GRIDZ;k++) {
                Ex[i*GRIDYZ+0*GRIDZ+k] = ExAtY0[i*GRIDZ+k]-abccoef*(Ex[i*GRIDYZ+1*GRIDZ+k]-Ex[i*GRIDYZ+0*GRIDZ+k]);
                ExAtY0[i*GRIDZ+k] = Ex[i*GRIDYZ+1*GRIDZ+k];
            }
        }
        for(int i=0;i<GRIDX;i++) {
            for(int k=0;k<(GRIDZ-1);k++) {
                Ez[i*GRIDYZ+0*GRIDZ+k] = EzAtY0[i*GRIDZ+k]-abccoef*(Ez[i*GRIDYZ+1*GRIDZ+k]-Ez[i*GRIDYZ+0*GRIDZ+k]);
                EzAtY0[i*GRIDZ+k] = Ez[i*GRIDYZ+1*GRIDZ+k];
            }
        }
    }
    
    if(YmaxInterface) {
        //right
        for(int i=0;i<(GRIDX-1);i++) {
            for(int k=0;k<GRIDZ;k++) {
                Hx[i*GRIDYZ+(GRIDY-1)*GRIDZ+k] = HxAtYmax[i*GRIDZ+k]-abccoef*(Hx[i*GRIDYZ+(GRIDY-2)*GRIDZ+k]-Hx[i*GRIDYZ+(GRIDY-1)*GRIDZ+k]);
                HxAtYmax[i*GRIDZ+k] = Hx[i*GRIDYZ+(GRIDY-2)*GRIDZ+k];
            }
        }
        for(int i=0;i<GRIDX;i++) {
            for(int k=0;k<(GRIDZ-1);k++) {
                Hz[i*GRIDYZ+(GRIDY-1)*GRIDZ+k] = HzAtYmax[i*GRIDZ+k]-abccoef*(Hz[i*GRIDYZ+(GRIDY-2)*GRIDZ+k]-Hz[i*GRIDYZ+(GRIDY-1)*GRIDZ+k]);
                HzAtYmax[i*GRIDZ+k] = Hz[i*GRIDYZ+(GRIDY-2)*GRIDZ+k];
            }
        }
    }
    
    if(Z0Interface) {
        //bottom
        
        for(int i=0;i<(GRIDX-1);i++) {
            for(int j=0;j<GRIDY;j++) {
                Ex[i*GRIDYZ+j*GRIDZ+0] = ExAtZ0[i*GRIDY+j]-abccoef*(Ex[i*GRIDYZ+j*GRIDZ+1]-Ex[i*GRIDYZ+j*GRIDZ+0]);
                ExAtZ0[i*GRIDY+j] = Ex[i*GRIDYZ+j*GRIDZ+1];
            }
        }
        
        for(int i=0;i<GRIDX;i++) {
            for(int j=0;j<(GRIDY-1);j++) {
                Ey[i*GRIDYZ+j*GRIDZ+0] = EyAtZ0[i*GRIDY+j]-abccoef*(Ey[i*GRIDYZ+j*GRIDZ+1]-Ey[i*GRIDYZ+j*GRIDZ+0]);
                EyAtZ0[i*GRIDY+j] = Ey[i*GRIDYZ+j*GRIDZ+1];
            }
        }
        
        
    }
    
    if(ZmaxInterface) {
        //top
        
        for(int i=0;i<(GRIDX-1);i++) {
            for(int j=0;j<GRIDY;j++) {
                Hx[i*GRIDYZ+j*GRIDZ+GRIDZ-1] = HxAtZmax[i*GRIDY+j]-abccoef*(Hx[i*GRIDYZ+j*GRIDZ+GRIDZ-2]-Hx[i*GRIDYZ+j*GRIDZ+(GRIDZ-1)]);
                HxAtZmax[i*GRIDY+j] = Hx[i*GRIDYZ+j*GRIDZ+(GRIDZ-2)];
            }
        }
        
        for(int i=0;i<GRIDX;i++) {
            for(int j=0;j<(GRIDY-1);j++) {
                Hy[i*GRIDYZ+j*GRIDZ+GRIDZ-1] = HyAtZmax[i*GRIDY+j]-abccoef*(Hy[i*GRIDYZ+j*GRIDZ+GRIDZ-2]-Hy[i*GRIDYZ+j*GRIDZ+(GRIDZ-1)]);
                HyAtZmax[i*GRIDY+j] = Hy[i*GRIDYZ+j*GRIDZ+(GRIDZ-2)];
            }
        }
    }
    
    return 1;
}


void* Ekernel(void *x) {
    int tid;
    tid = *((int *) x);
    
    int start = bounds[tid];
    int end = bounds[tid+1];
    for (int i=start; i<end; i++) {
        for (int j=1; j<GRIDY-1; j++) {
            int IJ = i*GRIDYZ + j*GRIDZ;
            for (int k=1; k<GRIDZ-1; k++) {
                int IJK = IJ + k;
                int Ijk = IJK + GRIDYZ;;
                int iJk = IJK + GRIDZ;
                int ijK = IJK + 1;
                Ez[IJK] = getC(materialType[IJK], 0)*Ez[IJK] + getC(materialType[IJK], 1)*(Hy[Ijk]-Hy[IJK]-Hx[iJk]+Hx[IJK]);
                Ex[IJK] = getC(materialType[IJK], 0)*Ex[IJK] + getC(materialType[IJK], 1)*(Hz[iJk]-Hz[IJK]-Hy[ijK]+Hy[IJK]);
                Ey[IJK] = getC(materialType[IJK], 0)*Ey[IJK] + getC(materialType[IJK], 1)*(Hx[ijK]-Hx[IJK]-Hz[Ijk]+Hz[IJK]);
            }
        }
    }
    return NULL;
}
void* Hkernel(void *x) {
    int tid;
    tid = *((int *) x);
    int start = bounds[tid];
    int end = bounds[tid+1];
    for (int i=start; i<end; i++) {
        for (int j=1; j<GRIDY-1; j++) {
            int IJ = i*GRIDYZ + j*GRIDZ;
            for (int k=1; k<GRIDZ-1; k++) {
                int IJK = IJ + k;
                int iJK = IJK - GRIDYZ;
                int IjK = IJK - GRIDZ;
                int IJk = IJK - 1;
                Hx[IJK] = getC(materialType[IJK], 2)*Hx[IJK] + getC(materialType[IJK], 3)*(Ey[IJK]-Ey[IJk]-Ez[IJK]+Ez[IjK]);
                Hy[IJK] = getC(materialType[IJK], 2)*Hy[IJK] + getC(materialType[IJK], 3)*(Ez[IJK]-Ez[iJK]-Ex[IJK]+Ex[IJk]);
                Hz[IJK] = getC(materialType[IJK], 2)*Hz[IJK] + getC(materialType[IJK], 3)*(Ex[IJK]-Ex[IjK]-Ey[IJK]+Ey[iJK]);
            }
        }
    }
    return NULL;
}

void E_CPU_para() {
    // Calculates E and H fields at the next timestep
    pthread_t threads[NTHREADS];
    int rc,ti;
    int thread_args[NTHREADS];
    
    
    for (ti=0; ti<NTHREADS; ++ti)
    {
        thread_args[ti] = ti;
        rc = pthread_create(&threads[ti], NULL, &Ekernel,(void *)&thread_args[ti]);
    }
    for (ti=0; ti<NTHREADS; ++ti) {
        rc = pthread_join(threads[ti], NULL);
    }
}
void H_CPU_para() {
    // Calculates H fields at the next timestep
    pthread_t threads[NTHREADS];
    int rc,ti;
    int thread_args[NTHREADS];
    
    
    for (ti=0; ti<NTHREADS; ++ti)
    {
        thread_args[ti] = ti;
        rc = pthread_create(&threads[ti], NULL, &Hkernel,(void *)&thread_args[ti]);
    }
    for (ti=0; ti<NTHREADS; ++ti) {
        rc = pthread_join(threads[ti], NULL);
    }
}


void initArraysToZero() {
    for (int ijk=0; ijk<MEM; ijk++) {
        
        materialType[ijk]='v';
        Ex[ijk]=0;
        Ey[ijk]=0;
        Ez[ijk]=0;
        Hx[ijk]=0;
        Hy[ijk]=0;
        Hz[ijk]=0;
    }
}


void makesource() {
    DOUBLE twoCSq = 2*(FWHM/2.35482)*(FWHM/2.35482);
    printf("\nsource waveform:\n");
    for (int n=0; n<TSTEPS; n++) {
        DOUBLE t = n*deltaT;
        DOUBLE tx = (t-t0);
        sourceJz[n] = pow(2,-tx*tx/(twoCSq)) * sin(t*frq*2*M_PI + phase);
        if (t<4e-12 && n%3==0) {
            printf("%f; ",sourceJz[n]);
        }
    }
    printf("\n");
}
void source(int n) {
    DOUBLE source = sourceJz[n];
    
    int starty = max(((sourceY0/domainY)*GRIDX*(1-sourceY/domainY)),0);
    int endy = min(((sourceY0/domainY)*GRIDX*(1+sourceY/domainY)),GRIDY);
    int startz = max((0.5*GRIDZ*(1-sourceZ/domainZ)),0);
    int endz = min((0.5*GRIDZ*(1+sourceZ/domainZ)),GRIDZ);
    
    DOUBLE theta = 45; // degrees in yz relative to +z
    
    DOUBLE sintheta = sin(theta * M_PI/180);
    DOUBLE costheta = cos(theta);
    
    for (int i=(sourceX0/domainX)*GRIDX +1; i<(sourceX0/domainX)*GRIDX +3; i++) {
        for (int j=starty; j<endy; j++) {
            for (int k=startz; k<endz; k++) {
                int IJK = IJKf(i,j,k);
                int rsq = ((GRIDY/2)-j)*((GRIDY/2)-j) + ((GRIDZ/2)-k)*((GRIDZ/2)-k);
                int norm = (GRIDY/2)*(GRIDY/2) + (GRIDZ/2)*(GRIDZ/2);
                double mult=pow(2,-6*rsq/norm);
                Ez[IJK]+=source*costheta*mult;
                Ey[IJK]+=source*sintheta*mult;
            }
        }
    }
}


void makeplanar(double theta_l, double phi_l, double x0_l,double y0_l,double z0_l,double d_l, double thickness) {
    double cth=cos(theta_l*M_PI/180);
    double sth=sin(theta_l*M_PI/180);
    double cph=cos(phi_l*M_PI/180);
    double sph=sin(phi_l*M_PI/180);
    
    double yAng = sth*cph;
    double xAng=sth*sph;
    double zAng=cth;
    
    for (int xi=0; xi<GRIDX; xi++) {
        for (int yi=0; yi<GRIDY; yi++) {
            for (int zi=0; zi<GRIDZ; zi++) {
                double xq= xi*domainX/GRIDX;
                double yq= yi*domainY/GRIDY;
                double zq= zi*domainZ/GRIDZ;
                
                int conditionP = (xq*xAng + yq*yAng + zq*zAng<(d_l/2.)+x0_l*xAng + y0_l*yAng + z0_l*zAng);
                int conditionM = (xq*xAng + yq*yAng + zq*zAng>(-d_l/2.)+x0_l*xAng + y0_l*yAng + z0_l*zAng);
                int conditionLR = (fabs(xq-x0_l)<thickness/2.);
                
                if(conditionP && conditionM && conditionLR) {
                    materialType[IJKf(xi,yi,zi)] = 'p';
                }
            }
        }
    }
    
}

void makeWaveplate() {
    
    makeplanar(90.,0,domainX/2,0,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0003,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0006,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0009,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0012,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0015,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0018,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0021,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0024,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0027,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0030,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0033,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0036,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0039,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0042,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0045,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0048,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0051,domainX/2,0.00012,0.004);
    makeplanar(90.,0,domainX/2,0.0054,domainX/2,0.00012,0.004);
}

//_______Near To Far Field______
FLOAT ExTminusXmax[GRIDYZ];
FLOAT EyTminusXmax[GRIDYZ];
FLOAT EzTminusXmax[GRIDYZ];

FLOAT HxTminusXmax[GRIDYZ];
FLOAT HyTminusXmax[GRIDYZ];
FLOAT HzTminusXmax[GRIDYZ];

FLOAT dEx[TSTEPS+10]={0};
FLOAT dEy[TSTEPS+10]={0};
FLOAT dEz[TSTEPS+10]={0};
FLOAT dHx[TSTEPS+10]={0};
FLOAT dHy[TSTEPS+10]={0};
FLOAT dHz[TSTEPS+10]={0};

void initdEH() {
    for (int n=0; n<TSTEPS+10; n++) {
        dEx[n]=0;
        dEy[n]=0;
        dEz[n]=0;
        dHx[n]=0;
        dHy[n]=0;
        dHz[n]=0;
    }
}

void addToPotentialFields(int n, double l) {
    double r = (domainX/2.)+l;
    double f = ((r*r) - domainX/2.)/(speedC*deltaT);
    double nnFloatE = n+0.5+f;
    int nnE = nnFloatE;
    double aE = nnFloatE-nnE;
    
    double nnFloatH = n+f;
    int nnH = nnFloatH;
    double aH = nnFloatH-nnH;
    
    double prefactor = speedC*deltaT*CouS/(4*M_PI*r);
    for (int j=0; j<GRIDY; j++) {
        for (int k=0; k<GRIDZ; k++) {
            int IJKmax=IJKf(GRIDX-2, j, k);
            int wallJK=IJKf(0, j, k);
            
            dEy[nnE-1]+=(1-aE)*prefactor*(Ey[IJKmax]-EyTminusXmax[wallJK]);
            dEy[nnE] += aE*prefactor*(Ey[IJKmax]-EyTminusXmax[wallJK]);
            dEz[nnE-1]+=(1-aE)*prefactor*(Ez[IJKmax]-EzTminusXmax[wallJK]);
            dEz[nnE] += aE*prefactor*(Ez[IJKmax]-EzTminusXmax[wallJK]);
            dEx[nnE-1]+=(1-aE)*prefactor*(Ex[IJKmax]-ExTminusXmax[wallJK]);
            dEx[nnE] += aE*prefactor*(Ex[IJKmax]-ExTminusXmax[wallJK]);
            
            dHy[nnH]-=(1-aH)*prefactor*(Hy[IJKmax]-HyTminusXmax[wallJK]);
            dHy[nnH-1] -= aH*prefactor*(Hy[IJKmax]-HyTminusXmax[wallJK]);
            dHz[nnH]-=(1-aH)*prefactor*(Hz[IJKmax]-HzTminusXmax[wallJK]);
            dHz[nnH-1] -= aH*prefactor*(Hz[IJKmax]-HzTminusXmax[wallJK]);
            dHx[nnH]-=(1-aH)*prefactor*(Hx[IJKmax]-HxTminusXmax[wallJK]);
            dHx[nnH-1] -= aH*prefactor*(Hx[IJKmax]-HxTminusXmax[wallJK]);
            
            ExTminusXmax[wallJK]=Ex[IJKmax];
            EyTminusXmax[wallJK]=Ey[IJKmax];
            EzTminusXmax[wallJK]=Ez[IJKmax];
            HxTminusXmax[wallJK]=Hx[IJKmax];
            HyTminusXmax[wallJK]=Hy[IJKmax];
            HzTminusXmax[wallJK]=Hz[IJKmax];
        }
    }
}

void saveField(FILE * file, FLOAT * field) {
    for (int n=0; n<TSTEPS; n++) {
        fprintf(file, "%.8f; ", field[n]);
    }
}



int main(int argc, const char * argv[]) {
    double tOld = current_timestamp();
    
    ExReferenceFile = fopen(ExReferenceFilename,"w");
    EyReferenceFile = fopen(EyReferenceFilename,"w");
    EzReferenceFile = fopen(EzReferenceFilename,"w");
    ExSampleFile = fopen(ExSampleFilename,"w");
    EySampleFile = fopen(EySampleFilename,"w");
    EzSampleFile = fopen(EzSampleFilename,"w");
    
    materialFieldOut = fopen(materialFile, "w");
    
    //ExPtReferenceFile = fopen(ExPtReferenceFilename,"w");
    //EyPtReferenceFile = fopen(EyPtReferenceFilename,"w");
    //EzPtReferenceFile = fopen(EzPtReferenceFilename,"w");
    //ExPtSampleFile = fopen(ExPtSampleFilename,"w");
    //EyPtSampleFile = fopen(EyPtSampleFilename,"w");
    //EzPtSampleFile = fopen(EzPtSampleFilename,"w");
    
    dExReferenceFile = fopen(dExReferenceFileName,"w");
    dEyReferenceFile = fopen(dEyReferenceFileName,"w");
    dEzReferenceFile = fopen(dEzReferenceFileName,"w");
    dHxReferenceFile = fopen(dHxReferenceFileName,"w");
    dHyReferenceFile = fopen(dHyReferenceFileName,"w");
    dHzReferenceFile = fopen(dHzReferenceFileName,"w");
    
    dExSampleFile = fopen(dExSampleFileName,"w");
    dEySampleFile = fopen(dEySampleFileName,"w");
    dEzSampleFile = fopen(dEzSampleFileName,"w");
    dHxSampleFile = fopen(dHxSampleFileName,"w");
    dHySampleFile = fopen(dHySampleFileName,"w");
    dHzSampleFile = fopen(dHzSampleFileName,"w");
    
    getBounds();
    getConsts(C_vac, vac);
    getConsts(C_paper, paper);
    
    // REFERENCE RUN
    
    initArraysToZero(); initABC(); makesource(); initdEH();
    
    printf("size of (t; x, y, z) = (%i; %i, %i, %i)\n", TSTEPS,GRIDX,GRIDY,GRIDZ);
    for (int n=0; n<TSTEPS; n++) {
        if(n%10==0) { printf("Reference: %i/%i\n",n/10,TSTEPS/10); }
        source(n);
        E_CPU_para();
        H_CPU_para();
        addToPotentialFields(n, 5*domainX);
        //E_CPU_seq();
        //H_CPU_seq();
        ABC(1,1,1,1,1,1);
        writeToFile(ExReferenceFile,Ex);
        writeToFile(EyReferenceFile,Ey);
        writeToFile(EzReferenceFile,Ez);
        
        //writePointToFile(ExPtReferenceFile,Ex);
        //writePointToFile(EyPtReferenceFile,Ey);
        //writePointToFile(EzPtReferenceFile,Ez);
        
        //takeWallToSingleValueAndWriteToFile(ExPtFSOut,Ex);
        //takeWallToSingleValueAndWriteToFile(EyPtFSOut,Ey);
        //takeWallToSingleValueAndWriteToFile(EzPtFSOut,Ez);
        
        
    }
    
    saveField(dExReferenceFile, dEx);
    saveField(dEyReferenceFile, dEy);
    saveField(dEzReferenceFile, dEz);
    saveField(dHxReferenceFile, dHx);
    saveField(dHyReferenceFile, dHy);
    saveField(dHzReferenceFile, dHz);
    
    // SAMPLE RUN
    
    initArraysToZero(); initABC(); makesource(); initdEH();
    makeWaveplate();
    writeMaterialFile(materialFieldOut, materialType);
    
    
    
    
    for (int n=0; n<TSTEPS; n++) {
        if(n%10==0) { printf("Sample: %i/%i\n",n/10,TSTEPS/10); }
        source(n);
        E_CPU_para();
        H_CPU_para();
        addToPotentialFields(n, 5*domainX);
        //E_CPU_seq();
        //H_CPU_seq();
        ABC(1,1,1,1,1,1);
        
        writeToFile(ExSampleFile,Ex);
        writeToFile(EySampleFile,Ey);
        writeToFile(EzSampleFile,Ez);
        
        
        //writePointToFile(ExPtSampleFile,Ex);
        //writePointToFile(EyPtSampleFile,Ey);
        //writePointToFile(EzPtSampleFile,Ez);
        
        //takeWallToSingleValueAndWriteToFile(ExPtWPOut,Ex);
        //takeWallToSingleValueAndWriteToFile(EyPtWPOut,Ey);
        //takeWallToSingleValueAndWriteToFile(EzPtWPOut,Ez);
    }
    
    saveField(dExSampleFile, dEx);
    saveField(dEySampleFile, dEy);
    saveField(dEzSampleFile, dEz);
    saveField(dHxSampleFile, dHx);
    saveField(dHySampleFile, dHy);
    saveField(dHzSampleFile, dHz);
    
    double tNew = current_timestamp();
    printf("%.2fs\n", (tNew-tOld)*0.001);
    
    return 0;
}