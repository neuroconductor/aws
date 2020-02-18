#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#include <Rinternals.h> // for SEXP
#include <R_ext/RS.h>
void F77_NAME(adsmse3m)( double* y, double* th, double* ni, double* sthi,
  int* pos, int* nv, int* ns, int* n1, int* n2, int* n3, int* ngrad, double* lambda,
  int* ncores, int* ind, double* w, int* n, double* thn, double* sw,
  double* swy, double* si, double* thi);
void F77_NAME(adsmse3p)( double* y, double* th, double* ni, int* pos, int* nv,
  int* n1, int* n2, int* n3, int* ngrad, double* lambda, int* ncoils,
  int* ncores, int* ind, double* w, int* n, double* thn, double* ldf,
  double* sw, double* swy, int* model);
void F77_NAME(adsmse3s)( double* y, double* y0, double* th, double* ni,
  double* th0, double* ni0, double* fsi2, double* fsi02, int* pos, int* nv, int* ns,
  int* n1, int* n2, int* n3, int* ngrad, double* lambda, double* ws0,
  int* ind, double* w, int* n, int* ind0, double* w0, int* n0, double* thn,
  double* nin, double* th0n, double* ni0n, double* sw, double* swy,
  double* thi, double* nii, double* fsi2i);
void F77_NAME(afmodem1)( double* y, int* n1, int* n2, int* n3, int* mask,
  double* h, double* vext, double* sigma);
void F77_NAME(afmodem2)( double* y, int* n1, int* n2, int* n3, int* mask,
  double* h, double* vext, double* sm);
void F77_NAME(afmodevn)(  double* y, int* n1, int* n2, int* n3, int* mask,
  double* h, double* vext, double* sigma);
void F77_NAME(asmse30p)( double* y, double* th, double* ni, int* pos, int* nv,
  int* n1, int* n2, int* n3, double* lambda, int* ncoils, int* ind, double* w,
  int* n, int* starts, int* nstarts, double* thn, double* ldf, double* swi,
  int* model);
void F77_NAME(awsadchi)( double* y, double* th, double* ni, double* fns,
  int* mask, int* n1, int* n2, int* n3, int* ind, double* w, int* nw,
  double* lambda, double* sigma, double* wad, int* nthreds, double* thn,
  double* sy);
void F77_NAME(awslchi2)( double* s, double* ksi, double* ni, double* sigma,
  double* vpar, double* L, int* mask, int* n1, int* n2, int* n3, int* ind,
  double* w, int* nw, double* minni, double* wad, double* sad, double* lambda,
  int* nthreds, int* iL, double* work, double* thn, double* sigman,
  double* ksin, double* flb, int* nfb);
void F77_NAME(awslgaus)( double* s, double* th, double* ni, double* sigma,
	int* mask, int* n1, int* n2, int* n3, int* ind, double* w, int* nw,
	double* minni, double* lambda, double* thn, double* sigman);
void F77_NAME(awsph1)(double* y, double* si, int* fix, int* nfix,
  int* n, int* degr, double* hw, double* hakt, double* hhom,
  double* lambda, double* theta, double* bi, double* bi2,
  double* bi0, double* ai, int* kern, double* spmin, double* lw,
  double* w, double* slw, double* sw, int* ind);
void F77_NAME(awsph2)(double* y, double* si, int* fix, int* nfix,
  int* n1, int* n2, int* degr, double* hw, double* hakt,
  double* hhom, double* lambda, double* theta, double* bi,
  double* bi2, double* bi0, double* ai, int* kern, double* spmin,
  double* lw, double* w, double* slw, double* sw, int* ind);
void F77_NAME(awsp1b)(double* y, int* fix, int* nfix, int* n,
  int* degr, double* hw, double* hakt, double* hhom, double* lambda,
  double* theta, double* bi, double* bi2, double* bi0, double* ai,
  int* kern, double* spmin, double* lw, double* w, double* slw,
  double* sw, int* ind);
void F77_NAME(awsp2)(double* y, int* fix, int* nfix, int* n1,
  int* n2, int* degr, double* hw, double* hakt, double* hhom, double* lambda,
  double* theta, double* bi, double* bi2, double* bi0, double* ai,
  int* kern, double* spmin, double* lw, double* w, double* slw,
  double* sw, int* ind);
void F77_NAME(awsvchi)( double* y, double* th, double* ni, double* fns,
  int* mask, int* n1, int* n2, int* n3, int* ind, double* w, int* nw,
  double* lambda, double* sigma, double* thn, double* sy);
void F77_NAME(bgstats)( double* g, int* n, double* bg, double* bghat);
void F77_NAME(caws)(double* y, int* pos, int* n1, int* n2, int* n3,
  double* hakt, double* lambda, double* theta,
  double* bi, double* bi2, double* bi0, double* ai, int* model,
  int* kern, double* spmin, double* lwght, double* wght);
void F77_NAME(cawsmask)(double* y, int* mask, int* ni, int* fix,
  int* n1, int* n2, double* hakt,
  double* bi, double* bi2, double* bi0, double* ai,
  int* kern, double* lwght, double* wght);
void F77_NAME(cawsw)(int* n1, int* n2, int* n3, double* hakt,
  double* lambda, double* theta, double* bi, int* model, int* kern,
  double* spmin, double* lwght, double* wght);
void F77_NAME(cawsw1)(int* n1, int* n2, int* n3, int* inx, int* iny,
  int* inz, int* anz, double* hakt, double* lambda, double* theta,
  double* bi, int* model, int* kern, double* spmin, double* lwght,
  double* wght);
void F77_NAME(caws6)(double* y, int* pos, int* n1, int* n2, int* n3,
  double* hakt, double* lambda, double* theta, double* fnc,
  double* bi, double* bi2, double* bi0, double* ai, int* kern,
  double* spmin, double* lwght, double* wght);
void F77_NAME(cgaws)(double* y, int* pos, double* si2,
  int* n1, int* n2, int* n3, double* hakt,
  double* lambda, double* theta, double* bi, double* bi2,
  double* bi0, double* gi, double* gi2, double* ai, int* kern,
  double* spmin, double* lwght, double* wght);
void F77_NAME(cgawsmas)(double* y, int* mask, int* ni, int* fix,
  double* si2, int* n1, int* n2, double* hakt, double* lambda,
  double* theta, double* bi, double* bi2, double* bi0, double* vred,
  double* ai, int* model, int* kern, double* spmin, double* lwght,
  double* wght);
void F77_NAME(chaws)(double* y, double* si2, int* pos, int* n1,
  int* n2, int* n3, double* hakt, double* lambda, double* theta,
  double* bi, double* bi2, double* bi0, double* ai,
  int* model, int* kern, double* spmin, double* lwght, double* wght);
void F77_NAME(chaws2)(double* y, double* si2, int* pos, int* wlse, int* n1,
	int* n2, int* n3, double* hakt, double* lambda,
	double* theta, double* bi, double* thn, int* kern, int* skern,
	double* spmin, double* spmax, double* lwght, double* wght);
void F77_NAME(chawsv)(double* y, double* res, double* si2, int* pos, int* wlse,
	int* n1, int* n2, int* n3, int* n4, double* hakt,
	double* lambda, double* theta, double* bi, double* resnew,
	double* thn, int* kern, int* skern, double* spmin,
	double* spmax, double* lwght, double* wght, double* resi);
void F77_NAME(exceed)(double* x, int* n, double* z, int* nz,
  double* exprob);
void F77_NAME(exceedm)(double* x, int* n, double* z, int* nz,
    double* exprob, int* mask);
void F77_NAME(gethani)(double* x, double* y, int* kern,
  double* value, double* wght, double* eps, double* bw);
void F77_NAME(getmsni0)( double* ni, int* n, int* lindi, double* msni);
void F77_NAME(getmsth0)( double* theta, int* n, int* lindi, double* msth);
void F77_NAME(getvofh)(double* bw, int* kern, double* wght,
  double* vol);
void F77_NAME(ghfse3i)( int* i4, int* kstar, double* k456, int* ng,
	double* kappa, double* vext, double* h, double* varred, int* n, int* dist);
void F77_NAME(hg1f1)( double* a, double* b, double* z, int* n, double* fz);
void F77_NAME(ihaws2)(double* y, double* si2, int* pos, int* wlse, int* n1,
  int* n2, int* n3, int* dv, double* hakt, double* lambda, double* theta,
	int* ncores, double* bi, double* thn, int* kern, int* skern, double* spmin,
	double* spmax, double* lwght, double* wght, double* swjy);
void F77_NAME(ipolsp)( double* theta, double* th0, double* ni, double* ni0,
  int* n, int* ng, int* gind, double* gw, int* nbv, int* nbvp1, double* msth,
  double* msni);
void F77_NAME(imcorr)(double* res, int* mask, int* n1, int* n2, int* n3, int* nv,
	double* scorr, int* l1, int* l2, int* l3);
void F77_NAME(ipolsp1)( double* theta, double* th0, double* ni, double* ni0,
  int* n, int* ng, int* gind, double* gw, int* nbv, int* nbvp1,
  double* msth, double* msni);
void F77_NAME(ivar)(double* res, double* resscale, int* nvoxel, int* nt, double* var);
void F77_NAME(k456krb)( double* par, double* b, double* matm, double* erg);
void F77_NAME(lkern1)(double* x, int* n, double* h, int* kern,
  int* m, double* khx);
void F77_NAME(lkfuls0)( double* h, double* vext, int* ind, double* wght,
  int* n);
void F77_NAME(lkfulse3)( double* h, double* kappa, double* k456, int* ng,
  double* vext, int* ind, double* wght, int* n, int* dist);
void F77_NAME(mask)(int* maskin, int* maskout, int* n1, int* n2, int* h);
void F77_NAME(median1d)(double* y, int* n, double* yhat);
void F77_NAME(median2d)(double* y, int* n1, int* n2, double* yhat);
void F77_NAME(median3d)(double* y, int* n1, int* n2, int* n3,
  double* yhat);
void F77_NAME(mediansm)(double* y, int* mask, int* n1, int* n2, int* n3,
  int* ind, int* nind, double* work, int* ncores, double* yout);
void F77_NAME(mpaws1)(int* n, int* dp1, int* dp2, double* ai,
  double* bi, double* theta, double* dmat, int* ind);
void F77_NAME(mpaws2)(int* n, int* dp1, int* dp2, double* ai,
  double* bi, double* theta, double* dmat, int* ind);
void F77_NAME(paramw3)(double* h, double* vext, int* indn, double* w, int* n);
void F77_NAME(pawswght)(int* n1, int* n2, int* n3, int* i1, int* i2, int* i3,
  double* hakt, double* lambda, double* theta, double* bi,
  int* model, int* kern, double* spmin, double* lwght, double* wght,
  int* npsize, double* wi);
void F77_NAME(pcaws)(double* y, int* pos, int* n1, int* n2, int* n3,
  double* hakt, double* lambda, double* theta, double* bi,
  double* bi2, double* bin, double* thnew, int* model, int* kern, double* spmin,
  double* lwght, double* wght, int* npsize);
void F77_NAME(pvaws)(double* y, int* pos, int* nv, int* n1, int* n2,
  int* n3, double* hakt, double* lambda, double* theta, double* bi,
  double* bin, double* thnew, int* ncores, double* spmin, double* lwght,
  double* wght, double* swjy, int* np1, int* np2, int* np3);
void F77_NAME(pvaws2)(double* y, int* pos, int* nv, int* nvd,
  int* n1, int* n2, int* n3, double* hakt, double* lambda,
  double* theta, double* bi, double*bin, double* thnew, double* invcov,
  int* ncores, double* spmin, double* lwght, double* wght,
  double* swjy, int* np1, int* np2, int* np3);
void F77_NAME(pvawsme)(double* y, double* yd, int* pos, int* nv, int* nvd, int* nd, int* n1, int* n2,
  int* n3, double* hakt, double* lambda, double* theta, double* bi, double* bin,
  double* thnew, double* ydnew, double* invcov, int* ncores, double* spmin, double* lwght,
  double* wght, double* swjy, double* swjd, int* np1, int* np2, int* np3);
void F77_NAME(sector)(double* x1, int* n1, double* x2, int* n2,
  int* nsect, int* sect, int* symm, double* insect);
void F77_NAME(segment)(double* y, int* pos, int* fix, double* level,
  double* delta, double* si2, int* n1, int* n2, int* n3, double* hakt,
  double* lambda, double* theta, double* bi, double* bi2,
  double* bi0, double* gi, double* vred, double* thetan, int* kern,
  double* spmin, double* lwght, double* wght, int* segm, int* segmn,
  double* beta, double* thresh, double* ext, double* fov, double* varest);
void F77_NAME(smooth3d)(double* y, double* si2, int* mask, int* wlse, int* nvox,
	int* n1, int* n2, int* n3, int* dv, double* hakt, double* thn, double* bi,
	int* kern, double* lwght, double* wght, double* swjy);
void F77_NAME(vaws)(double* y, int* pos, int* nv, int* n1, int* n2,
  int* n3, double* hakt, double* lambda, double* theta, double* bi,
  double* vred, double* thnew, int* ncores, double* spmin, double* lwght,
  double* wght, double* swjy);
void F77_NAME(vaws2cov)(double* y, int* pos, int* nv, int* nvd, int* n1,
  int* n2, int* n3, double* hakt, double* lambda, double* theta, double* bi,
  double* vred, double* thnew, double* invcov, int* ncores, double* spmin,
  double* lwght, double* wght, double* swjy, double* thi, double* invcovi);
void F77_NAME(cplxawss)(double* y, int* mask, int* nv, int* n1, int* n2, int* n3,
  double* hakt, double* lambda, double* theta, double* s2, double* bi,
  double* thnew, double* s2new, int* ncores, double* lwght, double* wght,
  double* swjy);
void F77_NAME(vpaws)(int* n, int* dp2, double* bi, double* bi2, double* var);
void F77_NAME(fillpat1)(double* x, int* n1, int* phw, int* psize, double* pmat);
void F77_NAME(fillpat2)(double* x, int* n1, int* n2, int* phw, int* psize,
  double* pmat);
void F77_NAME(fillpat3)(double* x, int* n1, int* n2, int* n3, int* phw,
  int* psize, double* pmat);
void F77_NAME(nlmeans)(double* x, int* n1, int* n2, int* n3, double* patch,
  int* pd, int* swd, double* tau, double* xhat);
static R_NativePrimitiveArgType adsmse3m_t[]={REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP,
	INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType adsmse3p_t[]={REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType adsmse3s_t[]={REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP};
static R_NativePrimitiveArgType afmodem1_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType afmodem2_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType afmodevn_t[]={REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType asmse30p_t[]={REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP,
  REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType awsadchi_t[]={REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType awslchi2_t[]={REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType awslgaus_t[]={REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP,
	REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType awsph1_t[]={REALSXP, REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType awsph2_t[]={REALSXP, REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,REALSXP, REALSXP, REALSXP, REALSXP,
  INTSXP};
static R_NativePrimitiveArgType awsp1b_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType awsp2_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType awsvchi_t[]={REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType bgstats_t[]={REALSXP, INTSXP, REALSXP,
  REALSXP};
static R_NativePrimitiveArgType caws_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType cawsmask_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType cawsw_t[]={INTSXP, INTSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType cawsw1_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType caws6_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType cgaws_t[]={REALSXP, INTSXP, REALSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType cgawsmas_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType chaws_t[]={REALSXP, REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType chaws2_t[]={REALSXP, REALSXP, INTSXP,
	INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType chawsv_t[]={REALSXP, REALSXP, REALSXP,
	INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType exceed_t[]={REALSXP, INTSXP, REALSXP, INTSXP,
  REALSXP};
static R_NativePrimitiveArgType exceedm_t[]={REALSXP, INTSXP, REALSXP, INTSXP,
  REALSXP, INTSXP};
static R_NativePrimitiveArgType gethani_t[]={REALSXP, REALSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType getmsni0_t[]={REALSXP, INTSXP, INTSXP,
  REALSXP};
static R_NativePrimitiveArgType getmsth0_t[]={REALSXP, INTSXP, INTSXP,
	 REALSXP};
static R_NativePrimitiveArgType getvofh_t[]={REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType ghfse3i_t[]={INTSXP, INTSXP, REALSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType hg1f1_t[]={REALSXP, REALSXP, REALSXP,
  INTSXP, REALSXP};
static R_NativePrimitiveArgType ihaws2_t[]={REALSXP, REALSXP, INTSXP, INTSXP,
	INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
	REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType imcorr_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
	INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType ipolsp_t[]={REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType ipolsp1_t[]={REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, INTSXP, REALSXP,
  REALSXP};
static R_NativePrimitiveArgType ivar_t[]={REALSXP, REALSXP, INTSXP, INTSXP,
  REALSXP};
static R_NativePrimitiveArgType k456krb_t[]={REALSXP, REALSXP,
  REALSXP, REALSXP};
static R_NativePrimitiveArgType lkern1_t[]={REALSXP, INTSXP, REALSXP, INTSXP,
  INTSXP, REALSXP};
static R_NativePrimitiveArgType lkfuls0_t[]={REALSXP, REALSXP,
  INTSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType lkfulse3_t[]={REALSXP, REALSXP, REALSXP,
  INTSXP, REALSXP, INTSXP, REALSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType mask_t[]={INTSXP, INTSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType median1d_t[]={REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType median2d_t[]={REALSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType median3d_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  REALSXP};
static R_NativePrimitiveArgType mediansm_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType mpaws1_t[]={INTSXP, INTSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType mpaws2_t[]={INTSXP, INTSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType paramw3_t[]={REALSXP, REALSXP, INTSXP, REALSXP,
  INTSXP};
static R_NativePrimitiveArgType pawswght_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType pcaws_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType pvaws_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType pvaws2_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType pvawsme_t[]={REALSXP, REALSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType sector_t[]={REALSXP, INTSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType segment_t[]={REALSXP, INTSXP, INTSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType smooth3d_t[]={REALSXP, REALSXP, INTSXP, INTSXP,
	INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP,
	REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType vaws_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType vaws2cov_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP};
static R_NativePrimitiveArgType cplxawss_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType vpaws_t[]={INTSXP, INTSXP, REALSXP, REALSXP,
  REALSXP};
static R_NativePrimitiveArgType fillpat1_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  REALSXP};
static R_NativePrimitiveArgType fillpat2_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, REALSXP};
static R_NativePrimitiveArgType fillpat3_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType nlmeans_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, INTSXP, INTSXP, REALSXP, REALSXP};

static const R_FortranMethodDef fmethods[] = {
						{"adsmse3m", (DL_FUNC) &adsmse3m_ , 21, adsmse3m_t},
						{"adsmse3p", (DL_FUNC) &adsmse3p_ , 20, adsmse3p_t},
						{"adsmse3s", (DL_FUNC) &adsmse3s_ , 32, adsmse3s_t},
						{"afmodem1", (DL_FUNC) &afmodem1_ , 8, afmodem1_t},
						{"afmodem2", (DL_FUNC) &afmodem2_ , 8, afmodem2_t},
						{"afmodevn", (DL_FUNC) &afmodevn_ , 8, afmodevn_t},
						{"asmse30p", (DL_FUNC) &asmse30p_ , 19, asmse30p_t},
						{"awsadchi", (DL_FUNC) &awsadchi_ , 17, awsadchi_t},
						{"awslchi2", (DL_FUNC) &awslchi2_ , 25, awslchi2_t},
						{"awslgaus", (DL_FUNC) &awslgaus_ , 15, awslgaus_t},
            {"awsph1", (DL_FUNC) &awsph1_ , 22, awsph1_t},
            {"awsph2", (DL_FUNC) &awsph2_ ,23, awsph2_t},
            {"awsp1b", (DL_FUNC) &awsp1b_ ,21, awsp1b_t},
            {"awsp2", (DL_FUNC) &awsp2_ ,22, awsp2_t},
						{"awsvchi", (DL_FUNC) &awsvchi_ , 15, awsvchi_t},
						{"bgstats", (DL_FUNC) &bgstats_ , 4, bgstats_t},
            {"caws", (DL_FUNC) &caws_ ,17, caws_t},
            {"cawsmask", (DL_FUNC) &cawsmask_ ,14, cawsmask_t},
            {"cawsw", (DL_FUNC) &cawsw_ ,12, cawsw_t},
            {"cawsw1", (DL_FUNC) &cawsw1_ ,16, cawsw1_t},
            {"caws6", (DL_FUNC) &caws6_ ,17, caws6_t},
            {"cgaws", (DL_FUNC) &cgaws_ ,19, cgaws_t},
            {"cgawsmas", (DL_FUNC) &cgawsmas_ ,20, cgawsmas_t},
            {"chaws", (DL_FUNC) &chaws_ ,18, chaws_t},
						{"chaws2", (DL_FUNC) &chaws2_ ,18, chaws2_t},
            {"chawsv", (DL_FUNC) &chawsv_ ,22, chawsv_t},
						{"exceed", (DL_FUNC) &exceed_ ,5, exceed_t},
            {"exceedm", (DL_FUNC) &exceedm_ ,6, exceedm_t},
            {"gethani", (DL_FUNC) &gethani_ ,7, gethani_t},
						{"getmsni0", (DL_FUNC) &getmsni0_ , 4, getmsni0_t},
						{"getmsth0", (DL_FUNC) &getmsth0_ , 4, getmsth0_t},
            {"getvofh", (DL_FUNC) &getvofh_ ,4, getvofh_t},
						{"ghfse3i", (DL_FUNC) &ghfse3i_ , 10, ghfse3i_t},
						{"hg1f1", (DL_FUNC) &hg1f1_ , 5, hg1f1_t},
						{"ihaws2", (DL_FUNC) &ihaws2_ ,21, ihaws2_t},
						{"imcorr", (DL_FUNC) &imcorr_ ,10, imcorr_t},
						{"ipolsp", (DL_FUNC) &ipolsp_ , 12, ipolsp_t},
            {"ipolsp1", (DL_FUNC) &ipolsp1_ , 12, ipolsp1_t},
						{"ivar", (DL_FUNC) &ivar_ ,5, ivar_t},
            {"k456krb", (DL_FUNC) &k456krb_ , 4, k456krb_t},
            {"lkern1", (DL_FUNC) &lkern1_ ,6, lkern1_t},
						{"lkfuls0", (DL_FUNC) &lkfuls0_ , 5, lkfuls0_t},
            {"lkfulse3", (DL_FUNC) &lkfulse3_ , 9, lkfulse3_t},
            {"mask", (DL_FUNC) &mask_ ,5, mask_t},
            {"median1d", (DL_FUNC) &median1d_ ,3, median1d_t},
            {"median2d", (DL_FUNC) &median2d_ ,4, median2d_t},
            {"median3d", (DL_FUNC) &median3d_ ,5, median3d_t},
            {"mediansm", (DL_FUNC) &mediansm_ , 10, mediansm_t},
            {"mpaws1", (DL_FUNC) &mpaws1_ ,8, mpaws1_t},
            {"mpaws2", (DL_FUNC) &mpaws2_ ,8, mpaws2_t},
            {"paramw3", (DL_FUNC) &paramw3_ , 5, paramw3_t},
            {"pawswght", (DL_FUNC) &pawswght_ ,17, pawswght_t},
            {"pcaws", (DL_FUNC) &pcaws_ ,18, pcaws_t},
            {"pvaws", (DL_FUNC) &pvaws_ ,20, pvaws_t},
            {"pvaws2", (DL_FUNC) &pvaws2_ ,22, pvaws2_t},
            {"pvawsme", (DL_FUNC) &pvawsme_ , 26, pvawsme_t},
            {"sector", (DL_FUNC) &sector_ ,8, sector_t},
            {"segment", (DL_FUNC) &segment_ ,29, segment_t},
						{"smooth3d", (DL_FUNC) &smooth3d_ ,16, smooth3d_t},
            {"vaws", (DL_FUNC) &vaws_ ,17, vaws_t},
            {"vaws2cov", (DL_FUNC) &vaws2cov_ ,21, vaws2cov_t},
            {"cplxawss", (DL_FUNC) &cplxawss_ ,17, cplxawss_t},
            {"vpaws", (DL_FUNC) &vpaws_ ,5, vpaws_t},
            {"fillpat1", (DL_FUNC) &fillpat1_ ,5, fillpat1_t},
            {"fillpat2", (DL_FUNC) &fillpat2_ ,6, fillpat2_t},
            {"fillpat3", (DL_FUNC) &fillpat3_ ,7, fillpat3_t},
            {"nlmeans", (DL_FUNC) &nlmeans_ ,9, nlmeans_t},
            {NULL, NULL, 0,NULL}
};

void R_init_aws(DllInfo *dll)
         {
             R_registerRoutines(dll, NULL, NULL, fmethods , NULL);
             R_useDynamicSymbols(dll,FALSE);
             R_forceSymbols(dll,TRUE);
         }
