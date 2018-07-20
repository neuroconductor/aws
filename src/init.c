#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}
#include <Rinternals.h> // for SEXP
#include <R_ext/RS.h>
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
void F77_NAME(caws)(double* y, int* fix, int* n1, int* n2, int* n3,
  double* hakt, double* hhom, double* lambda, double* theta,
  double* bi, double* bi2, double* bi0, double* ai, int* model,
  int* kern, double* spmin, double* lwght, double* wght);
void F77_NAME(cawsmask)(double* y, int* mask, int* ni, int* fix,
  int* n1, int* n2, double* hakt, double* lambda, double* theta,
  double* bi, double* bi2, double* bi0, double* ai, int* model,
  int* kern, double* spmin, double* lwght, double* wght);
void F77_NAME(cawsw)(int* n1, int* n2, int* n3, double* hakt,
  double* lambda, double* theta, double* bi, int* model, int* kern,
  double* spmin, double* lwght, double* wght);
void F77_NAME(cawsw1)(int* n1, int* n2, int* n3, int* inx, int* iny,
  int* inz, int* anz, double* hakt, double* lambda, double* theta,
  double* bi, int* model, int* kern, double* spmin, double* lwght,
  double* wght);
void F77_NAME(caws1)(double* y, int* n1, int* n2, int* n3,
  double* hakt, double* bi, double* bi2, double* bi0, double* ai,
  int* kern, double* lwght, double* wght);
void F77_NAME(caws6)(double* y, int* fix, int* n1, int* n2, int* n3,
  double* hakt, double* hhom, double* lambda, double* theta, double* fnc,
  double* bi, double* bi2, double* bi0, double* ai, int* kern,
  double* spmin, double* lwght, double* wght);
void F77_NAME(cgaws)(double* y, int* fix, int* mask, double* si2,
  int* n1, int* n2, int* n3, double* hakt, double* hhom,
  double* lambda, double* theta, double* bi, double* bi2,
  double* bi0, double* gi, double* vred, double* ai, int* kern,
  double* spmin, double* lwght, double* wght);
void F77_NAME(cgawsmas)(double* y, int* mask, int* ni, int* fix,
  double* si2, int* n1, int* n2, double* hakt, double* lambda,
  double* theta, double* bi, double* bi2, double* bi0, double* vred,
  double* ai, int* model, int* kern, double* spmin, double* lwght,
  double* wght);
void F77_NAME(chaws)(double* y, int* fix, double* si2, int* n1,
  int* n2, int* n3, double* hakt, double* lambda, double* theta,
  double* bi, double* bi2, double* bi0, double* vred, double* ai,
  int* model, int* kern, double* spmin, double* lwght, double* wght);
void F77_NAME(chaws1)(double* y, double* si2, int* n1, int* n2,
  int* n3, double* hakt, double* bi, double* bi2, double* bi0,
  double* vred, double* ai, int* kern, double* lwght, double* wght);
void F77_NAME(exceed)(double* x, int* n, double* z, int* nz,
  double* exprob);
void F77_NAME(gethani)(double* x, double* y, int* kern,
  double* value, double* wght, double* eps, double* bw);
void F77_NAME(getvofh)(double* bw, int* kern, double* wght,
  double* vol);
void F77_NAME(lkern1)(double* x, int* n, double* h, int* kern,
  int* m, double* khx);
void F77_NAME(mask)(int* maskin, int* maskout, int* n1, int* n2, int* h);
void F77_NAME(median1d)(double* y, int* n, double* yhat);
void F77_NAME(median2d)(double* y, int* n1, int* n2, double* yhat);
void F77_NAME(median3d)(double* y, int* n1, int* n2, int* n3,
  double* yhat);
void F77_NAME(mpaws1)(int* n, int* dp1, int* dp2, double* ai,
  double* bi, double* theta, double* dmat, int* ind);
void F77_NAME(mpaws2)(int* n, int* dp1, int* dp2, double* ai,
  double* bi, double* theta, double* dmat, int* ind);
void F77_NAME(pawswght)(int* n1, int* n2, int* n3, int* i1, int* i2, int* i3,
  double* hakt, double* lambda, double* theta, double* bi,
  int* model, int* kern, double* spmin, double* lwght, double* wght,
  int* npsize, double* wi);
void F77_NAME(pcaws)(double* y, int* n1, int* n2, int* n3,
  double* hakt, double* lambda, double* theta, double* bi,
  double* bi2, double* bi0, double(bin), double* ai, int* model, int* kern,
  double* spmin, double* lwght, double* wght, int* npsize);
void F77_NAME(pcaws2)(double* y, int* n1, int* n2, int* n3,
  double* hakt, double* lambda, double* theta, double* bi, double* bip,
  double* bi2, double* bi0, double(bin), double* ai, int* model, int* kern,
  double* spmin, double* lwght, double* wght, int* npsize);
void F77_NAME(pcaws3)(double* y, int* n1, int* n2, int* n3,
    double* hakt, double* lambda, double* theta, double* bi,
    double* bi2, double* bi0, double(bin), double* ai, int* model, int* kern,
    double* spmin, double* lwght, double* wght, int* npsize, int* qind);
void F77_NAME(pcawsm)(double* y, int* pos, int* n1, int* n2, int* n3,
  double* hakt, double* lambda, double* theta, double* bi,
  double* bi2, double* bin, double* thnew, int* model, int* kern, double* spmin,
  double* lwght, double* wght, int* npsize);
void F77_NAME(pvaws)(double* y, int* mask, int* nv, int* n1, int* n2,
  int* n3, double* hakt, double* lambda, double* theta, double* bi,
  double* bin, double* thnew, int* ncores, double* spmin, double* lwght,
  double* wght, double* swjy, int* np1, int* np2, int* np3);
void F77_NAME(pvaws2)(double* y, int* mask, int* nv, int* nvd,
  int* n1, int* n2, int* n3, double* hakt, double* lambda,
  double* theta, double* bi, double*bin, double* thnew, double* invcov,
  int* ncores, double* spmin, double* lwght, double* wght,
  double* swjy, int* np1, int* np2, int* np3);
void F77_NAME(sector)(double* x1, int* n1, double* x2, int* n2,
  int* nsect, int* sect, int* symm, double* insect);
void F77_NAME(segment)(double* y, int* fix, double* level,
  double* delta, double* si2, int* n1, int* n2, int* n3, double* hakt,
  double* lambda, double* theta, double* bi, double* bi2,
  double* bi0, double* gi, double* vred, double* thetan, int* kern,
  double* spmin, double* lwght, double* wght, int* segm, int* segmn,
  double* beta, double* thresh, double* ext, double* fov, double* varest);
void F77_NAME(vaws)(double* y, int* mask, int* nv, int* n1, int* n2,
  int* n3, double* hakt, double* lambda, double* theta, double* bi, double* bin,
  double* vred, double* thnew, int* ncores, double* spmin, double* lwght,
  double* wght, double* swjy);
  void F77_NAME(vaws2)(double* y, int* mask, int* nv, int* nvd, int* n1,
    int* n2, int* n3, double* hakt, double* lambda, double* theta, double* bi,
    double* vred, double* thnew, double* invcov, int* ncores, double* spmin,
    double* lwght, double* wght, double* swjy, double* thi, double* invcovi);
void F77_NAME(vpaws)(int* n, int* dp2, double* bi, double* bi2, double* var);
void F77_NAME(fillpat1)(double* x, int* n1, int* phw, int* psize, double* pmat);
void F77_NAME(fillpat2)(double* x, int* n1, int* n2, int* phw, int* psize,
  double* pmat);
void F77_NAME(fillpat3)(double* x, int* n1, int* n2, int* n3, int* phw,
  int* psize, double* pmat);
void F77_NAME(nlmeans)(double* x, int* n1, int* n2, int* n3, double* patch,
  int* pd, int* swd, double* tau, double* xhat);
static R_NativePrimitiveArgType awsph1_t[]={REALSXP, REALSXP, LGLSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType awsph2_t[]={REALSXP, REALSXP, LGLSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,REALSXP, REALSXP, REALSXP, REALSXP,
  INTSXP};
static R_NativePrimitiveArgType awsp1b_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType awsp2_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType caws_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType cawsmask_t[]={REALSXP, LGLSXP, INTSXP, LGLSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType cawsw_t[]={INTSXP, INTSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType cawsw1_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType caws1_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType caws6_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType cgaws_t[]={REALSXP, LGLSXP, LGLSXP, REALSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType cgawsmas_t[]={REALSXP, LGLSXP, INTSXP, LGLSXP,
  REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType chaws_t[]={REALSXP, LGLSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType chaws1_t[]={REALSXP, REALSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
  REALSXP};
static R_NativePrimitiveArgType exceed_t[]={REALSXP, INTSXP, REALSXP, INTSXP,
  REALSXP};
static R_NativePrimitiveArgType gethani_t[]={REALSXP, REALSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType getvofh_t[]={REALSXP, INTSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType lkern1_t[]={REALSXP, INTSXP, REALSXP, INTSXP,
  INTSXP, REALSXP};
static R_NativePrimitiveArgType mask_t[]={LGLSXP, LGLSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType median1d_t[]={REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType median2d_t[]={REALSXP, INTSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType median3d_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  REALSXP};
static R_NativePrimitiveArgType mpaws1_t[]={INTSXP, INTSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType mpaws2_t[]={INTSXP, INTSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType pawswght_t[]={INTSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, INTSXP, REALSXP};
static R_NativePrimitiveArgType pcaws_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType pcaws2_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType pcaws3_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType pcawsm_t[]={REALSXP, INTSXP, INTSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, INTSXP};
static R_NativePrimitiveArgType pvaws_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType pvaws2_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, INTSXP, INTSXP};
static R_NativePrimitiveArgType sector_t[]={REALSXP, INTSXP, REALSXP, INTSXP,
  INTSXP, INTSXP, LGLSXP, REALSXP};
static R_NativePrimitiveArgType segment_t[]={REALSXP, LGLSXP, REALSXP, REALSXP,
  REALSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, INTSXP,
  INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType vaws_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, INTSXP,
  REALSXP, REALSXP, REALSXP, REALSXP};
static R_NativePrimitiveArgType vaws2_t[]={REALSXP, LGLSXP, INTSXP, INTSXP,
  INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP, REALSXP, INTSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP,
  REALSXP};
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
            {"awsph1", (DL_FUNC) &awsph1_ , 22, awsph1_t},
            {"awsph2", (DL_FUNC) &awsph2_ ,23, awsph2_t},
            {"awsp1b", (DL_FUNC) &awsp1b_ ,21, awsp1b_t},
            {"awsp2", (DL_FUNC) &awsp2_ ,22, awsp2_t},
            {"caws", (DL_FUNC) &caws_ ,18, caws_t},
            {"cawsmask", (DL_FUNC) &cawsmask_ ,18, cawsmask_t},
            {"cawsw", (DL_FUNC) &cawsw_ ,12, cawsw_t},
            {"cawsw1", (DL_FUNC) &cawsw1_ ,16, cawsw1_t},
            {"caws1", (DL_FUNC) &caws1_ ,12, caws1_t},
            {"caws6", (DL_FUNC) &caws6_ ,18, caws6_t},
            {"cgaws", (DL_FUNC) &cgaws_ ,21, cgaws_t},
            {"cgawsmas", (DL_FUNC) &cgawsmas_ ,20, cgawsmas_t},
            {"chaws", (DL_FUNC) &chaws_ ,19, chaws_t},
            {"chaws1", (DL_FUNC) &chaws1_ ,14, chaws1_t},
            {"exceed", (DL_FUNC) &exceed_ ,5, exceed_t},
            {"gethani", (DL_FUNC) &gethani_ ,7, gethani_t},
            {"getvofh", (DL_FUNC) &getvofh_ ,4, getvofh_t},
            {"lkern1", (DL_FUNC) &lkern1_ ,6, lkern1_t},
            {"mask", (DL_FUNC) &mask_ ,5, mask_t},
            {"median1d", (DL_FUNC) &median1d_ ,3, median1d_t},
            {"median2d", (DL_FUNC) &median2d_ ,4, median2d_t},
            {"median3d", (DL_FUNC) &median3d_ ,5, median3d_t},
            {"mpaws1", (DL_FUNC) &mpaws1_ ,8, mpaws1_t},
            {"mpaws2", (DL_FUNC) &mpaws2_ ,8, mpaws2_t},
            {"pawswght", (DL_FUNC) &pawswght_ ,17, pawswght_t},
            {"pcaws", (DL_FUNC) &pcaws_ ,18, pcaws_t},
            {"pcaws2", (DL_FUNC) &pcaws2_ ,19, pcaws2_t},
            {"pcaws3", (DL_FUNC) &pcaws3_ ,19, pcaws3_t},
            {"pcawsm", (DL_FUNC) &pcawsm_ ,18, pcawsm_t},
            {"pvaws", (DL_FUNC) &pvaws_ ,20, pvaws_t},
            {"pvaws2", (DL_FUNC) &pvaws2_ ,22, pvaws2_t},
            {"sector", (DL_FUNC) &sector_ ,8, sector_t},
            {"segment", (DL_FUNC) &segment_ ,28, segment_t},
            {"vaws", (DL_FUNC) &vaws_ ,17, vaws_t},
            {"vaws2", (DL_FUNC) &vaws2_ ,21, vaws2_t},
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
