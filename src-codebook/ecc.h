#include<stdlib.h>
#include<assert.h>
#include<inttypes.h>
#include<stdbool.h>
#include<stdio.h>
#include<string.h>
#include<threads.h>
#include<gmp.h>

typedef struct {
  mpz_t A;
  mpz_t B;
  mpz_t p;
  mpz_t halfp;
  size_t hash_bytes;
} EC_Param;

typedef struct {
  EC_Param *curve;
  mpz_t coord_x;
  mpz_t coord_y;
  mpz_t coord_z;
  bool on_curve;
} EC_ProjPoint;

typedef struct {
  EC_ProjPoint *coord_xP;
  EC_ProjPoint *coord_yP;
  EC_ProjPoint *coord_zP;
} I_ProjPoint;

void init_ecc(mpz_t p);
void clear_ecc();

EC_Param parseCurveParam(char *p, char *A, char *B);
void printCurveParam(EC_Param *param);

EC_ProjPoint initPoint(EC_Param *param);
void ec_copy(EC_ProjPoint *dst, const EC_ProjPoint *src);
EC_ProjPoint parseAffineCurvePoint(EC_Param *param, char *x_raw, char *y_raw);
void printCurvePoint(char *name, EC_ProjPoint *point);

void scale(EC_ProjPoint *Q);
void add(EC_ProjPoint *dst, const EC_ProjPoint *op1, const EC_ProjPoint *op2);
void dbl(EC_ProjPoint *dst, const EC_ProjPoint *op1);
void mult(EC_ProjPoint *dst, const mpz_t n, const EC_ProjPoint *P);
void mult_ui(EC_ProjPoint *dst, uint64_t n, const EC_ProjPoint *P);

void liftPoint(I_ProjPoint *I_Q, const EC_ProjPoint *E_Q, const EC_ProjPoint *E_P);
void hashPoint(char *dst, const I_ProjPoint *I_Q);
