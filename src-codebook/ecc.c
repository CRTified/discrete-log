#include "ecc.h"

// Temporary variables used for computations. Preallocate with sufficient memory
thread_local mpz_t A, B, h, lhs, R, res, rhs, RR, s, ss, sss, u, uu, v, vv, vvv,
  w, X1Z2, XX, Y1Y1, Y1Z2, Y3, Z1Z2, ZZ;

void init_ecc(mpz_t p)
{
  size_t tgtsize = 32 * mpz_sizeinbase(p, 2);
  
  mpz_inits(A, B, h, lhs, R, res, rhs, RR, s, ss, sss, u, uu, v, vv, vvv, w, X1Z2, XX, Y1Y1, Y1Z2, Y3, Z1Z2, ZZ, NULL);
  mpz_realloc2(A, tgtsize);
  mpz_realloc2(B, tgtsize);
  mpz_realloc2(h, tgtsize);
  mpz_realloc2(lhs, tgtsize);
  mpz_realloc2(R, tgtsize);
  mpz_realloc2(res, tgtsize);
  mpz_realloc2(rhs, tgtsize);
  mpz_realloc2(RR, tgtsize);
  mpz_realloc2(s, tgtsize);
  mpz_realloc2(ss, tgtsize);
  mpz_realloc2(sss, tgtsize);
  mpz_realloc2(u, tgtsize);
  mpz_realloc2(uu, tgtsize);
  mpz_realloc2(v, tgtsize);
  mpz_realloc2(vv, tgtsize);
  mpz_realloc2(vvv, tgtsize);
  mpz_realloc2(w, tgtsize);
  mpz_realloc2(X1Z2, tgtsize);
  mpz_realloc2(XX, tgtsize);
  mpz_realloc2(Y1Y1, tgtsize);
  mpz_realloc2(Y1Z2, tgtsize);
  mpz_realloc2(Y3, tgtsize);
  mpz_realloc2(Z1Z2, tgtsize);
  mpz_realloc2(ZZ, tgtsize);

}

void clear_ecc()
{
  mpz_clears(A, B, h, lhs, R, res, rhs, RR, s, ss, sss, u, uu, v, vv, vvv, w, X1Z2, XX, Y1Y1, Y1Z2, Y3, Z1Z2, ZZ, NULL);
}

EC_Param parseCurveParam(char *p, char *A, char *B)
{
  EC_Param curve_param;
  mpz_init_set_str(curve_param.p, p, 10);
  assert(mpz_probab_prime_p(curve_param.p, 15));
  mpz_init(curve_param.halfp);
  mpz_tdiv_q_2exp(curve_param.halfp, curve_param.p, 1);
  
  mpz_init_set_str(curve_param.A, A, 10);
  mpz_mod(curve_param.A, curve_param.A, curve_param.p);

  mpz_init_set_str(curve_param.B, B, 10);
  mpz_mod(curve_param.B, curve_param.B, curve_param.p);

  curve_param.hash_bytes = mpz_sizeinbase(curve_param.p, 2) + 1;
  if(curve_param.hash_bytes & 0xF)
    // if it is not a multiple of 8, divide by 8 and add 1 for ceil
  {
    curve_param.hash_bytes >>= 3;
    curve_param.hash_bytes += 1;
  } else {
    curve_param.hash_bytes >>= 3;
  }
  
  return curve_param;
}

void printCurveParam(EC_Param *param)
{
  gmp_printf("Elliptic Curve y^2 = x^3 + Ax + B mod p with params\n"
	     "  p = %Zd\n"
	     "  A = %Zd\n"
	     "  B = %Zd\n",
	     param->p,
	     param->A,
	     param->B);
}

EC_ProjPoint initPoint(EC_Param *param)
{
  EC_ProjPoint point;
  point.curve = param;
  mpz_init_set_ui(point.coord_x, 0);
  mpz_init_set_ui(point.coord_y, 1);
  mpz_init_set_ui(point.coord_z, 0);
  point.on_curve = true;
  
  return point;
}

void ec_copy(EC_ProjPoint *dst, EC_ProjPoint const *src)
{
  dst->curve = src->curve;
  mpz_set(dst->coord_x, src->coord_x);
  mpz_set(dst->coord_y, src->coord_y);
  mpz_set(dst->coord_z, src->coord_z);
}

EC_ProjPoint parseAffineCurvePoint(EC_Param *param, char *x_raw, char *y_raw)
{
  EC_ProjPoint point;
  point.curve = param;
  
  mpz_init_set_str(point.coord_x, x_raw, 10);
  mpz_mod(point.coord_x, point.coord_x, point.curve->p);
  
  mpz_init_set_str(point.coord_y, y_raw, 10);
  mpz_mod(point.coord_y, point.coord_y, point.curve->p);
  
  mpz_init_set_ui(point.coord_z, 1);


  // y^2
  mpz_mul(lhs, point.coord_y, point.coord_y);
  mpz_mod(lhs, lhs, point.curve->p);

  // x^3
  mpz_mul(rhs, point.coord_x, point.coord_x);
  mpz_mul(rhs, rhs, point.coord_x);
  mpz_mod(rhs, rhs, point.curve->p);

  // x^3 + Ax
  mpz_addmul(rhs, point.curve->A, point.coord_x);
  mpz_mod(rhs, rhs, point.curve->p);

  // x^3 + Ax + B
  mpz_add(rhs, rhs, point.curve->B);

  mpz_mod(rhs, rhs, point.curve->p);
  
  point.on_curve = mpz_cmp(lhs, rhs) == 0 ? true : false;
  
  return point;
}

void printCurvePoint(char *name, EC_ProjPoint *point)
{
  gmp_printf("%s = (%Zd,\n"
	     "     %Zd,\n"
	     "     %Zd)\n"
	     "  %s\n",
	     name,
	     point->coord_x,
	     point->coord_y,
	     point->coord_z,
	     point->on_curve ? "is on curve" : "is not on curve");
}

void scale(EC_ProjPoint *Q)
{
  if(mpz_cmp_ui(Q->coord_z, 0) == 0)
  {
    mpz_set_ui(Q->coord_x, 0);
    mpz_set_ui(Q->coord_y, 1);
    mpz_set_ui(Q->coord_z, 0);
    return;
  }
  
  if(mpz_cmp_ui(Q->coord_z, 1) == 0)
  {
    return;
  }
  
  mpz_invert(A, Q->coord_z, Q->curve->p);
  mpz_mul(Q->coord_x, Q->coord_x, A);
  mpz_mod(Q->coord_x, Q->coord_x, Q->curve->p);

  mpz_mul(Q->coord_y, Q->coord_y, A);
  mpz_mod(Q->coord_y, Q->coord_y, Q->curve->p);

  mpz_set_ui(Q->coord_z, 1);
}


void add_mmadd_1998_cmo(EC_ProjPoint *dst, const EC_ProjPoint *op1, const EC_ProjPoint *op2)
{
  // u = Y2-Y1
  mpz_sub(u, op2->coord_y, op1->coord_y);
  
  // uu = u2
  mpz_mul(uu, u, u);
  //mpz_mod(uu, uu, p);

  // v = X2-X1
  mpz_sub(v, op2->coord_x, op1->coord_x);
  //mpz_mod(v, v, p);

  // vv = v2
  mpz_mul(vv, v, v);
  //mpz_mod(vv, vv, p);

  // vvv = v*vv
  mpz_mul(vvv, vv, v);
  //mpz_mod(vvv, vvv, p);

  // R = vv*X1
  mpz_mul(R, vv, op1->coord_x);
  //mpz_mod(R, R, p);

  // A = uu-vvv-2*R
  mpz_sub(A, uu, vvv);
  mpz_submul_ui(A, R, 2);
  //mpz_mod(A, A, p);

  // X3 = v*A
  mpz_mul(dst->coord_x, v, A);
  mpz_mod(dst->coord_x, dst->coord_x, op1->curve->p);
  
  // Y3 = u*(R-A)-vvv*Y1
  mpz_sub(Y3, R, A);
  mpz_mul(Y3, u, Y3);
  mpz_submul(Y3, vvv, op1->coord_y);
  mpz_mod(dst->coord_y, Y3, op1->curve->p);
  
  // Z3 = vvv*Z1
  mpz_mod(dst->coord_z, vvv, op1->curve->p);
}

void add_madd_1998_cmo(EC_ProjPoint *dst, const EC_ProjPoint *op1, const EC_ProjPoint *op2)
{
  // u = Y2*Z1-Y1
  mpz_mul(u, op2->coord_y, op1->coord_z);
  mpz_sub(u, u, op1->coord_y);
  //mpz_mod(u, u, p);
  
  // uu = u2
  mpz_mul(uu, u, u);
  //mpz_mod(uu, uu, p);

  // v = X2*Z1-X1
  
  mpz_mul(v, op2->coord_x, op1->coord_z);
  mpz_sub(v, v, op1->coord_x);
  //mpz_mod(v, v, p);

  // vv = v2
  mpz_mul(vv, v, v);
  //mpz_mod(vv, vv, p);

  // vvv = v*vv
  mpz_mul(vvv, vv, v);
  //mpz_mod(vvv, vvv, p);

  // R = vv*X1
  mpz_mul(R, vv, op1->coord_x);
  //mpz_mod(R, R, p);

  // A = uu*Z1-vvv-2*R
  mpz_mul(A, uu, op1->coord_z);
  mpz_sub(A, A, vvv);
  mpz_submul_ui(A, R, 2);
  //mpz_mod(A, A, p);

  // X3 = v*A
  mpz_mul(dst->coord_x, v, A);
  mpz_mod(dst->coord_x, dst->coord_x,op1->curve-> p);
  
  // Y3 = u*(R-A)-vvv*Y1
  mpz_sub(Y3, R, A);
  mpz_mul(Y3, u, Y3);
  mpz_submul(Y3, vvv, op1->coord_y);
  mpz_mod(dst->coord_y, Y3, op1->curve->p);
  
  // Z3 = vvv*Z1
  mpz_mul(dst->coord_z, vvv, op1->coord_z);
  mpz_mod(dst->coord_z, dst->coord_z, op1->curve->p);
}

void add_madd_1998_cmo2(EC_ProjPoint *dst, const EC_ProjPoint *op1, const EC_ProjPoint *op2)
{
  // Y1Z2 = Y1*Z2
  mpz_mul(Y1Z2, op1->coord_y, op2->coord_z);
  //mpz_mod(Y1Z2, Y1Z2, p);
  
  // X1Z2 = X1*Z2
  mpz_mul(X1Z2, op1->coord_x, op2->coord_z);
  //mpz_mod(X1Z2, X1Z2, p);
  
  // Z1Z2 = Z1*Z2
  mpz_mul(Z1Z2, op1->coord_z, op2->coord_z);
  //mpz_mod(Z1Z2, Z1Z2, p);
  
  // u = Y2*Z1-Y1Z2
  mpz_mul(u, op2->coord_y, op1->coord_z);
  mpz_sub(u, u, Y1Z2);
  //mpz_mod(u, u, p);
  
  // uu = u2
  mpz_mul(uu, u, u);
  //mpz_mod(uu, uu, p);
  
  // v = X2*Z1-X1Z2
  mpz_mul(v, op2->coord_x, op1->coord_z);
  mpz_sub(v, v, X1Z2);
  //mpz_mod(v, v, p);

  // vv = v2
  mpz_mul(vv, v, v);
  //mpz_mod(vv, vv, p);
  
  // vvv = v*vv
  mpz_mul(vvv, v, vv);
  //mpz_mod(vvv, vvv, p);
  
  // R = vv*X1Z2
  mpz_mul(R, vv, X1Z2);
  //mpz_mod(R, R, p);

  // A = uu*Z1Z2-vvv-2*R
  mpz_mul(A, uu, Z1Z2);
  mpz_sub(A, A, vvv);
  mpz_submul_ui(A, R, 2);
  //mpz_mod(A, A, p);

  // X3 = v*A
  mpz_mul(dst->coord_x, v, A);
  mpz_mod(dst->coord_x, dst->coord_x, op1->curve->p);
  
  // Y3 = u*(R-A)-vvv*Y1Z2
  mpz_sub(dst->coord_y, R, A);
  mpz_mul(dst->coord_y, u, dst->coord_y);
  mpz_submul(dst->coord_y, vvv, Y1Z2);
  mpz_mod(dst->coord_y, dst->coord_y,op1->curve-> p);
  
  // Z3 = vvv*Z1Z2
  mpz_mul(dst->coord_z, vvv, Z1Z2);
  mpz_mod(dst->coord_z, dst->coord_z, op1->curve->p);
}

void add(EC_ProjPoint *dst, const EC_ProjPoint *op1, const EC_ProjPoint *op2)
{
  if(mpz_cmp_ui(op1->coord_z, 0) == 0)
  {
    ec_copy(dst, op2);
    return;
  }
  
  if(mpz_cmp_ui(op2->coord_z, 0) == 0)
  {
    ec_copy(dst, op1);
    return;
  }
  
  
  if(mpz_cmp_ui(op2->coord_z, 1) == 0)
  {
    if(mpz_cmp_ui(op1->coord_z, 1) == 0)
    {
      add_mmadd_1998_cmo(dst, op1, op2);
    } else {
      add_madd_1998_cmo(dst, op1, op2);
    }
  } else {
    add_madd_1998_cmo2(dst, op1, op2);
  }
}

void mdbl_2007_bl(EC_ProjPoint *dst, const EC_ProjPoint *op1)
{
  // XX = X12
  mpz_mul(XX, op1->coord_x, op1->coord_x);
  //mpz_mod(XX, XX, p);
  // w = a+3*XX
  mpz_mul_ui(w, XX, 3);
  mpz_add(w, w, op1->curve->A);
  //mpz_mod(w, w, p);
  // Y1Y1 = Y12
  mpz_mul(Y1Y1, op1->coord_y, op1->coord_y);
  //mpz_mod(Y1Y1, Y1Y1, p);
  // R = 2*Y1Y1
  mpz_mul_ui(R, Y1Y1, 2);
  //mpz_mod(R, R, p);
  
  // sss = 4*Y1*R
  mpz_mul_ui(sss, op1->coord_y, 4);  
  mpz_mul(sss, sss, R);
  
  // RR = R2
  mpz_mul(RR, R, R);
  //mpz_mod(RR, RR, p);
  // B = (X1+R)2-XX-RR
  mpz_add(B, op1->coord_x, R);
  mpz_mul(B, B, B);
  mpz_sub(B, B, XX);
  mpz_sub(B, B, RR);
  //mpz_mod(B, B, p);
  // h = w2-2*B
  mpz_mul(h, w, w);
  mpz_submul_ui(h, B, 2);
  //mpz_mod(h, h, p);
  // X3 = 2*h*Y1
  mpz_mul(dst->coord_x, h, op1->coord_y);
  mpz_mul_ui(dst->coord_x, dst->coord_x, 2);
  mpz_mod(dst->coord_x, dst->coord_x, op1->curve->p);
  // Y3 = w*(B-h)-2*RR
  mpz_sub(dst->coord_y, B, h);
  mpz_mul(dst->coord_y, w, dst->coord_y);
  mpz_submul_ui(dst->coord_y, RR, 2);
  mpz_mod(dst->coord_y, dst->coord_y, op1->curve->p);
  // Z3 = sss
  mpz_mod(dst->coord_z, sss, op1->curve->p);

}

void dbl_2007_bl(EC_ProjPoint *dst, const EC_ProjPoint *op1)
{
  // XX = X12
  mpz_mul(XX, op1->coord_x, op1->coord_x);
  //mpz_mod(XX, XX, p);
  
  // ZZ = Z12
  mpz_mul(ZZ, op1->coord_z, op1->coord_z);
  //mpz_mod(ZZ, ZZ, p);
  
  // w = a*ZZ+3*XX
  mpz_mul_ui(w, XX, 3);
  mpz_addmul(w, op1->curve->A, ZZ);
  //mpz_mod(w, w, p);
  
  // s = 2*Y1*Z1
  mpz_mul(s, op1->coord_y, op1->coord_z);
  mpz_mul_ui(s, s, 2);
  //mpz_mod(s, s, p);
  
  // ss = s2
  mpz_mul(ss, s, s);
  //mpz_mod(ss, ss, p);

  // sss = s*ss
  mpz_mul(sss, s, ss);
  //mpz_mod(sss, sss, p);
  
  // R = Y1*s
  mpz_mul(R, op1->coord_y, s);
  //mpz_mod(R, R, p);
  
  // RR = R2
  mpz_mul(RR, R, R);
  //mpz_mod(RR, RR, p);
  
  // B = (X1+R)2-XX-RR
  mpz_add(B, op1->coord_x, R);
  mpz_mul(B, B, B);
  mpz_sub(B, B, XX);
  mpz_sub(B, B, RR);
  //mpz_mod(B, B, p);
  // h = w2-2*B
  mpz_mul(h, w, w);
  mpz_submul_ui(h, B, 2);
  //mpz_mod(h, h, p);
  
  // X3 = h*s
  mpz_mul(dst->coord_x, h, s);
  mpz_mod(dst->coord_x, dst->coord_x, op1->curve->p);
  
  // Y3 = w*(B-h)-2*RR
  mpz_sub(dst->coord_y, B, h);
  mpz_mul(dst->coord_y, w, dst->coord_y);
  mpz_submul_ui(dst->coord_y, RR, 2);
  mpz_mod(dst->coord_y, dst->coord_y, op1->curve->p);
  
  // Z3 = sss
  mpz_mod(dst->coord_z, sss, op1->curve->p);
}

void dbl(EC_ProjPoint *dst, const EC_ProjPoint *op1)
{
  if(mpz_cmp_ui(op1->coord_z, 1) == 0)
  {
    mdbl_2007_bl(dst, op1);
  } else {
    dbl_2007_bl(dst, op1);
  }
}

void mult(EC_ProjPoint *dst, const mpz_t n, const EC_ProjPoint *P)
{
  EC_ProjPoint tmp = initPoint(P->curve);
  mpz_set(tmp.coord_x, P->coord_x);
  mpz_set(tmp.coord_y, P->coord_y);
  mpz_set(tmp.coord_z, P->coord_z);

  mpz_set_ui(dst->coord_x, 0);
  mpz_set_ui(dst->coord_y, 1);
  mpz_set_ui(dst->coord_z, 0);
  
  size_t e_len = mpz_sizeinbase(n, 2);
  for(size_t i = 0; i < e_len; i++)
  {
    if(mpz_tstbit(n, i))
    {
      add(dst, dst, &tmp);
    }
    dbl(&tmp, &tmp);
  }
  mpz_clears(tmp.coord_x, tmp.coord_y, tmp.coord_z, NULL);
}
void mult_ui(EC_ProjPoint *dst, uint64_t n, const EC_ProjPoint *P)
{
  EC_ProjPoint tmp = initPoint(P->curve);
  mpz_set(tmp.coord_x, P->coord_x);
  mpz_set(tmp.coord_y, P->coord_y);
  mpz_set(tmp.coord_z, P->coord_z);

  mpz_set_ui(dst->coord_x, 0);
  mpz_set_ui(dst->coord_y, 1);
  mpz_set_ui(dst->coord_z, 0);

  while(n)
  {
    if(n & 1)
    {
      add(dst, dst, &tmp);
    }
    dbl(&tmp, &tmp);
    n >>= 1;
  }
  mpz_clears(tmp.coord_x, tmp.coord_y, tmp.coord_z, NULL);
}

void liftPoint(I_ProjPoint *I_Q, const EC_ProjPoint *E_Q, const EC_ProjPoint *E_P)
{
  mult(I_Q->coord_xP, E_Q->coord_x, E_P);
  scale(I_Q->coord_xP);
  // mult(I_Q->coord_yP, E_Q->coord_y, E_P);
  // scale(I_Q->coord_yP);
  // mult(I_Q->coord_zP, E_Q->coord_z, E_P);
  // scale(I_Q->coord_zP);
}

void hashPoint(char *dst, const I_ProjPoint *I_Q)
{
  EC_ProjPoint *tmp = I_Q->coord_xP;
  if(0 == mpz_cmp_ui(tmp->coord_z, 0))
  {
    memset(dst, 0, tmp->curve->hash_bytes);
    return;
  }
  mpz_mul_2exp(res, tmp->coord_x, 1);
  if(mpz_tstbit(tmp->coord_y, 0))
  {
    mpz_add_ui(res, res, 1);
  }
  mpz_export(dst, NULL, 1, tmp->curve->hash_bytes, 1, 0, res);
}
