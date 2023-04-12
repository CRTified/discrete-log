#include<stdio.h>
#include<stdlib.h>
#include<stdarg.h>
#include<assert.h>
#include<inttypes.h>
#include<stdbool.h>
#include<pthread.h>
#include<gmp.h>
#include <endian.h>

#include <sys/stat.h>
#include <sys/types.h>

#include "ecc.h"

#define MAX_THREADS 7
#define NUM_FILEBITS 10
#define NUM_FILES (1 << NUM_FILEBITS)

typedef struct {
  EC_ProjPoint *P_expl;
  EC_ProjPoint *R_expl;
  uint64_t start;
  uint64_t end;

  size_t hash_size;
  size_t d_size;
  size_t entry_size;
  char *buf;
  FILE **files;

  uint64_t nullbits;
} thread_args;

int compare(const void *a_raw, const void *b_raw)
{
  uint8_t *a = (uint8_t *)a_raw;
  uint8_t *b = (uint8_t *)b_raw;

  for(size_t i = 0; ; i++)
  {
    if(a[i] < b[i]) return -1;
    if(a[i] > b[i]) return  1;
  }
}

void *thread(void *raw_args)
{
  thread_args *args = (thread_args *)raw_args;


  if(args->start >= args->end)
  {
    pthread_exit(NULL);
  }
  init_ecc(args->P_expl->curve->p);

  EC_ProjPoint res = initPoint(args->R_expl->curve);
  ec_copy(&res, args->R_expl);
  mult_ui(&res, args->start, args->R_expl);
  scale(&res);
  //printCurvePoint("res", &res);

  I_ProjPoint resI;
  EC_ProjPoint i_xP = initPoint(args->P_expl->curve);
  resI.coord_xP = &i_xP;
  EC_ProjPoint i_yP = initPoint(args->P_expl->curve);
  resI.coord_yP = &i_yP;
  EC_ProjPoint i_zP = initPoint(args->P_expl->curve);
  resI.coord_zP = &i_zP;

  uint8_t *buf = malloc(args->entry_size);
  uint64_t findex = 0;

  liftPoint(&resI, &res, args->P_expl);
  hashPoint((char*)buf, &resI);

  for(size_t j = 0; j < args->d_size; j++)
  {
    buf[args->hash_size + j] = (args->start >> ((args->d_size - 1 - j) * 8)) & 0xFF;
  }
  findex = be64toh(((uint64_t *)buf)[0]) >> (64 - args->nullbits - NUM_FILEBITS);
  fwrite(buf, args->entry_size, 1, args->files[findex]);

  if(args->start <= 1)
  {
    dbl(&res, &res);
  } else {
    add(&res, &res, args->R_expl);
  }

  for(uint64_t i = args->start + 1; i < args->end; i++)
  {
    //printf("%lu\n", i);
    scale(&res);
    //printCurvePoint("res", &res);
    liftPoint(&resI, &res, args->P_expl);
    //printCurvePoint("resI-x", resI.coord_xP);
    hashPoint((char *)buf, &resI);

    for(size_t j = 0; j < args->d_size; j++)
    {
      buf[args->hash_size + j] = (i >> ((args->d_size - j - 1) * 8)) & 0xFF;
    }
    findex = be64toh(((uint64_t *)buf)[0]) >> (64 - args->nullbits - NUM_FILEBITS);
    fwrite(buf, args->entry_size, 1, args->files[findex]);

    add(&res, &res, args->R_expl);
  }

  mpz_clears(
    res.coord_x,
    res.coord_y,
    res.coord_z,

    i_xP.coord_x,
    i_xP.coord_y,
    i_xP.coord_z,

    i_yP.coord_x,
    i_yP.coord_y,
    i_yP.coord_z,

    i_zP.coord_x,
    i_zP.coord_y,
    i_zP.coord_z,

    NULL
    );

  free(buf);

  clear_ecc();
  pthread_exit(NULL);
}

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

int main(int argc, char *argv[])
{
  if(argc != 1 + 13)
  {
    printf("Usage:\n");
    printf("  %s p E_a4 E_a6 Px Py q I_a4 I_a6 Rx Ry order_R factor name\n", argv[0]);
    return 1;
  }


  mkdir(argv[13], 0777);

  EC_Param explicit_curve = parseCurveParam(argv[1], argv[2], argv[3]);
  printCurveParam(&explicit_curve);

  init_ecc(explicit_curve.p);

  EC_ProjPoint P = parseAffineCurvePoint(&explicit_curve, argv[4], argv[5]);
  printCurvePoint("P", &P);

  EC_Param implicit_curve = parseCurveParam(argv[6], argv[7], argv[8]);
  printCurveParam(&implicit_curve);

  EC_ProjPoint R = parseAffineCurvePoint(&implicit_curve, argv[9], argv[10]);
  printCurvePoint("R", &R);

  mpz_t order_R;
  mpz_init_set_str(order_R, argv[11], 10);
  gmp_printf("Order: %Zd\n", order_R);

  uint64_t factor;
  if(1 != sscanf(argv[12], "%" SCNu64 "", &factor))
  {
    return 1;
  }
  printf("Factor: %" PRIu64 "\n", factor);

  mpz_t Ndivf;
  mpz_init(Ndivf);

  if(mpz_tdiv_q_ui(Ndivf, order_R, factor))
  {
    printf("Factor does not divide order\n");
    return 1;
  }
  gmp_printf("Ndivf: %Zd\n", Ndivf);

  EC_ProjPoint NdivfR = initPoint(&implicit_curve);

  mult(&NdivfR, Ndivf, &R);
  scale(&NdivfR);
  printCurvePoint("NdivfR", &NdivfR);

  size_t entries = (factor + 1) >> 1;
  size_t d_size = (64 - __builtin_clzll(entries)) / 8 + 1;
  size_t entry_size = explicit_curve.hash_bytes + d_size;
  printf("Entries: %lu\n", entries);
  printf("Size for coordinate: %lu\n", explicit_curve.hash_bytes);
  printf("Size for d: %lu\n", d_size );
  printf("Size per entry: %lu\n", entry_size);

  printf("Generating\n");

  FILE *files[NUM_FILES];


  char filepath[256];

  for(uint64_t i = 0; i < NUM_FILES; i++)
  {
    snprintf(filepath, 256, "./%s/%lu-%04" PRIx64 ".bin", argv[13], factor, i);
    files[i] = fopen(filepath, "wb+");
    if(files[i] == NULL){
      char errmsg[1024];
      perror(errmsg);
      printf("%s\n", errmsg);
      return 3;
    };
  }


  pthread_t generator_threads[MAX_THREADS] = { 0 };
  thread_args args[MAX_THREADS];
  size_t stepSize = 1 + (entries - 1) / MAX_THREADS;
  printf("Step Size: %lu\n", stepSize);
  for(int i = 0; i < MAX_THREADS; i++)
  {
    args[i].P_expl = &P;
    args[i].R_expl = &NdivfR;
    args[i].start  = MAX(i * stepSize, 1);
    args[i].end = MIN((i + 1) * stepSize, entries);
    args[i].hash_size = explicit_curve.hash_bytes;
    args[i].d_size = d_size;
    args[i].entry_size = entry_size;
    args[i].files = files;
    args[i].nullbits = (8 * explicit_curve.hash_bytes)
      - (mpz_sizeinbase(explicit_curve.p, 2) + 1);

    //printf("%2lu %2lu\n", args[i].start, args[i].end);
    pthread_create(generator_threads + i, NULL, thread, args + i);
  }

  for(int i = 0; i < MAX_THREADS; i++)
  {
    pthread_join(generator_threads[i], NULL);
  }

  printf("Writing\n");


  FILE *f;
  snprintf(filepath, 256, "./%s/%lu.bin", argv[13], factor);
  f = fopen(filepath, "wb");
  char *buf = calloc(entry_size, 1);
  fwrite(buf, 1, entry_size, f);

  long int fsize;
  for(uint64_t i = 0; i < NUM_FILES; i++)
  {
    fseek(files[i], 0, SEEK_END);
    fsize = ftell(files[i]);
    rewind(files[i]);

    buf = realloc(buf, fsize);

    if(fsize != fread(buf, 1, fsize, files[i]))
    {
      printf("RIP, %lu\n", fsize);
      return 3;
    }

    qsort(buf, fsize / entry_size, entry_size, compare);

    fwrite(buf, fsize / entry_size, entry_size, f);
    fclose(files[i]);

    snprintf(filepath, 256, "./%s/%lu-%04" PRIx64 ".bin", argv[13], factor, i);
    remove(filepath);
  }

  free(buf);
  fclose(f);

  mpz_clears(explicit_curve.A
	     , explicit_curve.B
	     , explicit_curve.p
	     , explicit_curve.halfp
	     , implicit_curve.A
	     , implicit_curve.B
	     , implicit_curve.p
	     , implicit_curve.halfp
	     , P.coord_x
	     , P.coord_y
	     , P.coord_z
	     , R.coord_x
	     , R.coord_y
	     , R.coord_z
	     , order_R
	     , Ndivf
	     , NdivfR.coord_x
	     , NdivfR.coord_y
	     , NdivfR.coord_z
	     , NULL
    );
  clear_ecc();
  return 0;

}
