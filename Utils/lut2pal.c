/* Adapted from http://nsweb.tn.tudelft.nl/~gsteele/spyview/ 
 * Converts a .lut file from the Lut directory of WsXM to be used as a palette in gnuplot. 
 * Instructions:
 *  Compile as lut2pal
 *  run lut2pal <palette>.lut
 *  use in gnuplot as: load '<palette>.lut.pal'
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define BUFSIZE 1024*5

void usage()
{
  fprintf(stderr, "lut2pal file.lut (outputs file.lut.pal file.lut.control_points)\n");
  fprintf(stderr, " to be loaded in gnuplot using \"load 'file.lut.pal'\".\n");
  exit(0);
}

int main(int argc, char **argv)
{
  int *red_xcp, *grn_xcp, *bl_xcp;
  int *red_ycp, *grn_ycp, *bl_ycp;
  int nrcp, ngcp, nbcp;
  char buf[BUFSIZE];
  int n, i;

  if (strcmp(argv[1], "-h") == 0 || argc < 1)
    usage();

  FILE *fp_pal, *fp_cp;
  snprintf(buf, 256, "%s.pal", argv[1]);
  fp_pal = fopen(buf, "w");
  snprintf(buf, 256, "%s.control_points", argv[1]);
  fp_cp = fopen(buf, "w");
    
  FILE *fp = fopen(argv[1], "r");
  if (fp == NULL) usage();
  n = fread(buf, 1, BUFSIZE, fp);
  fclose(fp);
  
  char *p;
  p = buf;
  p = strstr(p, "Control Points");
  sscanf(p, "Control Points: %d", &nbcp);
  p++;
  p = strstr(p, "Control Points");
  sscanf(p, "Control Points: %d", &ngcp);
  p++;
  p = strstr(p, "Control Points");
  sscanf(p, "Control Points: %d", &nrcp);
  
  bl_xcp = malloc(sizeof(int)*nbcp);
  bl_ycp = malloc(sizeof(int)*nbcp);
  grn_xcp = malloc(sizeof(int)*ngcp);
  grn_ycp = malloc(sizeof(int)*ngcp);
  red_xcp = malloc(sizeof(int)*nrcp);
  red_ycp = malloc(sizeof(int)*nrcp);

  int pntnum;
  p = strstr(buf, "[Blue Info]");
  for (n=0; n<nbcp; n++)
    {
      p = strstr(p, "Control Point");
      sscanf(p, "Control Point %d", &pntnum);
      if (sscanf(p, "Control Point %*d: (%d , %d)", &bl_xcp[pntnum], &bl_ycp[pntnum]) != 2)
	{fprintf(fp_cp, "error reading blue point %d\n", n); exit(-1);}
      bl_ycp[pntnum] = 255-bl_ycp[pntnum];
      p++;
    }

  p = strstr(buf, "[Green Info]");
  for (n=0; n<ngcp; n++)
    {
      p = strstr(p, "Control Point");
      sscanf(p, "Control Point %d", &pntnum);
      if (sscanf(p, "Control Point %*d: (%d , %d)", &grn_xcp[pntnum], &grn_ycp[pntnum]) != 2)
	{fprintf(fp_cp, "error reading green point %d\n", n); exit(-1);}
      grn_ycp[pntnum] = 255-grn_ycp[pntnum];
      p++;
    }

  p = strstr(buf, "[Red Info]");
  for (n=0; n<nrcp; n++)
    {
      p = strstr(p, "Control Point");
      sscanf(p, "Control Point %d", &pntnum);
      if (sscanf(p, "Control Point %*d: (%d , %d)", &red_xcp[pntnum], &red_ycp[pntnum]) != 2)
	{fprintf(fp_cp, "error reading red point %d\n", n); exit(-1);}
      red_ycp[pntnum] = 255-red_ycp[pntnum];
      p++;
    }

  int r[256];
  int g[256];
  int b[256];
  double slope;
  
  for (n=0; n<nrcp-1; n++)
    {
      slope = (double)(red_ycp[n+1]-red_ycp[n])/(double)(red_xcp[n+1]-red_xcp[n]);
      for (int i=red_xcp[n]; i<=red_xcp[n+1]; i++)
	r[i] = red_ycp[n] + slope*(double)(i-red_xcp[n]);
      fprintf(fp_cp, "%d %d \n", red_xcp[n],red_ycp[n]);
    }
  fprintf(fp_cp, "%d %d \n", red_xcp[n],red_ycp[n]);
  fprintf(fp_cp, "\n\n");

  for (n=0; n<ngcp-1; n++)
    {
      slope = (double)(grn_ycp[n+1]-grn_ycp[n])/(double)(grn_xcp[n+1]-grn_xcp[n]);
      for (int i=grn_xcp[n]; i<=grn_xcp[n+1]; i++)
	g[i] = grn_ycp[n] + slope*(double)(i-grn_xcp[n]);
      fprintf(fp_cp, "%d %d \n", grn_xcp[n],grn_ycp[n]);
    }
  fprintf(fp_cp, "%d %d \n", grn_xcp[n],grn_ycp[n]);
  fprintf(fp_cp, "\n\n");

  for (n=0; n<nbcp-1; n++)
    {
      slope = (double)(bl_ycp[n+1]-bl_ycp[n])/(double)(bl_xcp[n+1]-bl_xcp[n]);
      for (int i=bl_xcp[n]; i<=bl_xcp[n+1]; i++)
	b[i] = bl_ycp[n] + slope*(double)(i-bl_xcp[n]);
      fprintf(fp_cp, "%d %d\n", bl_xcp[n], bl_ycp[n]);
    }
  fprintf(fp_cp, "%d %d\n", bl_xcp[n], bl_ycp[n]);

  fprintf(fp_pal, "set palette defined(\\\n");
  for (int i=0; i<255; i++)
    fprintf(fp_pal, "%i %f %f %f,\\\n", i, r[i]/255., g[i]/255., b[i]/255.);
  i=255;
  fprintf(fp_pal, "%i %f %f %f)\n", i, r[i]/255., g[i]/255., b[i]/255.);

  return 0;
}
  
