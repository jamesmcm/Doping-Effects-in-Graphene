#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>


#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

typedef struct header
{
   int randomisations;
   double vac_potential;
   double vac_concentration[2];
   char current;
   char gauge;
   int size[2];
   int wrap[2];
   int no_energy_vals;
   double flux;
   double seed;
}
Header;

#define MAX_FILE_NAME 60

void *xmalloc(size_t n);
Header *readfile(FILE *infile, int i, double *data[]);
int readoneline(char line[], int maxbytes, FILE *stream);
Header *readheader(FILE *infile);
void readdata(FILE *infile, double **data, int vals);
int checkdata(Header *head_a, double **data_a,
              Header *head_b, double **data_b);
void new_average(int no_energy_vals, double **data_a, int n_a,
                 double **data_b, int n_b, double **newdata);
double newmean(double mean_a, int n_a, double mean_b, int n_b);
double newerror(double mean_a, double err_a, int n_a,
                double mean_b, double err_b, int n_b);
void printnewfile(Header head, double **data, FILE *outfile);

int main(void)
{
   Header *inputhead[2], outputhead;
   double *inputdataone[4], *inputdatatwo[4], *outputdata[4];
   FILE *infile[2], *outfile;
   int i;
   int exitcode;

   printf("  Dataset merging program."
        "\n"
        "\n  Will merge two datasets provided they were generated with"
        "\n  the same parameters and have the same energy values."
        "\n\n");

   for (i = 0; i < 2; ++i) inputhead[i] = xmalloc(sizeof *inputhead[i]);

   inputhead[0] = readfile(infile[0], 0, inputdataone);
   inputhead[1] = readfile(infile[1], 1, inputdatatwo);

   exitcode=checkdata(inputhead[0], inputdataone, inputhead[1], inputdatatwo);
   if (exitcode != 0)
   {
      printf("  ERROR\n"
             "  Data files do not have matching parameters!\n Exit code: %i\n", exitcode);
      exit(0);
   }

   for (i = 0; i < 3; ++i)
   {
      outputdata[i] = xmalloc(inputhead[0]->no_energy_vals * 
                              sizeof(double));
   }

   new_average(inputhead[0]->no_energy_vals,
               inputdataone, inputhead[0]->randomisations,
               inputdatatwo, inputhead[1]->randomisations,
               outputdata);

   outputhead.randomisations       = inputhead[0]->randomisations + 
                                     inputhead[1]->randomisations;
   outputhead.vac_potential        = inputhead[0]->vac_potential;
   outputhead.vac_concentration[0] = inputhead[0]->vac_concentration[0];
   outputhead.vac_concentration[1] = inputhead[0]->vac_concentration[1];
   outputhead.current              = inputhead[0]->current;
   outputhead.gauge                = inputhead[0]->gauge;
   outputhead.size[0]              = inputhead[0]->size[0];
   outputhead.size[1]              = inputhead[0]->size[1];
   outputhead.wrap[0]              = inputhead[0]->wrap[0];
   outputhead.wrap[1]              = inputhead[0]->wrap[1];
   outputhead.no_energy_vals       = inputhead[0]->no_energy_vals;
   outputhead.flux                 = inputhead[0]->flux;
   outputhead.seed                 = inputhead[0]->seed;

   printnewfile(outputhead, outputdata, outfile);

   return 0;
}

void *xmalloc(size_t n)
{
   void *p = malloc(n);
   if (p == NULL)
   {
      fprintf(stderr, "Out of memory!\n");
      exit(1);
   }
	
   return p;
}

Header *readfile(FILE *infile, int i, double *data[])
{
   char name[MAX_FILE_NAME];
   int j, vals;

   Header *head = xmalloc(sizeof *head);

   printf("  Enter name of file %d: ", i+1);
   readoneline(name, MAX_FILE_NAME, stdin);
   printf("\n  Attempting to open \"%s\" ... ", name);

   if ((infile = fopen(name, "r")) == NULL)
   {
      printf("failed!\n");
      readfile(infile, i, data);
   }
   else
   {
      printf("success!\n");
      head = readheader(infile);
      vals = head->no_energy_vals;
      for (j = 0; j < 4; ++j)
      {
         data[j] = xmalloc(vals * sizeof(double));
      }
      readdata(infile, data, vals);
      return head;
   }
}

int readoneline(char line[], int maxbytes, FILE *stream)
{
   int length, j;
   while (1 == 1)
   {
      if (fgets(line, maxbytes, stream) == NULL)
      {
         fprintf(stderr, "Input NULL!\n");
         exit(1);
      }
		
      length = strlen(line) - 1;
      if (line[length] == '\n') line[length] = '\0';
		
      for (j = 0; line[j] != '\0'; ++j)
                               if (isspace(line[j]) == 0) return length;
   }
}

Header *readheader(FILE *infile)
{
   int initial, randomisations, size[2], wrap[2], no_energy_vals;
   char dummy[10], current, gauge;
   double vac_potential, vac_concentration[2], flux, seed;

   Header *head = xmalloc(sizeof *head); 

   initial = fgetc(infile);

   if (initial != '#')
   {
      printf("\n  Unexpected header, should start with \"#\"!"
             "\n  Exiting\n");
      exit(1);
   }

   fscanf(infile, "%d %lf %lf %lf %c %c %d %d %d %d %d %s %s %lf %lf",
          &randomisations, &vac_potential, &vac_concentration[0],
          &vac_concentration[1], &current, &gauge, &size[0], &size[1],
          &wrap[0], &wrap[1], &no_energy_vals, &dummy, &dummy, &flux,
          &seed);

   head->randomisations       = randomisations;
   head->vac_potential        = vac_potential;
   head->vac_concentration[0] = vac_concentration[0];
   head->vac_concentration[1] = vac_concentration[1];
   head->current              = current;
   head->gauge                = gauge;
   head->size[0]              = size[0];
   head->size[1]              = size[1];
   head->wrap[0]              = wrap[0];
   head->wrap[1]              = wrap[1];
   head->no_energy_vals       = no_energy_vals;
   head->flux                 = flux;
   head->seed                 = seed;

   return head;
}

void readdata(FILE *infile, double **data, int vals)
{
   int i, j;
   double dummy;

   for (j=0; j<vals; ++j)
   {
      for (i=0; i<4; ++i) fscanf(infile, "%lf", &data[i][j]);
      //fscanf(infile, "%lf", &dummy);
   }
   return;
}

int checkdata(Header *head_a, double **data_a,
              Header *head_b, double **data_b)
{
   int i, exitcode = 0;

   if (head_a->no_energy_vals == head_b->no_energy_vals)
   {
      for (i = 0; i < head_a->no_energy_vals; ++i)
      {
	if (data_a[0][i] != data_b[0][i]) {exitcode = 2;}
      }
   }
   else exitcode = 1;
   if (head_a->vac_potential != head_b->vac_potential ||
       head_a->vac_concentration[0] != head_b->vac_concentration[0] ||
       head_a->vac_concentration[1] != head_b->vac_concentration[1] ||
       head_a->current != head_b->current ||
       head_a->size[0] != head_b->size[0] ||
       head_a->size[1] != head_b->size[1] ||
       head_a->wrap[0] != head_b->wrap[0] ||
       head_a->wrap[1] != head_b->wrap[1] ||
       head_a->seed == head_b->seed ||
       head_a->flux != head_b->flux) exitcode = 3;
   if (head_a->flux != 0 && head_a->gauge != head_b->gauge) exitcode = 4;
   return exitcode;
}

void new_average(int no_energy_vals, double **data_a, int n_a,
                 double **data_b, int n_b, double **newdata)
{
   int i;

   for (i = 0; i < no_energy_vals; ++i)
   {
      newdata[0][i] = data_a[0][i];
      newdata[1][i] = newmean(data_a[1][i], n_a, data_b[1][i], n_b);
      newdata[2][i] = newerror(data_a[1][i], data_a[2][i], n_a,
                               data_b[1][i], data_b[2][i], n_b);
      newdata[3][i] = max(data_a[3][i], data_b[3][i]);

   }
}

double newmean(double mean_a, int n_a, double mean_b, int n_b)
{
   double mean;
   int n;

   mean_a *= n_a;
   mean_b *= n_b;
   n = n_a + n_b;
   mean = mean_a + mean_b;
   mean /= n;

   return mean;
}

double newerror(double mean_a, double err_a, int n_a,
                double mean_b, double err_b, int n_b)
{
   int n;
   double mean, err, sigma, sigma_a, sigma_b, d_fin, d_a, d_b;

   mean = newmean(mean_a, n_a, mean_b, n_b);
   n = n_a + n_b;

   sigma_a = err_a * sqrt(n_a);
   sigma_b = err_b * sqrt(n_b);

   d_a = ((sigma_a * sigma_a) + (mean_a * mean_a)) * n_a;
   d_b = ((sigma_b * sigma_b) + (mean_b * mean_b)) * n_b;

   d_fin = (d_a + d_b) / n;

   sigma = sqrt(d_fin - (mean * mean));
   err = sigma / sqrt(n);

   return err;
}

void printnewfile(Header head, double **data, FILE *outfile)
{
   char filename[MAX_FILE_NAME];
   int i, j;

   printf("  Please enter name for new file: ");
   readoneline(filename, MAX_FILE_NAME, stdin);

   printf("  Opening \"%s\" for writing... ", filename);
   if ((outfile = fopen(filename, "w")) == NULL)
   {
      printf("failed!\n");
      exit(1);
   }
   else printf("success!\n");

   fprintf(outfile,
           "# %d %g %g %g %c %c %d %d %d %d %d Custom E %g %g\n",
           head.randomisations, head.vac_potential,
           head.vac_concentration[0], head.vac_concentration[1],
           head.current, head.gauge, head.size[0], head.size[1],
           head.wrap[0], head.wrap[1], head.no_energy_vals, head.flux,
           head.seed);
   for (j = 0; j < head.no_energy_vals; ++j)
   {
      for (i = 0; i < 4; ++i) fprintf(outfile, "\t\t\t%g", data[i][j]);
      fprintf(outfile, "\n");
   }
   return;
}

