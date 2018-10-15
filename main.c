#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <math.h>
#include <pthread.h>

#define THRESH 0.001
#define NUM_COLORS 10
#define NUM_ROOTS 4 //for testing
#define PIXEL_LEN 12
#define COLOR_DEPTH 51
#define M_PI 3.14159265358979323846

struct compute_runner_struct {
	short int ** as;
	short int ** it;
	double ** roots_list;
	int exponent;
	int param_l;
	int ix_start;
	int ix_stop;
	int ix_step;
};

void power_im(double * a_re, double * a_im, double * t_re, double * t_im, int  n)
{
    double temp;

    if (n == 0)
    {
        * a_re = 1;
        * a_im = 0;
    }
    else if (n % 2 == 0)
    {

        power_im(a_re, a_im, t_re, t_im, n / 2);
        temp = *a_re * *a_re - *a_im * *a_im;
        *a_im = *a_re * *a_im + *a_im * *a_re;
        *a_re = temp;
    }
    else
    {
        power_im(a_re, a_im, t_re, t_im, n-1);
        temp = *a_re * *t_re - *a_im * *t_im;
        *a_im = *a_re * *t_im + *a_im * *t_re;
        *a_re = temp;
    }
}





double abs_val2(double* re, double* im, double* c_re, double* c_im)
{
    return (*re - *c_re)*(*re - *c_re) + (*im - *c_im)*(*im - *c_im);
}

void precomputed_roots(int d, double ** roots_list)
{
    roots_list[0][0] = 0;
    roots_list[0][1] = 0;
    for(size_t ix=1; ix <= d; ++ix)
    {
        roots_list[ix][0] = cos(2* M_PI * ix / (double)d);
        roots_list[ix][1] = sin(2 * M_PI * ix / (double)d);
    }
}

void* compute_runner(void* arg)
{
    struct compute_runner_struct *arg_struct =
			(struct compute_runner_struct*) arg;

    double re, im;
    double d_re, d_im;
    double t1, t2, t3, t4;
    int converged;
    int root;
    double c_1 = 1 - 1 / (double)arg_struct->exponent;
    double c_2 = arg_struct->exponent;
    int c_3 = arg_struct->exponent - 1;
    for ( int ix=arg_struct->ix_start; ix <= arg_struct->ix_stop; ix += arg_struct->ix_step ) {
        for ( int jx=0; jx < arg_struct->param_l; ++jx ){
            re = 4.0*(ix - arg_struct->param_l/2.0)/(double)arg_struct->param_l;
            im = 4.0*(jx - arg_struct->param_l/2.0)/(double)arg_struct->param_l;

            converged = 0;
            root = 0;

            while(!converged)
            {
                arg_struct->it[ix][jx]++;
                d_re = c_1 * re;
                d_im = c_1 * im;
                t1 = re;
                t2 = im;

                power_im(&re, &im, &t1, &t2, c_3);
                t3 = (re * re + im * im);
                t4 = t3 * c_2;
                re =  re / t4 ;
                im = -im /t4;
                re = re + d_re;
                im = im + d_im;

                if( re*re > 10000000000 || im*im > 10000000000 || t3 < 0.00000000000001 ){
                    converged = 1;
                    root = 0;
                    break;
                }

                for(int i = 0; i <= c_2; i++)
                {
                    double abs = abs_val2(&re, &im, &arg_struct->roots_list[i][0], &arg_struct->roots_list[i][1]);
                    if(abs < THRESH*THRESH)
                    {
                        converged = 1;
                        root = i;
                        break;
                    }
                }
            }
            arg_struct->as[ix][jx] = root;
        }
    }
    return NULL;
}

int main(int argc, char *argv[])
{
    int param_t, param_l;

    if(argc != 4)
    {
        printf("The number of arguments is invalid. Expected: -t# -l# #\n");
        exit(1);
    }

    for(int i=0; i<argc; i++)
    {
        if (argv[i][0] == '-' && argv[1][1] != '\0') {
            if(argv[i][1] == 't') param_t = atoi(&argv[i][2]);
            if(argv[i][1] == 'l') param_l = atoi(&argv[i][2]);
        }
    }

    int exponent = atoi(argv[argc-1]);



    //storing parameter values as string to write later
    char string_l[10];
    char string_t[5];
    char char_exponent[1];
    char header[30] = "P3\n";
    sprintf(string_l, "%d", param_l);
    sprintf(string_t, "%d", param_t);
    sprintf(char_exponent, "%d", exponent);

    char filename_attractors[30] = "newton_attractors_x";
    char filename_convergence[30] = "newton_convergence_x";
    strcat(filename_attractors, char_exponent);
    strcat(filename_attractors, ".ppm");
    strcat(filename_convergence, char_exponent);
    strcat(filename_convergence, ".ppm");
    strcat(header, string_l);
    strcat(header, " ");
    strcat(header, string_l);
    strcat(header, "\n255\n");

    printf("t = %d, l = %d, exponent is %d\n", param_t, param_l, exponent);
    printf("filename1: %s\n", filename_attractors);
    printf("filename2: %s\n", filename_convergence);

    double ** roots_list = (double**) malloc(sizeof(double*) * (exponent+1)*2); //array of roots
    for ( size_t ix = 0; ix < exponent+1; ++ix )
        roots_list[ix] = (double*) malloc(sizeof(double) * 2);

    precomputed_roots(exponent, roots_list);

    short int ** as = (short int**) malloc(sizeof(short int*) * param_l); //array of roots
    short int * as_p = (short int*) malloc(sizeof(short int) * param_l*param_l);
    for ( size_t ix = 0, jx = 0; ix < param_l; ++ix, jx += param_l )
        as[ix] = as_p + jx;

    short int ** it = (short int**) malloc(sizeof(short int*) * param_l); //array of roots
    short int * it_p = (short int*) malloc(sizeof(short int) * param_l*param_l);
    for ( size_t ix = 0, jx = 0; ix < param_l; ++ix, jx += param_l )
        it[ix] = it_p + jx;

    //-----------------------------------------initialization:
    for ( int ix=0; ix < param_l; ++ix ) {
        for ( int jx=0; jx < param_l; jx+=2 ){
            it[ix][jx] = 0;
            it[ix][jx+1] = 0; //should it?
        }
    }
    //-----------------------------------------color lookup table:
    int ** color_table = (int**) malloc(sizeof(int*) * NUM_COLORS);
    for ( size_t ix = 0; ix < NUM_COLORS; ++ix )
        color_table[ix] = (int*) malloc(sizeof(int) * 3);

    color_table[0][0] = 130;
    color_table[0][1] = 20;
    color_table[0][2] = 130;

    color_table[1][0] = 135;
    color_table[1][1] = 250;
    color_table[1][2] = 250;

    color_table[2][0] = 135;
    color_table[2][1] = 150;
    color_table[2][2] = 250;

    color_table[3][0] = 200;
    color_table[3][1] = 200;
    color_table[3][2] = 250;

    color_table[4][0] = 130;
    color_table[4][1] = 0;
    color_table[4][2] = 130;

    color_table[5][0] = 100;
    color_table[5][1] = 0;
    color_table[5][2] = 50;

    color_table[6][0] = 180;
    color_table[6][1] = 0;
    color_table[6][2] = 90;

    color_table[7][0] = 255;
    color_table[7][1] = 165;
    color_table[7][2] = 0;

    color_table[8][0] = 220;
    color_table[8][1] = 225;
    color_table[8][2] = 0;

    color_table[9][0] = 100;
    color_table[9][1] = 255;
    color_table[9][2] = 150;

    //-----------------------------------------color lookup table:
    char ** char_lookup_table = (char**) malloc(sizeof(char*) * NUM_COLORS);
    for ( size_t ix = 0; ix < NUM_COLORS; ++ix )
        char_lookup_table[ix] = (char*) malloc(sizeof(char) * PIXEL_LEN);

    char_lookup_table[0] = "130 020 130 ";

    char_lookup_table[1] = "135 250 250 ";

    char_lookup_table[2] = "135 250 150 ";

    char_lookup_table[3] = "200 200 250 ";

    char_lookup_table[4] = "130 000 130 ";

    char_lookup_table[5] = "100 000 050 ";

    char_lookup_table[6] = "180 000 090 ";

    char_lookup_table[7] = "255 165 000 ";

    char_lookup_table[8] = "220 225 000 ";

    char_lookup_table[9] = "100 255 150 ";

    //-----------------------------------------silly grey lookup table:
    char ** grey_lookup = (char**) malloc(sizeof(char*) * COLOR_DEPTH);
    for ( size_t ix = 0; ix < COLOR_DEPTH; ++ix )
        grey_lookup[ix] = (char*) malloc(sizeof(char) * PIXEL_LEN);

    grey_lookup[0] = "000 000 000 ";
    grey_lookup[1] = "005 005 005 ";
    grey_lookup[2] = "010 010 010 ";
    grey_lookup[3] = "015 015 015 ";
    grey_lookup[4] = "020 020 020 ";
    grey_lookup[5] = "025 025 025 ";
    grey_lookup[6] = "030 030 030 ";
    grey_lookup[7] = "035 035 035 ";
    grey_lookup[8] = "040 040 040 ";
    grey_lookup[9] = "045 045 045 ";
    grey_lookup[10] = "050 050 050 ";
    grey_lookup[11] = "055 055 055 ";
    grey_lookup[12] = "060 060 060 ";
    grey_lookup[13] = "065 065 065 ";
    grey_lookup[14] = "070 070 070 ";
    grey_lookup[15] = "075 075 075 ";
    grey_lookup[16] = "080 080 080 ";
    grey_lookup[17] = "085 085 085 ";
    grey_lookup[18] = "090 090 090 ";
    grey_lookup[19] = "095 095 095 ";
    grey_lookup[20] = "100 100 100 ";
    grey_lookup[21] = "105 105 105 ";
    grey_lookup[22] = "110 110 110 ";
    grey_lookup[23] = "115 115 115 ";
    grey_lookup[24] = "120 120 120 ";
    grey_lookup[25] = "125 125 125 ";
    grey_lookup[26] = "130 130 130 ";
    grey_lookup[27] = "135 135 135 ";
    grey_lookup[28] = "140 140 140 ";
    grey_lookup[29] = "145 145 145 ";
    grey_lookup[30] = "150 150 150 ";
    grey_lookup[31] = "155 155 155 ";
    grey_lookup[32] = "160 160 160 ";
    grey_lookup[33] = "165 165 165 ";
    grey_lookup[34] = "170 170 170 ";
    grey_lookup[35] = "175 175 175 ";
    grey_lookup[36] = "180 180 180 ";
    grey_lookup[37] = "185 185 185 ";
    grey_lookup[38] = "195 190 190 ";
    grey_lookup[39] = "200 200 200 ";
    grey_lookup[40] = "210 210 210 ";
    grey_lookup[41] = "215 215 215 ";
    grey_lookup[42] = "225 225 225 ";
    grey_lookup[43] = "230 230 230 ";
    grey_lookup[44] = "235 235 235 ";
    grey_lookup[45] = "240 240 240 ";
    grey_lookup[46] = "240 240 240 ";
    grey_lookup[47] = "240 240 240 ";
    grey_lookup[48] = "245 245 245 ";
    grey_lookup[49] = "250 250 250 ";
    grey_lookup[50] = "250 250 250 ";
    grey_lookup[51] = "255 255 255 ";


    //-----------------------------------------COMPUTATION 1:

    struct compute_runner_struct args[param_t];
    pthread_t thread_id[param_t];
    for(size_t i=0; i<param_t; ++i)
    {
        args[i].as = as;
        args[i].it = it;
        args[i].roots_list = roots_list;
        args[i].exponent = exponent;
        args[i].param_l = param_l;
        args[i].ix_start = i;
        args[i].ix_step = param_t;
        if(param_l%param_t == 0)
            args[i].ix_stop = param_l - param_t + i;
        else
            args[i].ix_stop = param_l%param_t <= i ? param_l - param_l%param_t - param_t + i :  param_l - param_l%param_t + i;

        /* launch thread */
        pthread_attr_t attr;
		pthread_attr_init(&attr);
        pthread_create(&thread_id[i], &attr, compute_runner, &args[i]);

    }

    for(size_t i=0; i<param_t; ++i)
        pthread_join(thread_id[i],NULL);
    /* Writing */
    FILE * pFile;
    pFile = fopen(filename_attractors, "w");
    printf("FIRST FILE:\n");
    //fprintf(pFile, "P3\n%ld %ld\n255\n", param_l, param_l);
    fwrite(&header, sizeof(char), strlen(header), pFile);
    for ( size_t ix=0; ix < param_l; ++ix ) {

        char pixels[15 * param_l + 1];
        char* p = pixels;
        for ( size_t jx=0; jx < param_l; ++jx ){

            char* color = char_lookup_table[as[ix][jx]];

            for(size_t j = 0; j<PIXEL_LEN; j++)
            {
                *(p) = color[j];
                p++;
            }

        }
        *(p) = '\n';
        *(p+1) = '\0';
        fwrite(&pixels, sizeof(char), strlen(pixels), pFile);

    }

    fclose (pFile);
    printf("SECOND FILE:\n");
    pFile = fopen(filename_convergence, "w");

    //fprintf(pFile, "P3\n%ld %ld\n255\n", param_l, param_l);
    fwrite(&header, sizeof(char), strlen(header), pFile);

    for ( int ix=0; ix < param_l; ++ix ) {
        char pixels[15 * param_l + 1];
        char* p = pixels;
        for ( int jx=0; jx < param_l; ++jx ){
            int v = it[ix][jx]; //0-255
            char* grey = grey_lookup[v<50? v : 50];
            for(int i = 0; i<PIXEL_LEN; i++)
            {
                *(p) = grey[i];
                p++;
            }
        }
        *(p) = '\n';
        *(p+1) = '\0';
        fwrite(&pixels, sizeof(char), strlen(pixels), pFile);
    }
    free(roots_list);
    free(as);
    free(it);
    fclose (pFile);
}
