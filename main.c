#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#define THRESH 0.001
#define NUM_COLORS 10
#define NUM_ROOTS 4 //for testing
#define PIXEL_LEN 12
#define COLOR_DEPTH 51

void mul_im(double* a_re, double* a_im, double* b_re, double* b_im, double* c_re, double* c_im)
{
    *c_re = *a_re * *b_re - *a_im * *b_im;
    *c_im = *a_re * *b_im + *a_im * *b_re;
}

void div_im(double* a_re, double* a_im, double* b_re, double* b_im, double* c_re, double* c_im)
{
    if((*b_re * *b_re + *b_im * *b_im) < THRESH)
    {
        *c_re = 1000;
        *c_im = 1000;
    }
    else
    {
        *c_re = (*a_re * *b_re + *a_im * *b_im) / (*b_re * *b_re + *b_im * *b_im);
        *c_im = (*b_re * *a_im - *a_re * *b_im) / (*b_re * *b_re + *b_im * *b_im);
    }
}

void derivative(int d, double* re, double* im, double* d_re, double* d_im, double * t_re, double * t_im)
{
    if(d == 3) //3 x^2
    {
        mul_im(re, im, re, im, d_re, d_im);
        *d_re = *d_re * 3;
        *d_im = *d_im * 3;
    }

    else if(d == 5) //5 x^4
    {
        mul_im(re, im, re, im, d_re, d_im);
        mul_im(d_re, d_im, re, im, t_re, t_im);
        mul_im(t_re, t_im, re, im, d_re, d_im);
        *d_re = *d_re * 5;
        *d_im = *d_im * 5;
    }
}

void fun_val(int d, double* re, double* im, double* f_re, double* f_im, double* t1_re, double* t1_im, double* t2_re, double* t2_im)
{
    // it has several arguments, because the values have to be stored somewhere on the way, could be improved
    if(d == 3) //x^3 - 1
    {
        *f_re = *re;
        *f_im = *im;
        mul_im(re, im, re, im, t1_re, t1_im);
        mul_im(t1_re, t1_im, f_re, f_im, t2_re, t2_im);
        *f_re = *t2_re - 1;
        *f_im = *t2_im;
    }

    if(d == 5) //x^5 - 1
    {
        *f_re = *re;
        *f_im = *im;
        mul_im(re, im, re, im, t1_re, t1_im);
        mul_im(t1_re, t1_im, t1_re, t1_im, t2_re, t2_im);
        mul_im(t2_re, t2_im, f_re, f_im, t1_re, t1_im);

        *f_re = *t1_re - 1;
        *f_im = *t1_im;
    }
}

double abs_val2(double* re, double* im, double* c_re, double* c_im)
{
    return (*re - *c_re)*(*re - *c_re) + (*im - *c_im)*(*im - *c_im);
}

double** precomputed_roots(int d)
{
    if(d == 3)
    {
        double ** as = (double**) malloc(sizeof(double*) * 4);
            for ( size_t ix = 0; ix < 4; ++ix )
                as[ix] = (double*) malloc(sizeof(double) * 2);
        as[0][0] = 0;
        as[0][1] = 0;
        as[1][0] = 1;
        as[1][1] = 0;
        as[2][0] = -0.5;
        as[2][1] = 0.866;
        as[3][0] = -0.5;
        as[3][1] = -0.866;
        return as;
    }
    if(d == 5)
    {
        double ** as = (double**) malloc(sizeof(double*) * 6);
            for ( size_t ix = 0; ix < 4; ++ix )
                as[ix] = (double*) malloc(sizeof(double) * 2);
        as[0][0] = 0;
        as[0][1] = 0;
        as[1][0] = 1;
        as[1][1] = 0;

        as[2][0] = -0.809;
        as[2][1] = 0.58779;
        as[3][0] = -0.809;
        as[3][1] = -0.58779;

        as[4][0] = 0.309;
        as[4][1] = 0.951;
        as[5][0] = 0.309;
        as[5][1] = -0.951;
        return as;
    }
    else
    {
        double ** as = (double**) malloc(sizeof(double*) * 3);
            for ( size_t ix = 0; ix < 3; ++ix )
                as[ix] = (double*) malloc(sizeof(double) * 2);
        return as;
    }
}

void* thread_placeholder(void)
{
    //
}

int main(int argc, char *argv[])
{
    int param_t, param_l;
    char *endptr;

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

    if(param_t > 1)
    {
        pthread_t thread_id[param_t];
        for(int i=0; i<param_t; ++i)
        {
            pthread_create(&thread_id[i], NULL, &thread_placeholder, NULL);
        }

        for(int i=0; i<param_t; ++i)
            pthread_join(thread_id[i],NULL);
    }

    //storing parameter values as string to write later
    char string_l[10];
    char string_t[5];
    char char_exponent[1];
    char header[30] = "P3\n";
    sprintf(string_l, "%ld", param_l);
    sprintf(string_t, "%ld", param_t);
    sprintf(char_exponent, "%ld", exponent);

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

    printf("t = %ld, l = %ld, exponent is %ld\n", param_t, param_l, exponent);
    printf("filename1: %s\n", filename_attractors);
    printf("filename2: %s\n", filename_convergence);

    double** roots_list = precomputed_roots(exponent);

    int ** as = (int**) malloc(sizeof(int*) * param_l); //array of roots
    for ( size_t ix = 0; ix < param_l; ++ix )
        as[ix] = (int*) malloc(sizeof(int) * param_l);

    int ** it = (int**) malloc(sizeof(int*) * param_l); //array of iterations
    for ( size_t ix = 0; ix < param_l; ++ix )
        it[ix] = (int*) malloc(sizeof(int) * param_l);

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
    for ( int ix=0; ix < param_l; ++ix ) {
        for ( int jx=0; jx < param_l; ++jx ){
            double re = 4.0*(ix - param_l/2.0)/(double)param_l;
            double im = 4.0*(jx - param_l/2.0)/(double)param_l;

            double d_re, d_im, f_re, f_im;
            double t1, t2, t3, t4;
            int converged = 0;
            int root = 0;

            while(!converged)
            {
                it[ix][jx]++;
                fun_val(exponent, &re, &im, &f_re, &f_im, &t1, &t2, &t3, &t4);
                derivative(exponent, &re, &im, &d_re, &d_im, &t3, &t4);
                div_im(&f_re, &f_im, &d_re, &d_im, &t1, &t2);
                re = re - t1;
                im = im - t2;
                if(re*re > 1000 || im*im > 1000){
                    converged = 1;
                    root = exponent;
                    break;
                }
                if(it[ix][jx] > 50){
                    converged = 1;
                    root = 0;
                    break;
                }
                for(int i = 0; i < NUM_ROOTS; i++)
                {
                    double abs = abs_val2(&re, &im, &roots_list[i][0], &roots_list[i][1]);
                    if(abs < THRESH*THRESH)
                    {
                        converged = 1;
                        root = i;
                        break;
                    }
                }
            }
            as[ix][jx] = root;
        }
    }

    FILE * pFile;
    pFile = fopen(filename_attractors, "w");
    printf("FIRST FILE:\n");
    //fprintf(pFile, "P3\n%ld %ld\n255\n", param_l, param_l);
    fwrite(&header, sizeof(char), strlen(header), pFile);

    for ( int ix=0; ix < param_l; ++ix ) {
        char pixels[15 * param_l + 1];
        char* p = pixels;
        for ( int jx=0; jx < param_l; ++jx ){
            char* color = char_lookup_table[as[ix][jx]];
            for(int j = 0; j<PIXEL_LEN; j++)
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
            char* grey = grey_lookup[v];
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

    free(as);
    free(it);
    fclose (pFile);
}
