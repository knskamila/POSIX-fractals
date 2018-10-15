# Assignment 2: Threads
### Kamila Kowalska, Sandra Viknander
---
## Introduction
The goal of the assignment is to provide a program utilizing POSIX threads to compute and represent convergence of all the points in a complex plane. The `newton` executable outputs two `.ppm` files and can be run with three parameters setting the number of threads, size of the images and the degree of the polynomial that is analyzed. The main difficulty is to fullfil the performance requirements, which makes synchronization between threads crucial.

## Design choices
While it would be possible to hardcode all the roots for all the polynomials considered, we decided against it. The computation of the roots is not time-consuming and performed only once at the start of the program, and such solution enables polynomials of degree 10 and higher to be analyzed. The function `precomputed_roots` handles computation of the right values given the exponent in the polynomial.

Similarly, there are no formulas for the value nor the derivative for a given exponent. Instead, the function `power_im` runs recursively to find the sought power of a given complex number.

## Parsing command line arguments
The condition is simply that there must be three arguments given. Arguments starting with respective letters are parsed using `atoi` function.

## Synchronization of compute and write threads
???

## Data transfer between compute and write threads
???

## Checking the convergence and divergence conditions
The points are assumed to have converged if their absolute distance to one of the roots is smaller than the threshold of 10^-3. Furthermore, the center of the coordinate space is treated as a root, since Newton method would never reach convergence for points close to it. Similarly, the iteration does not continue if any part of the number reaches 10^10. 
These conditions are checked in order every iteration, before doing the computations.
```
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
```
## Computation of of x_n in the iteration step.
As described above, the formulas aren't hardcoded and the program can calculate the iteration steps for any given exponent value. Each time step a corresponding value in an iteration matrix – storing the number of steps it took to converge – is incremented.
```
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

          (...)
}
```
## Writing to file
There are two separate loops which read values from `as` and `it` matrix, storing the root id and iteration number respectively, transforming them to an RGB color representation. All lines in the `.ppm` file are saved using `fwrite` function. To avoid slow performance, the use of string functions such as `strcat` was reduced to minimum.
Each line of pixels represented by R, G, B values is created by first allocating memory for a char array, then a pointer is used to iterate through memory addresses, so that an entire line can be saved at once:
```
char pixels[PIXEL_LEN * param_l + 1];
        char* p = pixels;
        for ( size_t jx=0; jx < param_l; ++jx ){

            char* color = char_lookup_table[as[ix][jx]];

            for(size_t j = 0; j<PIXEL_LEN; j++)
            {
                *(p) = color[j];
                p++;
            }
}
```
The `char_lookup_table` visible here is a predefined array of char values representing all the colors the roots can be visualized by. In the case of iteration number, the representation was limited to iterations in a [0, 51] range, where 51 iterations are represented by pure white. A corresponding lookup table is also predefined.

## Conclusions

While working on the assignment we could se the impact the use of the threads can make. Furthermore, the first version of the code utilized standard string functions, which made the code for a single thread and exponent 3 run for 17 seconds. Replacing these functions by pointer operations and using a lookup table of color values as char characters reduced this time to around 1 second.
