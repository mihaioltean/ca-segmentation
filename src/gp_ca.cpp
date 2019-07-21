//--------------------------------------------------------------------
//   Compiled with Microsoft Visual C++ 2013
//   Requires C++11 or newer (for thread support)
//--------------------------------------------------------------------
//   version 2016.08.29.1

#define USE_THREADS // comment this line if you don't have a >= C++11 compiler

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <limits.h>
#include <time.h>
#include <cstdlib>
#include <inttypes.h>

#include <random>

#ifdef USE_THREADS
#include <thread>
#include <mutex>
#endif

#include "qdbmp.h" // for bitmap read / write

#define CA_radius 3


#define num_variables (2 * CA_radius + 1) * (2 * CA_radius + 1) // starts from center then goes CCW from 0 degrees
#define num_operators 9

// +   min
// -   max
// *   if a then b else c
// /   if a < b then c else d

#define O_MIN -1
#define O_MAX -2
#define O_MEAN -3
#define O_DIFF -4
#define O_IFABCD -5
#define O_MAJORITAR_C -6
#define O_MINORITAR_C -7
#define O_MIN_C -8
#define O_MAX_C -9


char operators_string[9][20] = { "MIN", "MAX", "MEAN", "DIFF", "IFABCD", "MAJORITAR_COLOR", "MINORITAR_COLOR", "MIN_COMPONENT", "MAX_COMPONENT" };

#define color_tolerance 10
#define max_num_CA_iterations_with_no_improvements 10

//---------------------------------------------------------------------------
struct t_clip_region{
	unsigned long top_x, top_y;
	unsigned long width, height;
};
//---------------------------------------------------------------------------
struct t_rgb{
	unsigned char red, green, blue;
	t_rgb(void)
	{
		red = green = blue = 0;
	}
	t_rgb(unsigned char _r, unsigned char _g, unsigned char _b)
	{
		red = _r;
		green = _g;
		blue = _b;
	}
	void operator=(unsigned char c)
	{
		red = green = blue = c;
	}
	void replace_with_dominant(void)
	{
		if (red > green) {
			if (red > blue) {
				red = 255;
				green = 0;
				blue = 0;
			}
			else {
				red = 0;
				green = 0;
				blue = 255;
			}
		}
		else
			if (green > blue) {
				red = 0;
				green = 255;
				blue = 0;
			}
			else {
				red = 0;
				green = 0;
				blue = 255;
			}
	}
};
//---------------------------------------------------------------------------
t_rgb min_rgb(t_rgb &c1, t_rgb &c2)
{
	if (c1.red <= c2.red && c1.green <= c2.green && c1.blue <= c2.blue)
		return c1;
	else
		return c2;
}
//---------------------------------------------------------------------------
t_rgb max_rgb(t_rgb &c1, t_rgb &c2)
{
	if (c1.red >= c2.red && c1.green >= c2.green && c1.blue >= c2.blue)
		return c1;
	else
		return c2;
}
//---------------------------------------------------------------------------
t_rgb mean_rgb(t_rgb &c1, t_rgb &c2)
{
	return t_rgb((c1.red + c2.red) / 2, (c1.green + c2.green) / 2, (c1.blue + c2.blue) / 2);
}
//---------------------------------------------------------------------------
t_rgb diff_rgb(t_rgb &c1, t_rgb &c2)
{
	return t_rgb(abs(c1.red - c2.red), abs(c1.green - c2.green), abs(c1.blue - c2.blue));
}
//---------------------------------------------------------------------------
bool identical(t_rgb &a, t_rgb &b)
{
	return abs(a.red - b.red) <= color_tolerance && abs(a.green - b.green) <= color_tolerance && abs(a.blue - b.blue) <= color_tolerance;
}
//---------------------------------------------------------------------------
t_rgb majoritar_rgb(t_rgb &a, t_rgb &b, t_rgb &c)
{
	if (identical(a, b))
		return a;
	if (identical(a, c))
		return a;
	if (identical(b, c))
		return b;
	// one of them, random
	return a;
}
//---------------------------------------------------------------------------
t_rgb minoritar_rgb(t_rgb &a, t_rgb &b, t_rgb &c)
{
	if (identical(a, b))
		return c;
	if (identical(a, c))
		return b;
	if (identical(b, c))
		return a;
	// one of them, random
	return a;
}
//---------------------------------------------------------------------------
t_rgb ifabcd_rgb(t_rgb &a, t_rgb &b, t_rgb &c, t_rgb &d)
{
	if (a.red <= b.red && a.green <= b.green && a.blue <= b.blue)
		return c;
	else
		return d;
}
//---------------------------------------------------------------------------
t_rgb min_c(t_rgb &c)
{
	if (c.red <= c.green)
		if (c.red <= c.blue)
			return t_rgb(c.red, 0, 0);
		else
			return t_rgb(0, 0, c.blue);
	else
		if (c.green <= c.blue)
			return t_rgb(0, c.green, 0);
		else
			return t_rgb(0, 0, c.blue);
}
//---------------------------------------------------------------------------
t_rgb max_c(t_rgb &c)
{
	if (c.red >= c.green)
		if (c.red >= c.blue)
			return t_rgb(c.red, 0, 0);
		else
			return t_rgb(0, 0, c.blue);
	else
		if (c.green >= c.blue)
			return t_rgb(0, c.green, 0);
		else
			return t_rgb(0, 0, c.blue);
}
//---------------------------------------------------------------------------

struct t_code3{
	int op;				// either a variable, operator or constant; 
	// variables are indexed from 0: 0,1,2,...;
	// constants are indexed from num_variables
	// operators are -1, -2, -3...
	int adr1, adr2, adr3, adr4;    // pointers to arguments
};
//---------------------------------------------------------------------------
struct t_chromosome{
	t_code3 *prg;        // the program - a string of genes
	t_code3 *simplified_prg; // simplified program; unused instructions removed

	t_rgb *constants; // an array of constants

	int num_utilized_instructions; // num_utilized_instructions by the simplified program

	unsigned long fitness;        // the fitness (or the error)

	int num_steps_ca; // after how many steps the CA has stopped improving

	void mark(int k, bool* marked)
	{
		if ((prg[k].op < 0) && !marked[k]) {
			mark(prg[k].adr1, marked);

			switch (prg[k].op) {
			case O_MIN:
				mark(prg[k].adr2, marked);
				break;
			case O_MAX:
				mark(prg[k].adr2, marked);
				break;
			case O_MEAN:
				mark(prg[k].adr2, marked);
				break;
			case O_DIFF:
				mark(prg[k].adr2, marked);
				break;
			case O_MAJORITAR_C:
				mark(prg[k].adr2, marked);
				mark(prg[k].adr3, marked);
				break;
			case O_MINORITAR_C:
				mark(prg[k].adr2, marked);
				mark(prg[k].adr3, marked);
				break;
			case O_IFABCD:
				mark(prg[k].adr2, marked);
				mark(prg[k].adr3, marked);
				mark(prg[k].adr4, marked);
				break;
			}
		}
		marked[k] = true;
	}
	//---------------------------------------------------------------------------
	void simplify(int code_length)
	{
		bool *marked = new bool[code_length];
		for (int i = 0; i < code_length; marked[i++] = false);
		mark(code_length - 1, marked);

		// how many are skipped until a given instruction
		int *skipped = new int[code_length];
		if (!marked[0])
			skipped[0] = 1;
		else
			skipped[0] = 0;
		for (int i = 1; i < code_length; i++)
			if (!marked[i])
				skipped[i] = skipped[i - 1] + 1;
			else
				skipped[i] = skipped[i - 1];

		if (simplified_prg)
			delete[] simplified_prg;
		simplified_prg = new t_code3[code_length];

		num_utilized_instructions = 0;
		for (int i = 0; i < code_length; i++)
			if (marked[i]) {
				simplified_prg[num_utilized_instructions] = prg[i];
				if (prg[i].op < 0) {// operator
					simplified_prg[num_utilized_instructions].adr1 -= skipped[prg[i].adr1];
					simplified_prg[num_utilized_instructions].adr2 -= skipped[prg[i].adr2];
					simplified_prg[num_utilized_instructions].adr3 -= skipped[prg[i].adr3];
					simplified_prg[num_utilized_instructions].adr4 -= skipped[prg[i].adr4];
				}
				num_utilized_instructions++;
			}

		delete[] skipped;
		delete[] marked;
	}
	//------------------------------------------------------------------------------
	void to_string(char * s_dest, int code_length)
	{
		char tmp_s[100];
		s_dest[0] = 0;
		for (int i = 0; i < code_length; i++) {
			sprintf(tmp_s, "%d ", prg[i].op);
			strcat(s_dest, tmp_s);
			sprintf(tmp_s, "%d ", prg[i].adr1);
			strcat(s_dest, tmp_s);
			sprintf(tmp_s, "%d ", prg[i].adr2);
			strcat(s_dest, tmp_s);
			sprintf(tmp_s, "%d ", prg[i].adr3);
			strcat(s_dest, tmp_s);
			sprintf(tmp_s, "%d ", prg[i].adr4);
			strcat(s_dest, tmp_s);
		}

		sprintf(tmp_s, "%lg\n", fitness);
		strcat(s_dest, tmp_s);
	}
	//------------------------------------------------------------------------------
	void from_string(char* s_source, int code_length)
	{
		int num_consumed = 0;
		for (int i = 0; i < code_length; i++) {
			sscanf(s_source, "%d%d%d%n", &prg[i].op, &prg[i].adr1, &prg[i].adr2, &prg[i].adr3, &prg[i].adr4, &num_consumed);
			s_source += num_consumed;
		}
		sscanf(s_source, "%lf", &fitness);
	}
	//------------------------------------------------------------------------------


};
//---------------------------------------------------------------------------
struct t_parameters{
	int code_length;             // number of instructions in a chromosome
	int num_generations;
	int num_sub_populations;       // number of subpopulations
	int sub_population_size;                // subpopulation size
	double mutation_probability, crossover_probability;
	int num_constants;
	unsigned char constants_min, constants_max;   // the array for constants
	double variables_probability, operators_probability, constants_probability;

//	int num_CA_iterations;

#ifdef USE_THREADS
	int num_threads; // num threads. 
	//for best performances the number of subpopulations should be multiple of num_threads.
	// num_thread should no exceed the number of processor cores.
#endif
};
//---------------------------------------------------------------------------
struct t_seed{
	uint32_t z1, z2, z3, z4;
	t_seed(void)
	{
		z1 = z2 = z3 = z4 = 12345;
	}
	void init(uint32_t initial_seed, uint32_t seed)
	{
		z1 = z2 = z3 = z4 = initial_seed + seed;
	}
};
//---------------------------------------------------------------------------
uint32_t RNG(t_seed &seed)
{
	uint32_t b;
	b = ((seed.z1 << 6) ^ seed.z1) >> 13;
	seed.z1 = ((seed.z1 & 4294967294U) << 18) ^ b;
	b = ((seed.z2 << 2) ^ seed.z2) >> 27;
	seed.z2 = ((seed.z2 & 4294967288U) << 2) ^ b;
	b = ((seed.z3 << 13) ^ seed.z3) >> 21;
	seed.z3 = ((seed.z3 & 4294967280U) << 7) ^ b;
	b = ((seed.z4 << 3) ^ seed.z4) >> 12;
	seed.z4 = ((seed.z4 & 4294967168U) << 13) ^ b;
	return (seed.z1 ^ seed.z2 ^ seed.z3 ^ seed.z4);
}
//---------------------------------------------------------------------------
uint32_t my_int_rand(t_seed &seed, int _min, int _max)
{
	return RNG(seed) % (_max - _min + 1) + _min;
}
//---------------------------------------------------------------------------
double my_real_rand(t_seed &seed, double _min, double _max)
{
	return RNG(seed) * 2.3283064365386963e-10 * (_max - _min) + _min;
}
//---------------------------------------------------------------------------
void allocate_chromosome(t_chromosome &c, t_parameters &params)
{
	c.prg = new t_code3[params.code_length];
	c.simplified_prg = NULL;
	if (params.num_constants)
		c.constants = new t_rgb[params.num_constants];
	else
		c.constants = NULL;
}
//---------------------------------------------------------------------------
void delete_chromosome(t_chromosome &c)
{
	if (c.prg) {
		delete[] c.prg;
		c.prg = NULL;
	}
	if (c.simplified_prg) {
		delete[] c.simplified_prg;
		c.simplified_prg = NULL;
	}
	if (c.constants) {
		delete[] c.constants;
		c.constants = NULL;
	}
}
//---------------------------------------------------------------------------
void allocate_matrix(t_rgb **&work_matrix, unsigned long image_width, unsigned long image_height)
{
	work_matrix = new t_rgb*[image_height];
	for (unsigned long i = 0; i < image_height; i++)
		work_matrix[i] = new t_rgb[image_width];
}
//---------------------------------------------------------------------------
void allocate_matrices(t_rgb ***&matrices, unsigned long image_width, unsigned long image_height, int num_matrices)
{
	// for each thread we have a separate matrix
	matrices = new t_rgb**[num_matrices];
	for (int t = 0; t < num_matrices; t++)
		allocate_matrix(matrices[t], image_width, image_height);
}
//---------------------------------------------------------------------------
void delete_matrix(t_rgb **&work_matrix, unsigned long image_height)
{
	if (work_matrix) {
		for (unsigned long i = 0; i < image_height; i++)
			delete[] work_matrix[i];
		delete[] work_matrix;
	}
}
//---------------------------------------------------------------------------
void delete_matrices(t_rgb ***&matrices, unsigned long image_height, int num_matrices)
{
	if (matrices) {
		for (int t = 0; t < num_matrices; t++)
			delete_matrix(matrices[t], image_height);
		delete[] matrices;
	}
}
//---------------------------------------------------------------------------
void copy_individual(t_chromosome& dest, const t_chromosome& source, t_parameters &params)
{
	for (int i = 0; i < params.code_length; i++)
		dest.prg[i] = source.prg[i];
	for (int i = 0; i < params.num_constants; i++)
		dest.constants[i] = source.constants[i];
	dest.fitness = source.fitness;
	dest.num_steps_ca = source.num_steps_ca;
}
//---------------------------------------------------------------------------
void generate_random_chromosome(t_chromosome &a_chromosome, t_parameters &params, t_seed &seed) // randomly initializes the individuals
{
	// generate constants first
	for (int c = 0; c < params.num_constants; c++)
		a_chromosome.constants[c] = my_int_rand(seed, params.constants_min, params.constants_max);

	// on the first position we can have only a variable or a constant
	double sum = params.variables_probability + params.constants_probability;
	double p = my_real_rand(seed, 0, sum);

	if (p <= params.variables_probability)
		a_chromosome.prg[0].op = my_int_rand(seed, 0, num_variables - 1);
	else
		a_chromosome.prg[0].op = num_variables + my_int_rand(seed, 0, params.num_constants - 1);

	// for all other genes we put either an operator, variable or constant
	for (int i = 1; i < params.code_length; i++) {
		double p = my_real_rand(seed, 0, 1);

		if (p <= params.operators_probability)
			a_chromosome.prg[i].op = -1 - my_int_rand(seed, 0, num_operators - 1);        // an operator
		else
			if (p <= params.operators_probability + params.variables_probability)
				a_chromosome.prg[i].op = my_int_rand(seed, 0, num_variables - 1);     // a variable
			else
				a_chromosome.prg[i].op = num_variables + my_int_rand(seed, 0, params.num_constants - 1); // index of a constant

		a_chromosome.prg[i].adr1 = my_int_rand(seed, 0, i - 1);
		a_chromosome.prg[i].adr2 = my_int_rand(seed, 0, i - 1);
		a_chromosome.prg[i].adr3 = my_int_rand(seed, 0, i - 1);
		a_chromosome.prg[i].adr4 = my_int_rand(seed, 0, i - 1);
	}
}
//---------------------------------------------------------------------------
bool matrix_to_bmp_file(t_rgb **matrix, unsigned long image_width, unsigned long image_height, const char* file_name)
{
	BMP*    bmp;

	bmp = BMP_Create(image_width, image_height, 24);

	for (unsigned long x = 0; x < image_width; x++)
		for (unsigned long y = 0; y < image_height; y++)
			/* Invert RGB values */
			BMP_SetPixelRGB(bmp, x, y, matrix[y + CA_radius][x + CA_radius].red, matrix[y + CA_radius][x + CA_radius].green, matrix[y + CA_radius][x + CA_radius].blue);

	/* Save result */
	BMP_WriteFile(bmp, file_name);
	if (BMP_GetError() != BMP_OK) {
		BMP_Free(bmp);
		return false;
	}

	/* Free all memory allocated for the image */
	BMP_Free(bmp);
	return true;
}
//---------------------------------------------------------------------------
t_rgb evaluate_chromosome(t_code3 *prg, int code_length, t_rgb* constants, t_rgb **work_matrix, unsigned long x, unsigned long y, t_rgb *eval_array, t_rgb *variables)
{
	
	int count_vars = 0; 
	
	for (int i = -CA_radius; i <= CA_radius; i++)
		for (int j = -CA_radius; j <= CA_radius; j++)
			variables[count_vars++] = work_matrix[y + i][x + j];
	
	t_rgb *p = eval_array;
	for (int i = 0; i < code_length; i++)
		switch (prg[i].op) {
			case O_MIN:
				*p++ = min_rgb(eval_array[prg[i].adr1], eval_array[prg[i].adr2]);
				break;
			case O_MAX:
				*p++ = max_rgb(eval_array[prg[i].adr1], eval_array[prg[i].adr2]);
				break;
			case O_MEAN:
				*p++ = mean_rgb(eval_array[prg[i].adr1], eval_array[prg[i].adr2]);
				break;
			case O_DIFF:
				*p++ = diff_rgb(eval_array[prg[i].adr1], eval_array[prg[i].adr2]);
				break;
			case O_MAJORITAR_C:
				*p++ = majoritar_rgb(eval_array[prg[i].adr1], eval_array[prg[i].adr2], eval_array[prg[i].adr3]);
				break;
			case O_MINORITAR_C:
				*p++ = minoritar_rgb(eval_array[prg[i].adr1], eval_array[prg[i].adr2], eval_array[prg[i].adr3]);
				break;
			case O_IFABCD:
				*p++ = ifabcd_rgb(eval_array[prg[i].adr1], eval_array[prg[i].adr2], eval_array[prg[i].adr3], eval_array[prg[i].adr4]);
				break;
			case O_MIN_C:
				*p++ = min_c(eval_array[prg[i].adr1]);
				break;
			case O_MAX_C:
				*p++ = max_c(eval_array[prg[i].adr1]);
				break;
			default:
				if (prg[i].op < num_variables)
					*p++ = variables[prg[i].op];
				else
					*p++ = constants[prg[i].op - num_variables];
		}

	return eval_array[code_length - 1]; // last gene has the result	
}
//---------------------------------------------------------------------------
unsigned long matrix_quality(t_rgb **work_matrix, t_rgb **mask_matrix, t_clip_region & clip)
{
	long mean_background_red = 0;
	long mean_background_green = 0;
	long mean_background_blue = 0;

	long mean_object_red = 0;
	long mean_object_green = 0;
	long mean_object_blue = 0;

	int count_object_pixels = 0;
	int count_background_pixels = 0;

	for (unsigned long y = 0; y < clip.height; y++)
		for (unsigned long x = 0; x < clip.width; x++)
			if (!mask_matrix[y + clip.top_y][x + clip.top_x].red) {// black - background
				t_rgb tmp_rgb = work_matrix[y + CA_radius][x + CA_radius];
				mean_background_red += tmp_rgb.red;
				mean_background_green += tmp_rgb.green;
				mean_background_blue += tmp_rgb.blue;
				count_background_pixels++;
			}
			else {// white, object
				t_rgb tmp_rgb = work_matrix[y + CA_radius][x + CA_radius];
				mean_object_red += tmp_rgb.red;
				mean_object_green += tmp_rgb.green;
				mean_object_blue += tmp_rgb.blue;
				count_object_pixels++;
			}

			mean_background_red /= count_background_pixels;
			mean_background_green /= count_background_pixels;
			mean_background_blue /= count_background_pixels;

			mean_object_red /= count_object_pixels;
			mean_object_green /= count_object_pixels;
			mean_object_blue /= count_object_pixels;

			unsigned long variance_background = 0;
			unsigned long variance_object = 0;

			for (unsigned long y = 0; y < clip.height; y++)
				for (unsigned long x = 0; x < clip.width; x++)
					if (!mask_matrix[y + clip.top_y][x + clip.top_x].red)// black - background
						variance_background += abs((int)work_matrix[y + CA_radius][x + CA_radius].red - (int)mean_background_red) + abs((int)work_matrix[y + CA_radius][x + CA_radius].green - (int)mean_background_green) + abs((int)work_matrix[y + CA_radius][x + CA_radius].blue - (int)mean_background_blue);
					else// white, object
						variance_object += abs((int)work_matrix[y + CA_radius][x + CA_radius].red - (int)mean_object_red) + abs((int)work_matrix[y + CA_radius][x + CA_radius].green - (int)mean_object_green) + abs((int)work_matrix[y + CA_radius][x + CA_radius].blue - (int)mean_object_blue);;

			long extra_fitness = 0;
			if (abs((int)mean_background_red - (int)mean_object_red) + abs((int)mean_background_green - (int)mean_object_green) + abs((int)mean_background_blue - (int)mean_object_blue) < 50)
				extra_fitness = clip.width * clip.height * 256 * 3;

	return variance_background + variance_object + extra_fitness;
}
//---------------------------------------------------------------------------
void extend_matrix(t_rgb **work_matrix1, t_clip_region &clip)
{
	for (int r = 1; r <= CA_radius; r++) {
		for (unsigned long y = 0; y < clip.height + 2 * (r - 1); y++) {
			// first column
			work_matrix1[y + CA_radius - r + 1][CA_radius - r] = work_matrix1[y + CA_radius - r + 1][CA_radius - r + 1];
			// last column
			work_matrix1[y + CA_radius - r + 1][CA_radius + clip.width - 1 + r] = work_matrix1[y + CA_radius - r + 1][CA_radius + clip.width - 1 + r - 1];
		}
		// top, bottom row
		for (unsigned long x = 0; x < clip.width + 2 * (r - 1); x++) {
			// first row
			work_matrix1[CA_radius - r][x + CA_radius - r + 1] = work_matrix1[CA_radius - r + 1][x + CA_radius - r + 1];
			// last row
			work_matrix1[CA_radius + clip.height + r - 1][x + CA_radius - r + 1] = work_matrix1[CA_radius + clip.height - 1 + r - 1][x + CA_radius - r + 1];
		}
		// corners
		// top left
		work_matrix1[CA_radius - r][CA_radius - r] = work_matrix1[CA_radius - r + 1][CA_radius - r + 1];
		// bottom right
		work_matrix1[CA_radius + clip.height + r - 1][CA_radius + clip.width + r - 1] = work_matrix1[CA_radius + clip.height - 1 + r - 1][CA_radius + clip.width - 1 + r - 1];
		// top right
		work_matrix1[CA_radius - r][CA_radius + clip.width + r - 1] = work_matrix1[CA_radius - r + 1][CA_radius + clip.width - 1 + r - 1];
		// bottom left
		work_matrix1[CA_radius + clip.height + r - 1][CA_radius - r] = work_matrix1[CA_radius + clip.height - 1 + r - 1][CA_radius - r + 1];
	}
}
//---------------------------------------------------------------------------
unsigned long compute_fitness_1_image(t_chromosome &c, t_rgb **original_matrix, t_rgb **mask_matrix, t_clip_region &clip, t_rgb **work_matrix1, t_rgb **work_matrix2, int &num_steps_ca)
{
	// init
	for (unsigned long y = 0; y < clip.height; y++) {
		for (unsigned long x = 0; x < clip.width; x++)
			work_matrix1[y + CA_radius][x + CA_radius] = original_matrix[y + clip.top_y][x + clip.top_x];
	}
	extend_matrix(work_matrix1, clip);

//	matrix_to_bmp_file(work_matrix1, clip.width, clip.height, "c:\\temp\\test_ca.bmp");

	t_rgb *eval_array = new t_rgb[c.num_utilized_instructions];

	unsigned long fitness_1_image = ULONG_MAX;
	
	t_rgb variables[num_variables];
	int t = 0;
	int num_iterations_without_improvements = 0;

	num_steps_ca = 0;

	while (1) {
		if (t % 2) {
			for (unsigned long y = 0; y < clip.height; y++) {
				for (unsigned long x = 0; x < clip.width; x++)
					work_matrix1[y + CA_radius][x + CA_radius] = evaluate_chromosome(c.simplified_prg, c.num_utilized_instructions, c.constants, work_matrix2, x + CA_radius, y + CA_radius, eval_array, variables);
			}

			extend_matrix(work_matrix1, clip);
			unsigned long tmp_fitness = matrix_quality(work_matrix1, mask_matrix, clip);
			if (tmp_fitness < fitness_1_image) {
				fitness_1_image = tmp_fitness;
				num_iterations_without_improvements = 0;
				num_steps_ca = t;
			}
			else
				num_iterations_without_improvements++;
		}
		else {
			for (unsigned long y = 0; y < clip.height; y++) {
				for (unsigned long x = 0; x < clip.width; x++)
					work_matrix2[y + CA_radius][x + CA_radius] = evaluate_chromosome(c.simplified_prg, c.num_utilized_instructions, c.constants, work_matrix1, x + CA_radius, y + CA_radius, eval_array, variables);
			}

			extend_matrix(work_matrix2, clip);

			unsigned long tmp_fitness = matrix_quality(work_matrix2, mask_matrix, clip);
			if (tmp_fitness < fitness_1_image) {
				fitness_1_image = tmp_fitness;
				num_iterations_without_improvements = 0;
				num_steps_ca = t;
			}
			else
				num_iterations_without_improvements++;
		}

		t++;
		if (num_iterations_without_improvements >= max_num_CA_iterations_with_no_improvements)
			break; // exit while
	}

	/*
	if (num_CA_iterations % 2) {
		fitness_1_image = matrix_quality(work_matrix2, mask_matrix, clip);
//	if (image_file_name)
//	matrix_to_bmp_file(work_matrix1, clip.width, clip.height, image_file_name);
	}
	else {
		fitness_1_image = matrix_quality(work_matrix1, mask_matrix, clip);
//	if (image_file_name)
//	matrix_to_bmp_file(work_matrix2, clip.width, clip.height, image_file_name);
	}
	*/


	delete[] eval_array;

	return fitness_1_image;
}
//---------------------------------------------------------------------------
void compute_fitness(t_chromosome &c, int code_length, t_rgb ***original_matrices, t_rgb ***mask_matrices, int num_training_images, t_clip_region &clip, t_rgb **work_matrix1, t_rgb **work_matrix2)
{
	c.simplify(code_length);

	c.fitness = 0;
	c.num_steps_ca = 0;
	int CA_num_steps_1_image;
	for (int i = 0; i < num_training_images; i++) {
		c.fitness += compute_fitness_1_image(c, original_matrices[i], mask_matrices[i], clip, work_matrix1, work_matrix2, CA_num_steps_1_image);
		if (c.num_steps_ca < CA_num_steps_1_image)
			c.num_steps_ca = CA_num_steps_1_image;
	}
}
//---------------------------------------------------------------------------
void test_ca(t_chromosome &c, int code_length, t_rgb **original_matrix, unsigned long image_width, unsigned long image_height, char* image_file_name)
{
	c.simplify(code_length);

	t_clip_region full_image_clip;
	full_image_clip.top_x = 0;
	full_image_clip.top_y = 0;
	full_image_clip.width = image_width;
	full_image_clip.height = image_height;

	t_rgb **work_matrix1, **work_matrix2;
	allocate_matrix(work_matrix1, image_width + 2 * CA_radius, image_height + 2 * CA_radius);
	allocate_matrix(work_matrix2, image_width + 2 * CA_radius, image_height + 2 * CA_radius);

	t_rgb *eval_array = new t_rgb[c.num_utilized_instructions];
	t_rgb variables[num_variables];

	// init
	for (unsigned long y = 0; y < image_height; y++) {
		for (unsigned long x = 0; x < image_width; x++)
			work_matrix1[y + CA_radius][x + CA_radius] = original_matrix[y][x];
	}
	extend_matrix(work_matrix1, full_image_clip);

	for (int t = 0; t <= c.num_steps_ca; t++) {
		if (t % 2) {
			for (unsigned long y = 0; y < image_height; y++) {
				for (unsigned long x = 0; x < image_width; x++)
					work_matrix1[y + CA_radius][x + CA_radius] = evaluate_chromosome(c.simplified_prg, c.num_utilized_instructions, c.constants, work_matrix2, x + CA_radius, y + CA_radius, eval_array, variables);
			}
			extend_matrix(work_matrix1, full_image_clip);
		}
		else {
			for (unsigned long y = 0; y < image_height; y++)
				for (unsigned long x = 0; x < image_width; x++) {
					work_matrix2[y + CA_radius][x + CA_radius] = evaluate_chromosome(c.simplified_prg, c.num_utilized_instructions, c.constants, work_matrix1, x + CA_radius, y + CA_radius, eval_array, variables);
			}
			extend_matrix(work_matrix2, full_image_clip);
		}
	}
	if (c.num_steps_ca % 2) {
		if (image_file_name)
			matrix_to_bmp_file(work_matrix1, image_width, image_height, image_file_name);
	}
	else {
		if (image_file_name)
			matrix_to_bmp_file(work_matrix2, image_width, image_height, image_file_name);
	}

	delete_matrix(work_matrix1, image_height);
	delete_matrix(work_matrix2, image_height);

	delete[] eval_array;
}
//---------------------------------------------------------------------------
void mutation(t_chromosome &a_chromosome, t_parameters& params, t_seed& seed) // mutate the individual
{
	// mutate each symbol with the given probability
	// first gene must be a variable or constant
	double p = my_real_rand(seed, 0, 1);
	if (p < params.mutation_probability) {
		double sum = params.variables_probability + params.constants_probability;
		double p = my_real_rand(seed, 0, sum);

		if (p <= params.variables_probability)
			a_chromosome.prg[0].op = my_int_rand(seed, 0, num_variables - 1);
		else
			a_chromosome.prg[0].op = num_variables + my_int_rand(seed, 0, params.num_constants - 1);
	}
	// other genes
	for (int i = 1; i < params.code_length; i++) {
		p = my_real_rand(seed, 0, 1);      // mutate the operator
		if (p < params.mutation_probability) {
			// we mutate it, but we have to decide what we put here
			p = my_real_rand(seed, 0, 1);

			if (p <= params.operators_probability)
				a_chromosome.prg[i].op = -1 - my_int_rand(seed, 0, num_operators - 1);
			else
				if (p <= params.operators_probability + params.variables_probability)
					a_chromosome.prg[i].op = my_int_rand(seed, 0, num_variables - 1);
				else
					a_chromosome.prg[i].op = num_variables + my_int_rand(seed, 0, params.num_constants - 1); // index of a constant
		}

		p = my_real_rand(seed, 0, 1);      // mutate the first address  (adr1)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].adr1 = my_int_rand(seed, 0, i - 1);

		p = my_real_rand(seed, 0, 1);      // mutate the second address   (adr2)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].adr2 = my_int_rand(seed, 0, i - 1);

		p = my_real_rand(seed, 0, 1);      // mutate the second address   (adr2)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].adr3 = my_int_rand(seed, 0, i - 1);

		p = my_real_rand(seed, 0, 1);      // mutate the second address   (adr2)
		if (p < params.mutation_probability)
			a_chromosome.prg[i].adr4 = my_int_rand(seed, 0, i - 1);
	}
	// mutate the constants
	for (int c = 0; c < params.num_constants; c++) {
		p = my_real_rand(seed, 0, 1);
		if (p < params.mutation_probability)
			a_chromosome.constants[c] = my_int_rand(seed, params.constants_min, params.constants_max);
	}

	//a_chromosome.simplify(params.code_length);
}
//---------------------------------------------------------------------------
void one_cut_point_crossover(const t_chromosome &parent1, const t_chromosome &parent2, t_parameters &params, t_seed& seed, t_chromosome &offspring1, t_chromosome &offspring2)
{
	int cutting_pct = my_int_rand(seed, 0, params.code_length - 1);
	for (int i = 0; i < cutting_pct; i++) {
		offspring1.prg[i] = parent1.prg[i];
		offspring2.prg[i] = parent2.prg[i];
	}
	for (int i = cutting_pct; i < params.code_length; i++) {
		offspring1.prg[i] = parent2.prg[i];
		offspring2.prg[i] = parent1.prg[i];
	}
	// now the constants
	if (params.num_constants) {
		cutting_pct = my_int_rand(seed, 0, params.num_constants - 1);
		for (int i = 0; i < cutting_pct; i++) {
			offspring1.constants[i] = parent1.constants[i];
			offspring2.constants[i] = parent2.constants[i];
		}
		for (int i = cutting_pct; i < params.num_constants; i++) {
			offspring1.constants[i] = parent2.constants[i];
			offspring2.constants[i] = parent1.constants[i];
		}
	}
	//offspring1.simplify(params.code_length);
	//offspring2.simplify(params.code_length);
}
//---------------------------------------------------------------------------
void uniform_crossover(const t_chromosome &parent1, const t_chromosome &parent2, t_parameters &params, t_seed& seed, t_chromosome &offspring1, t_chromosome &offspring2)
{
	for (int i = 0; i < params.code_length; i++)
		if (my_int_rand(seed, 0, 1)) {
			offspring1.prg[i] = parent1.prg[i];
			offspring2.prg[i] = parent2.prg[i];
		}
		else {
			offspring1.prg[i] = parent2.prg[i];
			offspring2.prg[i] = parent1.prg[i];
		}

		// constants
		for (int i = 0; i < params.num_constants; i++)
			if (my_int_rand(seed, 0, 1)) {
				offspring1.constants[i] = parent1.constants[i];
				offspring2.constants[i] = parent2.constants[i];
			}
			else {
				offspring1.constants[i] = parent2.constants[i];
				offspring2.constants[i] = parent1.constants[i];
			}
			//offspring1.simplify(params.code_length);
			//offspring2.simplify(params.code_length);
}
//---------------------------------------------------------------------------
int sort_function(const void *a, const void *b)
{// comparator for quick sort
	if (((t_chromosome *)a)->fitness > ((t_chromosome *)b)->fitness)
		return 1;
	else
		if (((t_chromosome *)a)->fitness < ((t_chromosome *)b)->fitness)
			return -1;
		else
			return 0;
}
//---------------------------------------------------------------------------
bool print_chromosome(t_chromosome& a, t_parameters &params, char *filename)
{
	a.simplify(params.code_length);
	if (filename) {
		FILE *f = fopen(filename, "w");

		if (!f)
			return false;

		for (int i = 0; i < params.num_constants; i++)
			fprintf(f, "constants[%d] = (%d, %d, %d)\n", i, a.constants[i].red, a.constants[i].green, a.constants[i].blue);

		for (int i = 0; i < a.num_utilized_instructions; i++)
			if (a.simplified_prg[i].op < 0)
				fprintf(f, "%d: %s %d %d %d %d\n", i, operators_string[abs(a.simplified_prg[i].op) - 1], a.simplified_prg[i].adr1, a.simplified_prg[i].adr2, a.simplified_prg[i].adr3, a.simplified_prg[i].adr4);
			else
				if (a.simplified_prg[i].op < num_variables)
					fprintf(f, "%d: inputs[%d]\n", i, a.simplified_prg[i].op);
				else
					fprintf(f, "%d: constants[%d]\n", i, a.simplified_prg[i].op - num_variables);

		fprintf(f, "CA num interations = %ld\n", a.num_steps_ca);
		fprintf(f, "Fitness = %ld\n", a.fitness);

		fclose(f);

		return true;
	}
	else
		return false;
}
//---------------------------------------------------------------------------
int tournament_selection(t_chromosome *a_sub_pop, t_seed& seed, int sub_pop_size, int tournament_size)     // Size is the size of the tournament
{
	int r, p;
	p = my_int_rand(seed, 0, sub_pop_size - 1);
	for (int i = 1; i < tournament_size; i++) {
		r = my_int_rand(seed, 0, sub_pop_size - 1);
		p = a_sub_pop[r].fitness < a_sub_pop[p].fitness ? r : p;
	}
	return p;
}
//---------------------------------------------------------------------------
#ifdef USE_THREADS
void evolve_one_subpopulation(int *current_subpop_index, t_seed* seeds, std::mutex* mutex, t_chromosome ** sub_populations, int generation_index, t_parameters *params, t_rgb ***original_matrices, t_rgb ***mask_matrices, int num_training_images, t_clip_region *clip, t_rgb **work_matrix1, t_rgb **work_matrix2)
#else
void evolve_one_subpopulation(int *current_subpop_index, t_seed* seeds, t_chromosome ** sub_populations, int generation_index, t_parameters *params, t_rgb ***original_matrices, t_rgb ***mask_matrices, int num_training_images, t_clip_region *clip, t_rgb **work_matrix1, t_rgb **work_matrix2)

#endif
{
	int pop_index = 0;
	while (*current_subpop_index < params->num_sub_populations) {// still more subpopulations to evolve?
#ifdef USE_THREADS
		while (!mutex->try_lock()) {}// create a lock so that multiple threads will not evolve the same sub population
		pop_index = *current_subpop_index;
		(*current_subpop_index)++;
		mutex->unlock();
#else
		pop_index = *current_subpop_index;
		(*current_subpop_index)++;
#endif

		//srand(pop_index);
		// pop_index is the index of the subpopulation evolved by the current thread
		if (pop_index < params->num_sub_populations) {
			t_chromosome *a_sub_population = sub_populations[pop_index];

			t_chromosome offspring1, offspring2;
			allocate_chromosome(offspring1, *params);
			allocate_chromosome(offspring2, *params);

			if (generation_index == 0) {
				for (int i = 0; i < params->sub_population_size; i++) {
					//printf("pop index = %d  i = %d: ", pop_index, i);
					//printf("//////////////////////////////////\n");
					compute_fitness(a_sub_population[i], params->code_length, original_matrices, mask_matrices, num_training_images, *clip, work_matrix1, work_matrix2);
				}
				// sort ascendingly by fitness inside this population
				qsort((void *)a_sub_population, params->sub_population_size, sizeof(a_sub_population[0]), sort_function);
			}
			else // next generations
				for (int k = 0; k < params->sub_population_size; k += 2) {
					// we increase by 2 because at each step we create 2 offspring

					// choose the parents using binary tournament
					int r1 = tournament_selection(a_sub_population, seeds[pop_index], params->sub_population_size, 2);
					int r2 = tournament_selection(a_sub_population, seeds[pop_index], params->sub_population_size, 2);
					// crossover
					double p_0_1 = my_real_rand(seeds[pop_index], 0, 1);
					if (p_0_1 < params->crossover_probability)
						one_cut_point_crossover(a_sub_population[r1], a_sub_population[r2], *params, seeds[pop_index], offspring1, offspring2);
					else {// no crossover so the offspring are a copy of the parents
						copy_individual(offspring1, a_sub_population[r1], *params);
						copy_individual(offspring2, a_sub_population[r2], *params);
					}
					// mutate the result and compute fitness
					mutation(offspring1, *params, seeds[pop_index]);

					compute_fitness(offspring1, params->code_length, original_matrices, mask_matrices, num_training_images, *clip, work_matrix1, work_matrix2);

					// mutate the other offspring too
					mutation(offspring2, *params, seeds[pop_index]);

					compute_fitness(offspring2, params->code_length, original_matrices, mask_matrices, num_training_images, *clip, work_matrix1, work_matrix2);


					// replace the worst in the population
					if (offspring1.fitness < a_sub_population[params->sub_population_size - 1].fitness) {
						copy_individual(a_sub_population[params->sub_population_size - 1], offspring1, *params);
						qsort((void *)a_sub_population, params->sub_population_size, sizeof(a_sub_population[0]), sort_function);
					}
					if (offspring2.fitness < a_sub_population[params->sub_population_size - 1].fitness) {
						copy_individual(a_sub_population[params->sub_population_size - 1], offspring2, *params);
						qsort((void *)a_sub_population, params->sub_population_size, sizeof(a_sub_population[0]), sort_function);
					}
				}

			delete_chromosome(offspring1);
			delete_chromosome(offspring2);
		}
	}
}
//---------------------------------------------------------------------------
bool save_state(int num_sub_populations, int sub_population_size, int code_length, t_seed* seeds, t_chromosome** sub_populations, int current_generation_index, int initial_seed)
{
	char filename[100];
	sprintf(filename, "state_%d_%d.txt", current_generation_index, initial_seed);
	FILE* f = fopen(filename, "w");

	if (!f)
		return false;

	// save parameters
	fprintf(f, "%d %d %d %d\n", num_sub_populations, sub_population_size, code_length, current_generation_index);

	// save the population
	char* buffer = new char[code_length * 5 * 11];
	for (int p = 0; p < num_sub_populations; p++) {
		//save seeds
		fprintf(f, "%d %d %d %d\n", seeds[p].z1, seeds[p].z2, seeds[p].z3, seeds[p].z4);
		for (int c = 0; c < sub_population_size; c++) {
			sub_populations[p][c].to_string(buffer, code_length);
			fprintf(f, "%d\n", buffer);
		}
	}

	fclose(f);
}
//--------------------------------------------------------------------
bool load_state(const char* filename, t_seed* seeds, t_chromosome** sub_populations, int& current_generation_index)
{
	FILE* f = fopen(filename, "r");

	if (!f)
		return false;

	int num_sub_populations, sub_population_size, code_length;

	// load parameters
	fscanf(f, "%d %d %d %d\n", &num_sub_populations, &sub_population_size, &code_length, &current_generation_index);


	// load the population
	char* buffer = new char[code_length * 5 * 11];
	for (int p = 0; p < num_sub_populations; p++) {
		// load seeds
		fscanf(f, "%d %d %d %d\n", &seeds[p].z1, &seeds[p].z2, &seeds[p].z3, &seeds[p].z4);

		for (int c = 0; c < sub_population_size; c++) {
			fgets(buffer, code_length * 5 * 11, f);
			sub_populations[p][c].from_string(buffer, code_length);
		}
	}

	fclose(f);
}
//--------------------------------------------------------------------

void start_steady_state_gp(t_chromosome **sub_populations, t_parameters &params, uint32_t initial_seed, t_seed* seeds, t_rgb ***original_matrices, t_rgb ***mask_matrices, int num_images, t_clip_region &clip, int image_width, int image_height)
{
	// a steady state model -
	// Newly created inviduals replace the worst ones (if the offspring are better) in the same (sub) population.

#ifdef USE_THREADS
	// allocate memory for
	t_rgb *** work_matrices1, *** work_matrices2;
	allocate_matrices(work_matrices1, clip.width + 2 * CA_radius, clip.height + 2 * CA_radius, params.num_threads);
	allocate_matrices(work_matrices2, clip.width + 2 * CA_radius, clip.height + 2 * CA_radius, params.num_threads);

	// an array of threads. Each sub population is evolved by a thread
	std::thread **mep_threads = new std::thread*[params.num_threads];
	// we create a fixed number of threads and each thread will take and evolve one subpopulation, then it will take another one
	std::mutex mutex;
	// we need a mutex to make sure that the same subpopulation will not be evolved twice by different threads
#else
	t_rgb ** work_matrix1, ** work_matrix2;
	allocate_matrix(work_matrix1, clip.width + 2 * CA_radius, clip.height + 2 * CA_radius);
	allocate_matrix(work_matrix2, clip.width + 2 * CA_radius, clip.height + 2 * CA_radius);
#endif

	int num_training_images = 3;

	unsigned long best_fitness_so_far = LONG_MAX;
	int best_subpopulation_index = 0;
	// evolve for a fixed number of generations
	for (int generation = 0; generation < params.num_generations; generation++) { // for each generation

#ifdef USE_THREADS
		int current_subpop_index = 0;
		for (int t = 0; t < params.num_threads; t++)
			mep_threads[t] = new std::thread(evolve_one_subpopulation, &current_subpop_index, seeds, &mutex, sub_populations, generation, &params, original_matrices, mask_matrices, num_training_images, &clip, work_matrices1[t], work_matrices2[t]);

		for (int t = 0; t < params.num_threads; t++) {
			mep_threads[t]->join();
			delete mep_threads[t];
		}
#else
		int p = 0;
		evolve_one_subpopulation(&p, seeds, sub_populations, generation, &params, original_matrices, mask_matrices, num_training_images, &clip, work_matrix1, work_matrix2);
#endif
		// find the best individual
		best_subpopulation_index = 0; // the index of the subpopulation containing the best invidual
		unsigned long long sum_fitness = 0;
		for (int p = 0; p < params.num_sub_populations; p++) {
		//	printf("--------------------------\n");
			if (sub_populations[p][0].fitness < sub_populations[best_subpopulation_index][0].fitness)
				best_subpopulation_index = p;
			for (int q = 0; q < params.sub_population_size; q++) {
				sum_fitness += sub_populations[p][q].fitness;
		//		printf("%ld ", sub_populations[p][q].fitness);
			}
		}
		sum_fitness /= params.num_sub_populations * params.sub_population_size;
		printf("generation = %d, best fitness = %ld, average fitness = %ld\n", generation, sub_populations[best_subpopulation_index][0].fitness, sum_fitness);

		if (best_fitness_so_far > sub_populations[best_subpopulation_index][0].fitness) {
			best_fitness_so_far = sub_populations[best_subpopulation_index][0].fitness;

			//sub_populations[best_subpopulation_index][0].simplify(params.code_length);

			// debug only
#ifdef USE_THREADS
			//fitness(sub_populations[best_subpopulation_index][0], params.code_length, original_matrices, mask_matrices, num_training_images, clip, work_matrices1[0], work_matrices2[0]);
#else
			//fitness(sub_populations[best_subpopulation_index][0], params.code_length, original_matrices, mask_matrices, num_training_images, clip, work_matrix1, work_matrix2);
#endif

			char file_name[50];
			for (int c = 0; c < num_images; c++) {
				sprintf(file_name, "ca_%u_%d_%d.bmp", initial_seed, generation, c);
 				test_ca(sub_populations[best_subpopulation_index][0], params.code_length, original_matrices[c], image_width, image_height, file_name);
			}
			sprintf(file_name, "ca_%u_%d.txt", initial_seed, generation);
			print_chromosome(sub_populations[best_subpopulation_index][0], params, file_name);
		}
		// now copy one individual from one population to the next one.
		// the copied invidual will replace the worst in the next one (if is better)

		for (int p = 0; p < params.num_sub_populations; p++) {
			int  k = my_int_rand(seeds[p], 0, params.sub_population_size - 1);// the individual to be copied
			// replace the worst in the next population (p + 1) - only if is better
			int index_next_pop = (p + 1) % params.num_sub_populations; // index of the next subpopulation (taken in circular order)
			if (sub_populations[p][k].fitness < sub_populations[index_next_pop][params.sub_population_size - 1].fitness) {
				copy_individual(sub_populations[index_next_pop][params.sub_population_size - 1], sub_populations[p][k], params);
				qsort((void *)sub_populations[index_next_pop], params.sub_population_size, sizeof(sub_populations[0][0]), sort_function);
			}
		}

		save_state(params.num_sub_populations, params.sub_population_size, params.code_length, seeds, sub_populations, generation, initial_seed);
	}

#ifdef USE_THREADS
	delete[] mep_threads;
#endif

	// print best chromosome

	print_chromosome(sub_populations[best_subpopulation_index][0], params, NULL);

	// free memory

	for (int p = 0; p < params.num_sub_populations; p++) {
		for (int i = 0; i < params.sub_population_size; i++)
			delete_chromosome(sub_populations[p][i]);
		delete[] sub_populations[p];
	}
	delete[] sub_populations;

#ifdef USE_THREADS
	delete_matrices(work_matrices1, clip.height + 2 * CA_radius, params.num_threads);
	delete_matrices(work_matrices2, clip.height + 2 * CA_radius, params.num_threads);
#else
	delete_matrix(work_matrix1, clip.height + 2 * CA_radius);
	delete_matrix(work_matrix2, clip.height + 2 * CA_radius);
#endif
}
//--------------------------------------------------------------------
bool read_image(const char* file_name, t_rgb **&matrix, unsigned long& image_width, unsigned long &image_height)
{
	BMP *bmp;
	
	bmp = BMP_ReadFile(file_name);

	if (BMP_GetError() != BMP_OK) {
		//printf("An error has occurred: %s (code %d)\n", BMP_GetErrorDescription(), BMP_GetError());
		return false;
	}

	image_width = BMP_GetWidth(bmp);
	image_height = BMP_GetHeight(bmp);

	allocate_matrix(matrix, image_width, image_height);

	for (unsigned long x = 0; x < image_width; x++)
		for (unsigned long y = 0; y < image_height; y++) {
			BMP_GetPixelRGB(bmp, x, y, &matrix[y][x].red, &matrix[y][x].green, &matrix[y][x].blue);
		}
	BMP_Free(bmp);

	return true;
}
//--------------------------------------------------------------------
bool read_input_images(char *path_to_images, int &num_images, t_rgb ***&original_matrix, t_rgb ***& mask_matrix, unsigned long& image_width, unsigned long &image_height)
{
	num_images = 16;

	original_matrix = new t_rgb**[num_images];
	mask_matrix = new t_rgb**[num_images];
	char bmp_file_name[256];


	int t = 0;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "3063_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "3063_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "108004_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "108004_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "3096_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "3096_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "8068_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "8068_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;
	
	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "69040_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "69040_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;

	t++;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "69020_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "69020_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "100098_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "100098_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;
	
////////////////

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "106047_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "106047_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "130026_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "130026_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;


	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "134008_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "134008_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "135037_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "135037_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "304034_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "304034_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "43051_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "43051_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "69022_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "69022_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;


	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "70011_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "70011_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "80099_original.bmp");
	if (!read_image(bmp_file_name, original_matrix[t], image_width, image_height))
		return false;

	strcpy(bmp_file_name, path_to_images);
	strcat(bmp_file_name, "80099_mask.bmp");
	if (!read_image(bmp_file_name, mask_matrix[t], image_width, image_height))
		return false;
	t++;

	return true; // OK
}
//--------------------------------------------------------------------
int main(int argc, char** argv)
{
	t_parameters params;
	params.num_sub_populations = 6;
	params.sub_population_size = 10;				// the number of individuals in population  (must be an even number!)
	params.code_length = 30;
	params.num_generations = 100000;				// the number of generations
	params.mutation_probability = 0.01;             // mutation probability
	params.crossover_probability = 0.9;             // crossover probability

	params.variables_probability = 0.5;
	params.operators_probability = 0.5;
	params.constants_probability = 1 - params.variables_probability - params.operators_probability; // sum of variables_prob + operators_prob + constants_prob MUST BE 1 !

	params.num_constants = 0;                       // use constants from [0..255] interval
	params.constants_min = 0;
	params.constants_max = 255;

#ifdef USE_THREADS
	params.num_threads = 6;
#endif


	unsigned long image_width, image_height; // all images are assumed to have the same size.
	t_rgb ***original_matrix;
	t_rgb ***mask_matrix;

	int num_images = 0;

	char path_to_images[256];
	strcpy(path_to_images, "c:\\Mihai\\Dropbox\\ca-segmentation\\images\\"); // must be set by the user

	if (!read_input_images(path_to_images, num_images, original_matrix, mask_matrix, image_width, image_height)) {
		printf("Cannot find input files! Please specify the full path!");
		getchar();
		return 1;
	}

	// working on a clip, not on the full image
	t_clip_region clip;
	// clip MUST contain both object and background
	clip.top_x = 0;
	clip.top_y = 120;
	clip.width = image_width;
	clip.height = 30; // image_height;

	char save_state_filename[1000];
	save_state_filename[0] = 0;
	t_seed* seeds = new t_seed[params.num_sub_populations];


	time_t start_time;
	time(&start_time);

	// allocate memory for all sub populations
	t_chromosome **sub_populations; // an array of sub populations
	sub_populations = new t_chromosome*[params.num_sub_populations];
	for (int p = 0; p < params.num_sub_populations; p++) {
		sub_populations[p] = new t_chromosome[params.sub_population_size];
		for (int i = 0; i < params.sub_population_size; i++)
			allocate_chromosome(sub_populations[p][i], params); // allocate each individual in the subpopulation
	}

	int current_generation_index;

	uint32_t initial_seed = 0;
	if (argc == 1) {
		initial_seed = 10000; // must be greater than 127

		for (int p = 0; p < params.num_sub_populations; p++) {
			seeds[p].init(initial_seed, p);
			for (int c = 0; c < params.sub_population_size; c++)
				generate_random_chromosome(sub_populations[p][c], params, seeds[p]);
		}
	}
	else// 2 arguments
		// it can be either an int (initial seed)
		// or a string (state)
		if (argv[1][0] == 's') {// state
			if (!load_state("state_x_x.txt", seeds, sub_populations, current_generation_index)) {
				printf("Cannot load state. Press Enter!");
				getchar();
			}
		}
		else {// seed
			initial_seed = atoi(argv[1]);

			for (int p = 0; p < params.num_sub_populations; p++) {
				seeds[p].init(initial_seed, p);
				for (int c = 0; c < params.sub_population_size; c++) 
				generate_random_chromosome(sub_populations[p][c], params, seeds[p]);
			}
		}

	printf("evolving...\n");
	start_steady_state_gp(sub_populations, params, initial_seed, seeds, original_matrix, mask_matrix, num_images, clip, image_width, image_height);

	time_t end_time;
	time(&end_time);

	double running_time = difftime(end_time, start_time);

	printf("running time = %lf\n", running_time);

	delete_matrices(original_matrix, image_height, num_images);
	delete_matrices(mask_matrix, image_height, num_images);

	delete[] seeds;


	printf("Press enter ...");
	getchar();

	return 0;
}
//--------------------------------------------------------------------