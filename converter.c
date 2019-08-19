#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

struct color {
	double value[3];
};

struct color cartesian_to_polar(struct color source_color) {
	struct color target_color;
	target_color.value[0] = source_color.value[0];
	target_color.value[1] = hypot(source_color.value[1], source_color.value[2]);
	target_color.value[2] = atan2(source_color.value[2], source_color.value[1])*180.0/M_PI;

	if (target_color.value[2] < 0) {
		target_color.value[2] = target_color.value[2] + 360;
	}

	return target_color;
}

struct color polar_to_cartesian(struct color source_color) {
	struct color target_color;
	double angle;
	angle = source_color.value[2]/180*(M_PI);
	target_color.value[0] = source_color.value[0];
	target_color.value[1] = cos(angle)*source_color.value[1];
	target_color.value[2] = sin(angle)*source_color.value[1];
	return target_color;
}


/* XYZ = CIE 1931 XYZ */

struct color XYZ_values_to_OSA74(double X, double Y, double Z) {
	struct color target_color;
	double x = X/(X+Y+Z);
	double y = Y/(X+Y+Z);
	double Y0 = (4.4934*x*x + 4.3034*y*y - 4.276*x*y - 1.3744*x - 2.5643*y + 1.8103)*Y;
	double Lambda = 5.9*(cbrt(Y0) - 2.0/3 + 0.042*cbrt(Y0-30));
	double C = Lambda/(5.9*(cbrt(Y0) - 2.0/3));
	double R = 0.7990*X + 0.4194*Y - 0.1648*Z;
	double G = -0.4493*X + 1.3265*Y + 0.0927*Z;
	double B = -0.1149*X + 0.3394*Y + 0.7170*Z;
	/* target_color.value[0] = (Lambda - 14.3993)/sqrt(2); */
	target_color.value[0] = (Lambda - 14.4)/sqrt(2);
	target_color.value[1] = C*(1.7*cbrt(R) + 8*cbrt(G) - 9.7*cbrt(B));
	target_color.value[2] = C*(-13.7*cbrt(R) + 17.7*cbrt(G) - 4*cbrt(B));
	return target_color;
}

struct color XYZ_values_to_OSA90(double X, double Y, double Z) {
	struct color target_color;
	double x = X/(X+Y+Z);
	double y = Y/(X+Y+Z);
	double Y0 = (4.4934*x*x + 4.3034*y*y - 4.276*x*y - 1.3744*x - 2.5643*y + 1.8103)*Y;
	double Lambda = 5.9*(cbrt(Y0) - 2.0/3 + 0.06*cbrt(Y0-30));
	double C = Lambda/(5.9*(cbrt(Y0) - 2.0/3));
	double R = 0.9285*X + 0.3251*Y - 0.1915*Z;
	double G = -0.4493*X + 1.3265*Y + 0.0927*Z;
	double B = -0.2032*X + 0.6*Y + 0.5523*Z;
	/* target_color.value[0] = (Lambda - 14.3993)/sqrt(2); */
	target_color.value[0] = (Lambda - 14.4)/sqrt(2);
	target_color.value[1] = C*(-1.3*cbrt(R) + 17*cbrt(G) - 15.7*cbrt(B));
	target_color.value[2] = C*(-12.7*cbrt(R) + 19*cbrt(G) - 6.3*cbrt(B));
	return target_color;
}

struct color XYZ_to_OSA74(struct color source_color) {
	double X, Y, Z;
	X = source_color.value[0];
	Y = source_color.value[1];
	Z = source_color.value[2];
	return XYZ_values_to_OSA74(X, Y, Z);
}

struct color XYZ_to_OSA90(struct color source_color) {
	double X, Y, Z;
	X = source_color.value[0];
	Y = source_color.value[1];
	Z = source_color.value[2];
	return XYZ_values_to_OSA90(X, Y, Z);
}

int print_matrix(gsl_matrix *matrix) {
	int row, column;
	for(row = 0; row < matrix->size1; row++) {
		for(column = 0; column < matrix->size2; column++) {
			printf("%8.8f", gsl_matrix_get(matrix, row, column));
		}
		putchar('\n');
	}
	return 0;
}

struct color OSA74_to_XYZ(struct color source_color) {
	double L, j, g;
	L = source_color.value[0];
	j = source_color.value[1];
	g = source_color.value[2];
	struct color target_color;
	int max_iterations = 600;
	double X, Y, Z, previous_X, previous_Y, previous_Z, Delta;
	double L0, L1, L2, L3, j0, j1, j2, j3, g0, g1, g2, g3;
	struct color color0, color1, color2, color3;

	gsl_matrix *matrix0 = gsl_matrix_alloc(3, 3);
	gsl_matrix *inverse_matrix0 = gsl_matrix_alloc(3, 3);
	gsl_matrix *matrix1 = gsl_matrix_alloc(3, 1);
	gsl_matrix *matrix2 = gsl_matrix_alloc(3, 1);
	gsl_permutation *permutation = gsl_permutation_alloc(3);
	int sign = 0;

	/* start values */
	/*X = 28.4; Y = 30.0; Z = 32.2; Delta = 0.5;*/
	X = 28.4; Y = 30.0; Z = 32.2; Delta = 0.005;
	gsl_matrix_set(matrix2, 0, 0, X);
	gsl_matrix_set(matrix2, 1, 0, Y);
	gsl_matrix_set(matrix2, 2, 0, Z);

	int i = 0; do {
		/* calculate colors for influence matrix */
		color0 = XYZ_values_to_OSA74(X, Y, Z);
		color1 = XYZ_values_to_OSA74(X+Delta, Y, Z);
		color2 = XYZ_values_to_OSA74(X, Y+Delta, Z);
		color3 = XYZ_values_to_OSA74(X, Y, Z+Delta);
		L0 = color0.value[0]; j0 = color0.value[1]; g0 = color0.value[2];
		L1 = color1.value[0]; j1 = color1.value[1]; g1 = color1.value[2];
		L2 = color2.value[0]; j2 = color2.value[1]; g2 = color2.value[2];
		L3 = color3.value[0]; j3 = color3.value[1]; g3 = color3.value[2];

		gsl_matrix_set(matrix0, 0, 0, (L1-L0)/Delta);
		gsl_matrix_set(matrix0, 0, 1, (L2-L0)/Delta);
		gsl_matrix_set(matrix0, 0, 2, (L3-L0)/Delta);
		gsl_matrix_set(matrix0, 1, 0, (j1-j0)/Delta);
		gsl_matrix_set(matrix0, 1, 1, (j2-j0)/Delta);
		gsl_matrix_set(matrix0, 1, 2, (j3-j0)/Delta);
		gsl_matrix_set(matrix0, 2, 0, (g1-g0)/Delta);
		gsl_matrix_set(matrix0, 2, 1, (g2-g0)/Delta);
		gsl_matrix_set(matrix0, 2, 2, (g3-g0)/Delta);

		gsl_matrix_set(matrix1, 0, 0, L-L0);
		gsl_matrix_set(matrix1, 1, 0, j-j0);
		gsl_matrix_set(matrix1, 2, 0, g-g0);

		/* update matrix2:  matrix2 = inverse_matrix0  matrix1 + matrix2 */
		/* incremental term (inverse_matrix0 × matrix1) is divided by (5+cbrt(i)) */
		gsl_linalg_LU_decomp(matrix0, permutation, &sign);
		gsl_linalg_LU_invert(matrix0, permutation, inverse_matrix0);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1/(5+cbrt(i)), inverse_matrix0, matrix1, 1.0, matrix2);

		previous_X = X; previous_Y = Y; previous_Z = Z;
		X = gsl_matrix_get(matrix2, 0, 0);
		Y = gsl_matrix_get(matrix2, 1, 0);
		Z = gsl_matrix_get(matrix2, 2, 0);

		/* printf("\n%d\n", i);
		print_matrix(matrix2); */

		i++;
	} while ((i < max_iterations) && ((fabs(X - previous_X) > 0.0001) || (fabs(Y - previous_Y) > 0.0001) || (fabs(Z - previous_Z) > 0.0001)));

	/*printf("Iterations: %i\n", i);*/

	/*free memory*/
	gsl_matrix_free(matrix0);
	gsl_matrix_free(inverse_matrix0);
	gsl_matrix_free(matrix1);
	gsl_matrix_free(matrix2);
	/*gsl_matrix_free(matrix_product);*/
	gsl_permutation_free(permutation);

	if (i == max_iterations) {
		fprintf(stderr, "Error: Maximum iterations reached.\n");
	}

	target_color.value[0] = X;
	target_color.value[1] = Y;
	target_color.value[2] = Z;
	return target_color;
}

struct color OSA90_to_XYZ(struct color source_color) {
	double L, j, g;
	L = source_color.value[0];
	j = source_color.value[1];
	g = source_color.value[2];
	struct color target_color;
	int max_iterations = 600;
	double X, Y, Z, previous_X, previous_Y, previous_Z, Delta;
	double L0, L1, L2, L3, j0, j1, j2, j3, g0, g1, g2, g3;
	struct color color0, color1, color2, color3;

	gsl_matrix *matrix0 = gsl_matrix_alloc(3, 3);
	gsl_matrix *inverse_matrix0 = gsl_matrix_alloc(3, 3);
	gsl_matrix *matrix1 = gsl_matrix_alloc(3, 1);
	gsl_matrix *matrix2 = gsl_matrix_alloc(3, 1);
	gsl_permutation *permutation = gsl_permutation_alloc(3);
	int sign = 0;

	/* start values */
	/*X = 28.4; Y = 30.0; Z = 32.2; Delta = 0.5;*/
	X = 28.4; Y = 30.0; Z = 32.2; Delta = 0.005;
	gsl_matrix_set(matrix2, 0, 0, X);
	gsl_matrix_set(matrix2, 1, 0, Y);
	gsl_matrix_set(matrix2, 2, 0, Z);

	int i = 0; do {
		/* calculate colors for influence matrix */
		color0 = XYZ_values_to_OSA90(X, Y, Z);
		color1 = XYZ_values_to_OSA90(X+Delta, Y, Z);
		color2 = XYZ_values_to_OSA90(X, Y+Delta, Z);
		color3 = XYZ_values_to_OSA90(X, Y, Z+Delta);
		L0 = color0.value[0]; j0 = color0.value[1]; g0 = color0.value[2];
		L1 = color1.value[0]; j1 = color1.value[1]; g1 = color1.value[2];
		L2 = color2.value[0]; j2 = color2.value[1]; g2 = color2.value[2];
		L3 = color3.value[0]; j3 = color3.value[1]; g3 = color3.value[2];

		gsl_matrix_set(matrix0, 0, 0, (L1-L0)/Delta);
		gsl_matrix_set(matrix0, 0, 1, (L2-L0)/Delta);
		gsl_matrix_set(matrix0, 0, 2, (L3-L0)/Delta);
		gsl_matrix_set(matrix0, 1, 0, (j1-j0)/Delta);
		gsl_matrix_set(matrix0, 1, 1, (j2-j0)/Delta);
		gsl_matrix_set(matrix0, 1, 2, (j3-j0)/Delta);
		gsl_matrix_set(matrix0, 2, 0, (g1-g0)/Delta);
		gsl_matrix_set(matrix0, 2, 1, (g2-g0)/Delta);
		gsl_matrix_set(matrix0, 2, 2, (g3-g0)/Delta);

		gsl_matrix_set(matrix1, 0, 0, L-L0);
		gsl_matrix_set(matrix1, 1, 0, j-j0);
		gsl_matrix_set(matrix1, 2, 0, g-g0);

		/* update matrix2:  matrix2 = inverse_matrix0  matrix1 + matrix2 */
		/* incremental term (inverse_matrix0 × matrix1) is divided by (5+cbrt(i)) */
		gsl_linalg_LU_decomp(matrix0, permutation, &sign);
		gsl_linalg_LU_invert(matrix0, permutation, inverse_matrix0);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1/(5+cbrt(i)), inverse_matrix0, matrix1, 1.0, matrix2);

		previous_X = X; previous_Y = Y; previous_Z = Z;
		X = gsl_matrix_get(matrix2, 0, 0);
		Y = gsl_matrix_get(matrix2, 1, 0);
		Z = gsl_matrix_get(matrix2, 2, 0);

		/* printf("\n%d\n", i);
		print_matrix(matrix2); */

		i++;
	} while ((i < max_iterations) && ((fabs(X - previous_X) > 0.0001) || (fabs(Y - previous_Y) > 0.0001) || (fabs(Z - previous_Z) > 0.0001)));

	/*printf("Iterations: %i\n", i);*/

	/*free memory*/
	gsl_matrix_free(matrix0);
	gsl_matrix_free(inverse_matrix0);
	gsl_matrix_free(matrix1);
	gsl_matrix_free(matrix2);
	/*gsl_matrix_free(matrix_product);*/
	gsl_permutation_free(permutation);

	if (i == max_iterations) {
		fprintf(stderr, "Error: Maximum iterations reached.\n");
	}

	target_color.value[0] = X;
	target_color.value[1] = Y;
	target_color.value[2] = Z;
	return target_color;
}

double npow(double base, double exponent) {
	if (base < 0) {
		return -pow(fabs(base), exponent);
	} else {
		return pow(base, exponent);
	}
}

const double c1 = 3424/pow(2, 12);
const double c2 = 2413/pow(2, 7);
const double c3 = 2392/pow(2, 7);
const double n = 2610/pow(2, 14);
const double p = 1.7*2523/pow(2, 5);

double perceptual_quantizer(double cone_response) {
	return npow((c1 + c2*npow(cone_response/10000, n)) / (1 + c3*npow(cone_response/10000, n)), p);
}

double invert_perceptual_quantizer(double cone_response) {
	/*printf("input: %f\n", cone_response);*/
	/*printf("return: %f\n", 10000*npow((c1 - npow(cone_response, 1/p))/(c3*npow(cone_response, 1/p) - c2), 1/n));*/

	return 10000*npow((c1 - npow(cone_response, 1/p))/(c3*npow(cone_response, 1/p) - c2), 1/n);
}

const double d = -0.56;
const double d0 = 1.6295499532821566E-11;

struct color XYZ_to_Jab(struct color source_color) {
	double X = source_color.value[0];
	double Y = source_color.value[1];
	double Z = source_color.value[2];
	double adjusted_X, adjusted_Y, L, M, S, adjusted_L, adjusted_M, adjusted_S, I, J, a, b, C, h;
	struct color target_color;

	adjusted_X = 1.15*X - 0.15*Z;
	adjusted_Y = 0.66*Y + 0.34*X;

	L = 0.41478972*adjusted_X + 0.579999*adjusted_Y + 0.0146480*Z;
	M = -0.2015100*adjusted_X + 1.120649*adjusted_Y + 0.0531008*Z;
	S = -0.0166008*adjusted_X + 0.264800*adjusted_Y + 0.6684799*Z;

	adjusted_L = perceptual_quantizer(L);
	adjusted_M = perceptual_quantizer(M);
	adjusted_S = perceptual_quantizer(S);

	I = 0.5*adjusted_L + 0.5*adjusted_M;
	J = ((1+d)*I)/(1+d*I) - d0;

	a = 3.524000*adjusted_L - 4.066708*adjusted_M + 0.542708*adjusted_S;
	b = 0.199076*adjusted_L + 1.096799*adjusted_M - 1.295875*adjusted_S;

	target_color.value[0] = J;
	target_color.value[1] = a;
	target_color.value[2] = b;
	return target_color;
}

struct color Jab_to_XYZ(struct color source_color) {
	double J = source_color.value[0];
	double a = source_color.value[1];
	double b = source_color.value[2];
	double I, adjusted_L, adjusted_M, adjusted_S, L, M, S, adjusted_X, adjusted_Y, X, Y, Z;
	struct color target_color;

	I = ((J + d0)/(1 + d - d*(J + d0)));

	adjusted_L = I + 0.138605*a + 0.0580473*b;
	adjusted_M = I - 0.138605*a - 0.0580473*b;
	adjusted_S = I - 0.0960192*a - 0.811892*b;

	L = invert_perceptual_quantizer(adjusted_L);
	M = invert_perceptual_quantizer(adjusted_M);
	S = invert_perceptual_quantizer(adjusted_S);

	adjusted_X = 1.92423*L - 1.00479*M + 0.0376514*S;
	adjusted_Y = 0.350317*L + 0.726481*M - 0.0653844*S;
	Z = -0.0909828*L - 0.312728*M + 1.52277*S;

	X = (adjusted_X + 0.15*Z)/1.15;
	Y = (adjusted_Y - 0.34*X)/0.66;

	target_color.value[0] = X;
	target_color.value[1] = Y;
	target_color.value[2] = Z;
	return target_color;
}

const double J_scale = 0.167174631034780;

struct color scale_Jab(struct color source_color) {
	struct color target_color;

	target_color.value[0] = source_color.value[0]/J_scale*100;
	target_color.value[1] = source_color.value[1]/J_scale*100;
	target_color.value[2] = source_color.value[2];
	return target_color;
}

struct color unscale_Jab(struct color source_color) {
	struct color target_color;

	target_color.value[0] = source_color.value[0]*J_scale/100;
	target_color.value[1] = source_color.value[1]*J_scale/100;
	target_color.value[2] = source_color.value[2];
	return target_color;
}

double apply_gamma_correction(double value) {
	if (value > 0.0031308) {
		return 1.055*pow(value, (1.0/2.4)) - 0.055;
	} else {
		return 12.92*value;
	}
}

double remove_gamma_correction(double value) {
	if (value > 0.04045) {
		return pow(((value+0.055)/1.055), 2.4);
	} else {
		return value/12.92;
	}
}

struct color XYZ_to_sRGB(struct color source_color) {
	double X, Y, Z;
	X = source_color.value[0];
	Y = source_color.value[1];
	Z = source_color.value[2];
	struct color target_color;
	double linear_R, linear_G, linear_B, R, G, B;

	/* XYZ source values range = [0, 100] */
	X = X/100.0; Y = Y/100.0; Z = Z/100.0;

	/* 2° observer, D65 illuminant */
	linear_R = X* 3.2406 + Y*-1.5372 + Z*-0.4986;
	linear_G = X*-0.9689 + Y* 1.8758 + Z* 0.0415;
	linear_B = X* 0.0557 + Y*-0.2040 + Z* 1.0570;

	R = apply_gamma_correction(linear_R);
	G = apply_gamma_correction(linear_G);
	B = apply_gamma_correction(linear_B);

	if ((R<0.0) || (R>1.0) || (G<0.0) || (G>1.0) || (B<0.0) || (B>1.0)) {
		errno = ERANGE;
	}

	/* RGB target values range = [0, 255] */
	target_color.value[0] = R*255;
	target_color.value[1] = G*255;
	target_color.value[2] = B*255;

	return target_color;
}

struct color sRGB_to_XYZ(struct color source_color) {
	double R, G, B;
	R = source_color.value[0];
	G = source_color.value[1];
	B = source_color.value[2];
	struct color target_color;
	double linear_R, linear_G, linear_B, X, Y, Z;

	/* RGB source values range = [0, 255] */
	R = R/255.0; G = G/255.0; B = B/255.0;

	linear_R = remove_gamma_correction(R);
	linear_G = remove_gamma_correction(G);
	linear_B = remove_gamma_correction(B);

	/* 2° observer, D65 illuminant */
	X = linear_R*0.4124 + linear_G*0.3576 + linear_B*0.1805;
	Y = linear_R*0.2126 + linear_G*0.7152 + linear_B*0.0722;
	Z = linear_R*0.0193 + linear_G*0.1192 + linear_B*0.9505;

	/* XYZ target values range = [0, 100] */
	target_color.value[0] = 100*X;
	target_color.value[1] = 100*Y;
	target_color.value[2] = 100*Z;
	return target_color;
}

const double Y_threshold = pow((6/29.0), 3);
const double D65_X = 95.0156;
const double D65_Y = 100;
const double D65_Z = 108.8199;

/*const double D65_X = 97.31273;*/
/*const double D65_Y = 100.00;*/
/*const double D65_Z = 138.59590;*/


struct color XYZ_to_Luv(struct color source_color) {
	/* L, u, v stand for L*, u*, v* */
	double X, Y, Z, prime_u, prime_v, relative_Y, white_prime_u, white_prime_v, L, u, v;
	struct color target_color;

	X = source_color.value[0];
	Y = source_color.value[1];
	Z = source_color.value[2];

	prime_u = (4*X) / (X + 15*Y + 3*Z);
	prime_v = (9*Y) / (X + 15*Y + 3*Z);

	white_prime_u = (4*D65_X) / (D65_X + 15*D65_Y + 3*D65_Z);
	white_prime_v = (9*D65_Y) / (D65_X + 15*D65_Y + 3*D65_Z);

	relative_Y = Y/D65_Y;

	if (relative_Y > pow((6/29.0), 3)) {
		L = cbrt(relative_Y)*116 - 16;
	} else {
		L = relative_Y*pow((29/3.0), 3);
	}

	u = 13*L*(prime_u - white_prime_u);
	v = 13*L*(prime_v - white_prime_v);

	target_color.value[0] = L;
	target_color.value[1] = u;
	target_color.value[2] = v;
	return target_color;
}

struct color Luv_to_XYZ(struct color source_color) {
	/* L, u, v stand for L*, u*, v* */
	double L, u, v, white_prime_u, white_prime_v, prime_u, prime_v, X, Y, Z;
	struct color target_color;

	L = source_color.value[0];
	u = source_color.value[1];
	v = source_color.value[2];

	white_prime_u = (4*D65_X) / (D65_X + 15*D65_Y + 3*D65_Z);
	white_prime_v = (9*D65_Y) / (D65_X + 15*D65_Y + 3*D65_Z);

	prime_u = u/(13*L) + white_prime_u;
	prime_v = v/(13*L) + white_prime_v;

	if (L > 8.0) {
		Y = D65_Y*pow((L + 16)/116.0, 3);
	} else {
		Y = D65_Y*L*pow(3/29.0, 3);
	}

	X = Y * (9*prime_u)/(4.0*prime_v);
	Z = Y * (12 - 3*prime_u - 20*prime_v)/(4.0*prime_v);

	target_color.value[0] = X;
	target_color.value[1] = Y;
	target_color.value[2] = Z;
	return target_color;
}

/* OSA74 range: -19.1032--10.0938 (29.1970) */

struct color OSA74_to_pOSA74(struct color source_color) {
	double L, j, g, new_L, chroma, hue;
	struct color target_color;

	L = source_color.value[0];
	j = source_color.value[1];
	g = source_color.value[2];

	new_L = (L*sqrt(2) + 19.1032)*100/29.1970;
	chroma = hypot(j, g)*100/29.1970;
	hue = atan2(g, j)*180.0/M_PI;

	if (hue < 0) {
		hue = hue + 360;
	}

	target_color.value[0] = new_L;
	target_color.value[1] = chroma;
	target_color.value[2] = hue;
	return target_color;
}

struct color pOSA74_to_OSA74(struct color source_color) {
	double new_L, chroma, hue, L, angle, j, g;
	struct color target_color;

	new_L = source_color.value[0];
	chroma = source_color.value[1];
	hue = source_color.value[2];

	L = (new_L*29.1970/100 - 19.1032)/sqrt(2);

	angle = hue/180*(M_PI);
	j = cos(angle)*chroma*29.1970/100;
	g = sin(angle)*chroma*29.1970/100;

	target_color.value[0] = L;
	target_color.value[1] = j;
	target_color.value[2] = g;
	return target_color;
}

int RGB_to_hex(struct color RGB_color) {
	int R = (int)(RGB_color.value[0] + 0.5);
	int G = (int)(RGB_color.value[1] + 0.5);
	int B = (int)(RGB_color.value[2] + 0.5);
	return (R<<16) | (G<<8) | B;
}

struct color hex_to_RGB(int hex_value) {
	struct color RGB_color;
	RGB_color.value[0] = ((hex_value>>16) & 0xFF)/255.0;
	RGB_color.value[1] = ((hex_value>>8) & 0xFF)/255.0;
	RGB_color.value[2] = (hex_value & 0xFF)/255.0;
	return RGB_color;
}


int main(int argc, char *argv[]) {
	int i;
	struct color source_color;
	struct color intermediate_color;
	struct color target_color;
	char next_char;
	char *source_format = argv[1];
	char *target_format = argv[2];

	/* check arguments */

	/* check argument count */
	if (argc != 6) {
		fprintf(stderr, "Usage: %s source_format target_format value1 value2 value3\n", argv[0]);
		return EXIT_FAILURE;
	}

	if (strcmp(source_format, target_format) == 0) {
		fprintf(stderr, "Source and target formats are identical.\n");
		return EXIT_FAILURE;
	}

	/* check if color values are numbers */
	for (i = 3; i < argc; i++) {
		if (sscanf(argv[i], "%lf%c", &source_color.value[i-3], &next_char) != 1) {
			fprintf(stderr, "Invalid color values, values must be numbers.\n");
			return EXIT_FAILURE;
		}
	}

	/* OSA/polar conversions do not need intermediate conversion to XYZ */
	if ((strcmp(source_format, "pOSA74") == 0) && (strcmp(target_format, "OSA74") == 0)) {
		target_color = pOSA74_to_OSA74(source_color);
		printf("%.2f %.2f %.2f\n", target_color.value[0], target_color.value[1], target_color.value[2]);
		return EXIT_SUCCESS;
	} else if ((strcmp(source_format, "OSA74") == 0) && (strcmp(target_format, "pOSA74") == 0)) {
		target_color = OSA74_to_pOSA74(source_color);
		printf("%.2f %.2f %.2f\n", target_color.value[0], target_color.value[1], target_color.value[2]);
		return EXIT_SUCCESS;
	}

	/* convert from source format to XYZ */
	if (strcmp(source_format, "Jab") == 0) {
		intermediate_color = Jab_to_XYZ(polar_to_cartesian(unscale_Jab(source_color)));
		/*intermediate_color = Jab_to_XYZ(polar_to_cartesian(source_color));*/
		/*intermediate_color = Jab_to_XYZ(source_color);*/
	} else if (strcmp(source_format, "LCh") == 0) {
		intermediate_color = Luv_to_XYZ(polar_to_cartesian(source_color));
	} else if (strcmp(source_format, "Luv") == 0) {
		intermediate_color = Luv_to_XYZ(source_color);
	} else if (strcmp(source_format, "OSA74") == 0) {
		intermediate_color = OSA74_to_XYZ(source_color);
	} else if (strcmp(source_format, "OSA90") == 0) {
		intermediate_color = OSA90_to_XYZ(source_color);
	} else if (strcmp(source_format, "pOSA74") == 0) {
		/*intermediate_color = OSA90_to_XYZ(polar_to_OSA(source_color));*/
		intermediate_color = OSA74_to_XYZ(pOSA74_to_OSA74(source_color));
	} else if (strcmp(source_format, "sRGB") == 0) {
		intermediate_color = sRGB_to_XYZ(source_color);
	} else if (strcmp(source_format, "XYZ") == 0) {
		intermediate_color = source_color;
	} else {
		fprintf(stderr, "Invalid source format.");
		return EXIT_FAILURE;
	}

	/* convert from XYZ to target format */
	if (strcmp(target_format, "Jab") == 0) {
		target_color = scale_Jab(cartesian_to_polar(XYZ_to_Jab(intermediate_color)));
		/*target_color = cartesian_to_polar(XYZ_to_Jab(intermediate_color));*/
		/*target_color = XYZ_to_Jab(intermediate_color);*/
	} else if (strcmp(target_format, "LCh") == 0) {
		target_color = cartesian_to_polar(XYZ_to_Luv(intermediate_color));
	} else if (strcmp(target_format, "Luv") == 0) {
		target_color = XYZ_to_Luv(intermediate_color);
	} else if (strcmp(target_format, "OSA74") == 0) {
		target_color = XYZ_to_OSA74(intermediate_color);
	} else if (strcmp(target_format, "OSA90") == 0) {
		target_color = XYZ_to_OSA90(intermediate_color);
	} else if (strcmp(target_format, "pOSA74") == 0) {
		target_color = OSA74_to_pOSA74(XYZ_to_OSA74(intermediate_color));
	} else if (strcmp(target_format, "polarOSA74") == 0) {
		target_color = cartesian_to_polar(XYZ_to_OSA74(intermediate_color));
	} else if (strcmp(target_format, "sRGB") == 0) {
		target_color = XYZ_to_sRGB(intermediate_color);
	} else if (strcmp(target_format, "XYZ") == 0) {
		target_color = intermediate_color;
	} else {
		fprintf(stderr, "Invalid target format.");
		return EXIT_FAILURE;
	}
	if (errno == ERANGE) {
		fprintf(stderr, "The requested color lies outside of the %s gamut.\n", target_format);
		/*return EXIT_FAILURE;*/
		/*printf("%d %d %d\n", 149, 149, 148);*/
		return EXIT_SUCCESS;
	}

	/* round RGB values to integers */
	if (strcmp(target_format, "sRGB") == 0) {
		/*printf("%.0f %.0f %.0f\n", target_color.value[0], target_color.value[1], target_color.value[2]);*/
		printf("#%06x\n", RGB_to_hex(target_color));
	} else {
		printf("%.4f %.4f %.4f\n", target_color.value[0], target_color.value[1], target_color.value[2]);
	}

	return EXIT_SUCCESS;
}
