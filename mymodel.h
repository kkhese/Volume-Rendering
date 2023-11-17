/* IMPORTANT: After you download this file, you should rename it
	to "mymodel.h"
*/

/* Definition of image buffers */
#define SLCS 128
#define ROWS 128
#define COLS 128
unsigned char	CT[SLCS][ROWS][COLS];		/* a 3D array for CT data */
unsigned char	SHADING[SLCS][ROWS][COLS];	/* a 3D array for shading values */

#define IMG_ROWS 512
#define IMG_COLS 512
unsigned char	out_img[IMG_ROWS][IMG_COLS] = { 255, };

/* Camera parameters */
float VRP[3] = { 128.0, 64.0, 250.0 };
float VPN[3] = { -64.0, 0.0, -186.0 };
float VUP[3] = { 0.0, 1.0, 0.0 };

/* Image Plane Sizes */
float focal = 0.05f;	/* 50 mm lens */
float xmin = 0.0175f;	/* 35 mm "film" */
float ymin = 0.0175f;
float xmax = -0.0175f;
float ymax = -0.0175f;

// Saturation value of gradient
float epsillon = 7.25f;

/* Light direction (unit length vector) */
float L[3] = { 0.577f, -0.577f, -0.577f };
/* Light Intensity */
float Ip = 255.0;
// Reflection coeff.
float kd = 1.0f;
/* === transformation matrices (to be constructed) === */
/* Initialize Transformation from the world to the camera coordinates */
float Mwc[4][4] =
{ 1.0, 0.0, 0.0, 0.0,
 0.0, 1.0, 0.0, 0.0,
 0.0, 0.0, 1.0, 0.0,
 0.0, 0.0, 0.0, 1.0 };

/* Initialize Transformation from the camera to the world coordinates */
float Mcw[4][4] =
{ 1.0, 0.0, 0.0, 0.0,
 0.0, 1.0, 0.0, 0.0,
 0.0, 0.0, 1.0, 0.0,
 0.0, 0.0, 0.0, 1.0 };

////////////////////// Fundamental functions were added to header file ////////////////////////
// Obtain magnitude of vecter => returning to nor
float Length(float vector[3])
{
	float abs = vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2];
	float mag = sqrt(abs);
	return mag;
}

// |Vec2 X Vec1| = cross ( cross is normalized vector ) : Order is important!!!!!!!!
// Not |Vec1 X Vec2|, but |Vec2 X Vec1|
void crossX(float vec1[3], float vec2[3], float cross[3]) {
	// cross = Vec2 X Vec 1
	cross[0] = vec2[1] * vec1[2] - vec2[2] * vec1[1];
	cross[1] = (vec2[2] * vec1[0] - vec2[0] * vec1[2]);
	cross[2] = vec2[0] * vec1[1] - vec2[1] * vec1[0];
	// normalized cross, |cross|
	float mag = Length(cross);
	cross[0] = cross[0] / mag;
	cross[1] = cross[1] / mag;
	cross[2] = cross[2] / mag;
}

// out = A X B ( 4by4 X 4by4 matrix calculation )
void matrixcal(float A[4][4], float B[4][4], float out[4][4]) {
	float component = 0.0;
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			component = 0.0;
			for (int k = 0; k < 4; k++) {
				component += A[i][k] * B[k][j];
			}
			out[i][j] = component;
		}
	}
}

// out = A X B ( 4by4 X 1by4 matrix calculation )
void matrixcal2(float A[4][4], float B[4], float out[4]) {
	float component = 0.0;
	for (int i = 0; i < 4; i++) {
		component = 0.0;
		for (int j = 0; j < 4; j++) {
			component += A[i][j] * B[j];
		}
		out[i] = component;
	}
}

// out = transpose matrix of A, used for obtaining inverted matrix
void Transpose(float A[4][4], float out[4][4]) {
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			out[j][i] = A[i][j];
		}
	}
}

// Obtaining u,v,n unit vector & rotation, inverted rotation matrix from 2 vectors (VPN, VUP) 
void uvn(float VUP[3], float VPN[3], float R[4][4], float invR[4][4])
{
	float length = Length(VPN);
	float u[3], v[3], n[3];
	// n vector
	n[0] = VPN[0] / length;
	n[1] = VPN[1] / length;
	n[2] = VPN[2] / length;
	// u and v vectors
	crossX(VPN, VUP, u);
	crossX(u, n, v);
	// Rotation vetor
	for (int i = 0; i < 3; i++)
	{
		R[0][i] = u[i];
		R[1][i] = v[i];
		R[2][i] = n[i];
		R[3][i] = 0.0;
	}
	for (int j = 0; j < 3; j++)
	{
		R[j][3] = 0.0;
	}
	R[3][3] = 1.0;
	// Inverted rotation vector
	Transpose(R, invR);
}
////////////////////////////////// End of Header file //////////////////////////////////////

