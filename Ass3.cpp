// Ass3.cpp : This file contains the 'main' function. Program execution begins and ends there.
// Header Files
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "mymodel.h"
#include <cmath>

using namespace std;

// Compute Shadings for each point
void ComputeShading()
{
	for (int z = 1; z < (SLCS - 1); z++)			// boundary layers are skipped.
		for (int y = 1; y < (ROWS - 1); y++)
			for (int x = 1; x < (COLS - 1); x++) {
				// compute the partial derivative at [z, y, x]
				// Derivative at x, use two unit steps x-1 ~ x+1  
				float deltax = (CT[z][y][x+1] - CT[z][y][x-1])/2.0f;
				// Derivative at y, use two unit steps y-1 ~ y+1  
				float deltay = (CT[z][y+1][x] - CT[z][y-1][x])/2.0f;
				// Derivative at z, use two unit steps z-1 ~ z+1  
				float deltaz = (CT[z+1][y][x] - CT[z-1][y][x])/2.0f;
				// magnitude of the gradient
				float DD = deltax * deltax + deltay * deltay + deltaz * deltaz;
				float D = sqrt(DD);
				// if the magnitude of the gradient is less than a
				// pre-specified epsilon value, set shading to 0.
				if (D < epsillon) {
					SHADING[z][y][x] = (unsigned char) 0;
				}
				// otherwise, compute the diffuse shading.
				// save the result in the shading volume.
				else {
					float N[3] = { deltax / D, deltay / D, deltaz / D };
					SHADING[z][y][x] = (unsigned char) (Ip * kd *( (N[0] * (L[0]/Length(L))) + (N[1] * (L[1]/Length(L))) + (N[2] * (L[2]/Length(L)))) );
				}
			}

}

// Construct ray for each pixel
void RayConstruction(float i, float j, float P0[3], float V[3])
{
	// Xc, Yc Mapping Screen on film 
	float Xc, Yc;
	Xc = ((xmax - xmin) * (j / (IMG_ROWS - 1))) + xmin;
	Yc = ((ymax - ymin) * (i / (IMG_COLS - 1))) + ymin;
	// Origin of World Coordinate
	float Origin[4] = { 0.0, 0.0, 0.0, 1.0 };
	// rawP is world coordinate value for camera origin
	float rawP[4];
	matrixcal2(Mcw, Origin, rawP);
	// Homogenius to 3D vector, P0 is the camera location in world coordinates
	P0[0] = rawP[0]; P0[1] = rawP[1]; P0[2] = rawP[2];
	// one selected point of image plane(film) by i and j
	float rawP1[4] = { Xc, Yc, focal, 1.0 };
	float P1[4] = { 0.0, };
	// transform the selected point of film to world coordinates
	matrixcal2(Mcw, rawP1, P1);
	// Ray vector for selected point by i and j
	float rawV[3] = { 0.0, };
	for (int i = 0; i < 3; i++)
	{
		rawV[i] = P1[i] - P0[i];
	}
	// V is normalized ray vector
	V[0] = rawV[0] / Length(rawV);
	V[1] = rawV[1] / Length(rawV);
	V[2] = rawV[2] / Length(rawV);
}

// Intersection between Ray and Cube
int RayCubeIntersection(float P0[3], float V[3], float ts[2]) {	
	float tempts;
	int n = 0;
	// for plane x=0
	tempts = (float) (-P0[0] / V[0]);
	// check if its also inside y,z range
	if ((P0[1] + (V[1] * tempts) > 0.0f) && (P0[1] + V[1] * tempts < (float)ROWS) && (P0[2] + (V[2] * tempts) > 0.0f) && (P0[2] + V[2] * tempts < (float)SLCS))
	{
		ts[0] = tempts;
		n = n+1;
	}
	// for plane x=127
	tempts = (float)((127.0-P0[0]) / V[0]);
	// check if its also inside y,z range
	if ((P0[1] + (V[1] * tempts) > 0.0f) && (P0[1] + V[1] * tempts < (float)ROWS) && (P0[2] + (V[2] * tempts) > 0.0f) && (P0[2] + V[2] * tempts < (float)SLCS))
	{
		if (ts[0] == 0.0) {
			ts[0] = tempts;
			n = n + 1;
		}
		else {
			ts[1] = tempts;
			n = n + 1;
		}
	}
	// for plane y=0
	tempts = (float)(-P0[1] / V[1]);
	// check if its also inside x,z range
	if ((P0[0] + (V[0] * tempts) > 0.0f) && (P0[0] + V[0] * tempts < (float)COLS) && (P0[2] + (V[2] * tempts) > 0.0f) && (P0[2] + V[2] * tempts < (float)SLCS))
	{
		if (ts[0] == 0.0) {
			ts[0] = tempts;
			n = n + 1;
		}
		else {
			ts[1] = tempts;
			n = n + 1;
		}
	}
	// for plane y=127
	tempts = (float) ((127.0-P0[1]) / V[1]);
	// check if its also inside x,z range
	if ((P0[0] + (V[0] * tempts) > 0.0f) && (P0[0] + V[0] * tempts < (float)COLS) && (P0[2] + (V[2] * tempts) > 0.0f) && (P0[2] + V[2] * tempts < (float)SLCS))
	{
		if (ts[0] == 0.0) {
			ts[0] = tempts;
			n = n + 1;
		}
		else {
			ts[1] = tempts;
			n = n + 1;
		}
	}
	// for plane z=0
	tempts = (float)(-P0[2] / V[2]);
	// check if its also inside x,y range
	if ((P0[0] + (V[0] * tempts) > 0.0f) && (P0[0] + V[0] * tempts < (float)COLS) && (P0[1] + (V[1] * tempts) > 0.0f) && (P0[1] + V[1] * tempts < (float)ROWS))
	{
		if (ts[0] == 0.0) {
			ts[0] = tempts;
			n = n + 1;
		}
		else {
			ts[1] = tempts;
			n = n + 1;
		}
	}
	// for plane z=127
	tempts = (float)((127.0-P0[2]) / V[2]);
	// check if its also inside x,y range
	if ((P0[0] + (V[0] * tempts) > 0.0f) && (P0[0] + V[0] * tempts < (float)COLS) && (P0[1] + (V[1] * tempts) > 0.0f) && (P0[1] + V[1] * tempts < (float)ROWS))
	{
		if (ts[0] == 0.0) {
			ts[0] = tempts;
			n = n + 1;
		}
		else {
			ts[1] = tempts;
			n = n + 1;
		}
	}
	// Make sure ts[0] is smaller than ts[1].
	// If not, swap those.
	if (ts[0] > ts[1])
	{
		float tstemp = ts[1];
		ts[1] = ts[0];
		ts[0] = tstemp;
	}
	return (n);
}

// Function for trilinear value for selected point from 3D array [128][128][128]
float Trilinear(float point[3], unsigned char data[SLCS][ROWS][COLS]) {
	// for x, +y-z side
	float y1z0 = (float) data[(int)floor(point[2])][(int)ceil(point[1])][(int)ceil(point[0])] * (point[0]- floor(point[0]))
		+ (float) data[(int)floor(point[2])][(int)ceil(point[1])][(int)floor(point[0])] * (ceil(point[0]) - point[0]);
	// for x, -y-z side
	float y0z0 = (float) data[(int)floor(point[2])][(int)floor(point[1])][(int)ceil(point[0])] * (point[0] - floor(point[0]))
		+ (float) data[(int)floor(point[2])][(int)floor(point[1])][(int)floor(point[0])] * (ceil(point[0]) - point[0]);
	// for x, +y+z side
	float y1z1 = (float) data[(int)ceil(point[2])][(int)ceil(point[1])][(int)ceil(point[0])] * (point[0] - floor(point[0]))
		+ (float) data[(int)ceil(point[2])][(int)ceil(point[1])][(int)floor(point[0])] * (ceil(point[0]) - point[0]);
	// for x, -y+z side
	float y0z1 = (float) data[(int)ceil(point[2])][(int)floor(point[1])][(int)ceil(point[0])] * (point[0] - floor(point[0]))
		+ (float) data[(int)ceil(point[2])][(int)floor(point[1])][(int)floor(point[0])] * (ceil(point[0]) - point[0]);
	// x was cosidered now for y
	// -z side
	float z0 = y1z0 * (point[1] - floor(point[1])) + y0z0 * (ceil(point[1]) - point[1]);
	// +z side
	float z1 = y1z1 * (point[1] - floor(point[1])) + y0z1 * (ceil(point[1]) - point[1]);
	// x and y were considered now for z
	float result = z1 * (point[2] - floor(point[2])) + z0 * (ceil(point[2]) - point[2]);
	return (result);
}

// Volume rendering
unsigned char VolumeRayTracing(float VRP[3], float V[3], float ts[2])
{
	float Dt = 1.0;		// the interval for sampling along the ray
	float C = 0.0;		// for accumulating the shading value
	float T = 1.0;		// for accumulating the transparency
	/* Marching through the CT volume from t0 to t1 by step size Dt.*/
	if (ts[0] > 0.0) {
		for (float t = ts[0]; t <= ts[1]; t += Dt) {
			// front-to-back order
			// Compute the 3D coordinates of the current sample position in the volume:
			float Coorx = VRP[0] + V[0] * t;
			float Coory = VRP[1] + V[1] * t;
			float Coorz = VRP[2] + V[2] * t;
			float point[3] = { Coorx, Coory, Coorz };
			// Density trilinear interpolation
			float den = Trilinear(point, CT);
			// Shading trilinear interpolation
			float shad = Trilinear(point, SHADING);
			/* Accumulate the shading values in the front-to-back order.
			Note: You will accumulate the transparency. This value
			can be used in the for-loop for early termination. */
			// Density scaled from 0-255 to 0-1.0
			C += (float) (T * (den/255.0) * (shad));
			T *= (float) (1.0 - (den/255.0));
			// If transparency becomes below 0, keep it as 0 for preventing overflow issue
			if (T <= 0.0) {
				T = 0.0;
			}
			// If shading value exceed 255, keep it as 255.
			if (C > 255.0)
			{
				C = 255.0;
			}

		}
	}
	return ((unsigned char) (int) (C));
}


/////////////////////////////////////////////////////////////////////////////////
////////////				Main Program Starts				/////////////////////
/////////////////////////////////////////////////////////////////////////////////
int main()
{
	FILE* fin, * fout; /* input and output file id¡¯s */
	int n=0;

	/* Load the CT data into the array */
	if (fopen_s(&fin, "./smallHead.den", "rb") != NULL) {
		printf("Open CT DATA File Error.\n");
	}
	for (int i = 0; i < SLCS; i++) { /* Read one slice at a time. */
		n = fread(&CT[i][0][0], sizeof(char), ROWS * COLS, fin);
		if (n < ROWS * COLS * sizeof(char)) {
			printf("Read CT data slice %d error.\n", i);
		}
	}
	
	/* Compute Shadings	*/
	ComputeShading();
	///////////////////////////////////////////////
	// Initialize global data structures
	///////////////////////////////////////////////
	// Origin in Cameraview
	float P0[3] = { 0.0, };
	// Initialize ray vector
	float V[3] = { 0.0, };
	// Get T and invT matrix @ T[3], invT[3]
	float T[4][4] = { {1.0, 0.0, 0.0, -VRP[0]}, {0.0, 1.0, 0.0, -VRP[1]}, {0.0, 0.0, 1.0, -VRP[2]}, {0.0, 0.0, 0.0, 1.0} };
	float invT[4][4] = { {1.0, 0.0, 0.0, VRP[0]}, {0.0, 1.0, 0.0, VRP[1]}, {0.0, 0.0, 1.0, VRP[2]}, {0.0, 0.0, 0.0, 1.0} };
	// Get R and inverted R matrices from VUP and VPN
	float R[4][4];
	float invR[4][4];
	uvn(VUP, VPN, R, invR);
	// Modify Mwc and Mcw from R and T
	matrixcal(R, T, Mwc);
	matrixcal(invT, invR, Mcw);

	/* ================================================= */
	// The Main Ray-Tracing Volume Rendering Part.
	/* ================================================= */
	for (int i = 0; i < IMG_ROWS; i++) {
		for (int j = 0; j < IMG_COLS; j++) {
			int n = 0;
			// Get ray vectors for each pixel
			RayConstruction((float)i, (float)j, P0, V);
			// Save 2 intersection points P0 + t0*V & P0 + V1*V
			float ts[2] = { 0.0, };
			// Check intersectino points between ray and cube
			n=RayCubeIntersection(P0, V, ts);	// n is # of intersections
//			cout << ts[0] << " " << ts[1] << endl;
			if (n>1) {							// only if there are 2 intersections
				out_img[i][j] = VolumeRayTracing(P0, V, ts);
			}
		}
	}
	/* Save the output image */
	if (fopen_s(&fout, "./outimg.raw", "wb") != 0) {
		printf("Can't creat output image.\n");
	}
	printf(" ... Save the output image\n");
	n = fwrite(out_img, sizeof(char), IMG_ROWS * IMG_COLS, fout);
	if (n < IMG_ROWS * IMG_COLS * sizeof(char)) {
		printf("Write output image error.\n");
	}
	// Ready to write RAW image file
	fclose(fin);
	fclose(fout);
	return 0;
}