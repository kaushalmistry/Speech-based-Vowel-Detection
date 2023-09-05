// Durbins_Algo.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include "stdio.h"
#include "conio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

#define p 12
#define PI 22/7
#define frame 320
char *vowels[] = {"a", "e", "i", "o", "u"};


// void calculateCis(long double A[][], long double C[]) {
// 	for (int m = 1; m <= p; m++) {
// 		C[m] = A[m][p];
// 		// long double sum = 0;
// 		for (int k = 1; k < m; k++) {
// 			// printf("\nDebug: %lf", A[m-k][p]);
// 			C[m] += ((double)k/m) * C[k] * A[m-k][p];
// 		}
// 		// C[m] += sum;
// 	}
// }


void durbinsAlgo(long double R[], long double C[]) {
	long double E[p+1];
	long double K[p+1];
	long double A[p+1][p+1];

	E[0] = R[0];
	// for (int i = 0; i <= p; i++) {
	// 	printf("\nR[%d] = %f", i, R[i]);
	// }

	for (int i = 1; i <= p; i++) {
		// printf("\n I = %d", i);
		long double sum = 0;
		for (int j = 1; j <= i-1; j++) {
			// printf("\nInner loop J = %d", j);
			sum += A[j][i-1] * R[i-j];
		}

		// printf("\n Sum = %f", sum);
		K[i] = (R[i] - sum) / E[i-1];
		// printf("\nK[%d] = %f", i, K[i]);

		A[i][i] = K[i];

		for (int j = 1; j <= i-1; j++) {
			// printf("\nInner loop 2 J = %d", j);
			A[j][i] = A[j][i-1] - (K[i] * A[i-j][i-1]);
			// printf("\nA[%d][%d] =  %f  ", j, i, A[j][i]);
		}

		E[i] = (1-(K[i]* K[i])) * E[i-1];
		// printf("\nE[%d] = %f", i, E[i]);
	}


	// printf("\nA[1] = %f", A[1][2]);
	// printf("\nA[2] = %f", A[2][2]);

	for (int i = 1; i <= p; i++) {
		// for (int j = 1; j <= p; j++) {
		// 	if (i<=j)
				// printf("\nA[%d] = %lf", i, A[i][p]);
		// }
	}


	// long double C[p+1];
	// long double wm;
	for (int m = 1; m <= p; m++) {
		C[m] = A[m][p];
		long double sum = 0;
		for (int k = 1; k < m; k++) {
			// printf("\nDebug: %lf", A[m-k][p]);
			sum += ((double)k/m) * C[k] * A[m-k][p];
		}
		// printf("\nSined Window: %lf", wm);
		// printf("\n Sum = %lf", sum);

		// sum += A[m][p];
		// C[m] = sum * wm;
		C[m] += sum;
		// temp += sum;

		// printf("\nC[%d] 1 = %lf", m, C[m]);
		// sum = C[m];
		// long double temp = C[m];
		// printf("\nTemp = %lf", temp);

		// C[m] = temp * wm;
		// C[m] = C[m] * wm;

		// printf("\nC[%d] = %lf", m, C[m]);
	}

	// calculateCis(A, C);

	printf("\n\nCi's are as below:");
	for (int i = 1; i <= p; i++) {
		printf("\nC[%d] = %lf", i, C[i]);
	}
}

void applyRaisedSineWindow(long double C[]) {
	long double wm;
	for (int i = 1; i <= p; i++) {
		wm = 1 + (p / 2 * sin((long double)PI * i / p)); // Calculating raised sine window
		C[i] = C[i] * wm;
	}
}

void applyDCShift(long double samples[], int totalSamples) {
	long double mean = 0;
	for (int i = 0; i < 400; i++) {
		mean += samples[i];
	}

	mean = mean / 400;

	for (int i = 0; i < totalSamples; i++) {
		samples[i] = samples[i] - mean;
	}
}


long double findNormalizationValue(long double samples[], int totalSamples)
{

    long double maxAmp = -1;

	for (int i = 0; i < totalSamples; i++) {
        maxAmp = abs(samples[i]) > abs(maxAmp) ? samples[i] : maxAmp;
    }

    return maxAmp;
}

void normalizeValuesWithHammingWindow(long double samples[], int totalSamples, long double maxAmplitude) {
	for(int i = 0; i < totalSamples; i++)
    {
		long double wn = 0.54 - 0.46 * cos((2 * (long double)PI * i) / (frame - 1));
		// printf("\n\nwn[%d] = %lf\n", i, wn);
		long double n = samples[i];
		// printf("Before normalization sample[%d] = %lf\n", i, n);
		// samples[i] =  (n * 5000 / maxAmplitude) * wn;
		samples[i] =  (n * 5000 / maxAmplitude);
		
		// printf("After normalization sample[%d] = %lf\n", i, samples[i]);
    }
}

void calculateRis(long double samples[], int start, long double R[]) {

	for (int i = 0; i <= p; i++) {
		long double sum = 0;
		int count = 0;
		for (int j = start; j < start + frame - i; j++) {
			count++;
			sum += samples[j] * samples[j+i];
		}
		R[i] = sum;
		printf("\n R[%d] = %lf with count = %d", i, R[i], count);
	}

}

// Finds the middle frame with the highest energy and return the start point for the 1st frame
int findStartPointForFrame(long double samples[], int totalSamples) {
	int startPoint = 0;

	long double maxEnergy = -1;

	long double curEnergy = 0;

	int count = 1;
	for (int i = 0; i < totalSamples; i++) {
		if (i % 320 == 0) {
			if (curEnergy > maxEnergy) {
				startPoint = i - 320;
				maxEnergy = curEnergy;
			}
			printf("%lf\n", curEnergy);
			curEnergy = 0;
			// count = 1;
		}
		// curEnergy += (samples[i] * samples[i]);
		// count++;
	}
	// printf("Start point = %d\n", startPoint);
	return startPoint - 640;

}


void printRecongnizedVowel(long double distance[]) {
	long double min = distance[0];
	int index = 0;

	for (int i = 0; i < 5; i++) {
		printf("\nThe tDistance [%d] = %lf", i, distance[i]);
		if (min > distance[i]) {
			min = distance[i];
			index = i;
		}
	}

	switch (index)
	{
	case 0:
		printf("\nVowel recognized as A");
		break;
		
	case 1:
		printf("\nVowel recognized as E");
		break;

	case 2:
		printf("\nVowel recognized as I");
		break;

	case 3:
		printf("\nVowel recognized as O");
		break;
	
	
	case 4:
		printf("\nVowel recognized as U");
		break;

	default:
		break;
	}
	printf("\n\n\n");
}



void TestInput(char *testFile, char *codeBook) {
	long double tConstant[] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};

	long double samples[10000]; // Initializing array to store the samples

	FILE* fptr = fopen(testFile, "r"); // Opening the input file to read
	int totalSamples = 0; // Total samples

	long double n; // Reading the input

	fscanf_s(fptr, "%lf", &n);
    while (!feof(fptr))
    {
		samples[totalSamples] = n;
		totalSamples++;
        fscanf_s(fptr, "%lf", &n);
    }
    fclose(fptr);

	long double codeBookArr[5][12];

	fptr = fopen(codeBook, "r");
	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 12; j++) {
        	fscanf_s(fptr, "%lf", &n);
			codeBookArr[i][j] = n;
		}
	}

	applyDCShift(samples, totalSamples);

	long double maxAmplitude = findNormalizationValue(samples, totalSamples);

	normalizeValuesWithHammingWindow(samples, totalSamples, maxAmplitude);


	int startPoint = findStartPointForFrame(samples, totalSamples);


	long double tDistance[5];
	for (int frameNo = 0; frameNo < 5; frameNo++) {
		long double R[p+1];
		calculateRis(samples, startPoint + frameNo*frame, R);

		long double C[p+1];

		durbinsAlgo(R, C);
		applyRaisedSineWindow(C);

		// Calculate Tokhura's distance
		
		for (int i = 0; i < 5; i++) {
			long double sum = 0;
			for (int j = 1; j <= p; j++) {
				sum += tConstant[j-1] * (C[j] - codeBookArr[i][j-1])*(C[j] - codeBookArr[i][j-1]);
			}
			if (frameNo > 0) {
				tDistance[i] = sum < tDistance[i] ? sum : tDistance[i];
			} else {
				tDistance[i] = sum;
			}
		}
	}
	printRecongnizedVowel(tDistance);

}


/**
 * @brief It generates the CI's and stores in the output file
 */
void createCiFrames(char* fileName, char* outputFile) {
	long double samples[30000]; // Initializing array to store the samples

	FILE* fptr = fopen(fileName, "r"); // Opening the input file to read
	int totalSamples = 0; // Total samples

	long double n; // Reading the input

	fscanf_s(fptr, "%lf", &n);
    while (!feof(fptr))
    {
		samples[totalSamples] = n;
		totalSamples++;
        fscanf_s(fptr, "%lf", &n);
    }
    fclose(fptr);

	//applyDCShift(samples, totalSamples);

	long double maxAmplitude = findNormalizationValue(samples, totalSamples);

	// normalizeValuesWithHammingWindow(samples, totalSamples, maxAmplitude);


	int startPoint = findStartPointForFrame(samples, totalSamples);


	FILE* outputFilePtr = fopen(outputFile, "a");

	for (int frameNo = 0; frameNo < 5; frameNo++) {
		long double R[p+1];
		calculateRis(samples, startPoint + frameNo*frame, R);

		long double C[p+1];

		durbinsAlgo(R, C);
		applyRaisedSineWindow(C);

		for (int i = 1; i <= p; i++) {
			// printf("\nC[%d] = %lf", i, C[i]);
            fprintf(outputFilePtr, "%lf ", C[i]);
		}
		fprintf(outputFilePtr, "\n");

	}
	
	fclose(outputFilePtr);

}


void createCodeBook() {
	char ipFile[100];
	char outputFileName[100];
	char fileNo[5];


	for (int i = 0; i < 1; i++) {
		strcpy(outputFileName, "files/");
		strcat(outputFileName, vowels[i]);
		strcat(outputFileName, "_50_Frames.txt");

		// printf("\nOutput file Name = %s", outputFileName);

		for (int j = 1; j <= 1; j++) {	
			strcpy(ipFile, "files/vowels_files/");
			strcat(ipFile, vowels[i]);
			strcat(ipFile, "_Files/224101031_");
			strcat(ipFile, vowels[i]);
			sprintf(fileNo, "%d", j);
			strcat(ipFile, "_");
			strcat(ipFile, fileNo);
			strcat(ipFile, ".txt");
			printf("\n\n\nInput file Name = %s", ipFile);
			createCiFrames(ipFile, outputFileName);
		}
	}


	long double arr[5][12] = {0};
	int row = 0;
	long double C;

	for (int i = 0; i < 5; i++) {
		strcpy(ipFile, "files/");
		strcat(ipFile, vowels[i]);
		strcat(ipFile, "_50_Frames.txt");

		FILE* fptr = fopen(ipFile, "r");
		for (int i = 0; i < 50; i++) {
			for (int j = 0; j < 12; j++) {
				fscanf(fptr, "%lf", &C);
				arr[row][j] += C;
			}
		}
		fclose(fptr);
		row++;
	}

	char *finalFile = "files/vowelCodeBook.txt";
	FILE* codeBookPtr = fopen(finalFile, "w");

	for (int i = 0; i < 5; i++) {
		for (int j = 0; j < 12; j++) {
			fprintf(codeBookPtr, "%lf ", arr[i][j] / 50);
		}
		fprintf(codeBookPtr, "\n");
	}
	
	fclose(codeBookPtr);

}

int _tmain(int argc, _TCHAR* argv[])
{

	createCodeBook();

	// int i = 4;
	// char ipFile[100];
	// char fileNo[5];
	// char *codeBook = "files/vowelCodeBook.txt";

	// for (int j = 1; j <= 10; j++) {
	// 	// strcpy(ipFile, "files/vowels_files/");
	// 	// strcat(ipFile, vowels[i]);
	// 	// strcat(ipFile, "_Files/224101031_");
	// 	// strcat(ipFile, vowels[i]);
	// 	// sprintf(fileNo, "%d", j);
	// 	// strcat(ipFile, "_");
	// 	// strcat(ipFile, fileNo);
	// 	// strcat(ipFile, ".txt");
		
	// 		strcpy(ipFile, "files/vowel/");
	// 		strcat(ipFile, vowels[i]);
	// 		strcat(ipFile, "_test/");
	// 		strcat(ipFile, vowels[i]);
	// 		sprintf(fileNo, "%d", j);
	// 		strcat(ipFile, fileNo);
	// 		strcat(ipFile, ".txt");
	// 	printf("\nInput file Name = %s", ipFile);
	// 	TestInput(ipFile, codeBook);
	// }
	
	getchar();
	return 0;
}