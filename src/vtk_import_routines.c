#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _MSC_VER
    #define FSCANF fscanf_s
#else
    #define FSCANF fscanf
#endif

extern int whiteSpaceCounter(FILE*);

const char* cellString = "CELLS"; 

void readVTK(const char* filePath, double* points, int* numPoints, int* cells, int* numCells){
	int numLine = 0;
	int index = 0;
	int i;
	int count = 0;

	char buffer[100];
	
	FILE* file = NULL;

#ifdef _MSC_VER
    fopen_s(&file,filePath, "r");
#else
    file = fopen(filePath, "r");
#endif	

	// Trash first 4 lines on file
	fgets(buffer, 100, file);
	fgets(buffer, 100, file);
	fgets(buffer, 100, file);
	fgets(buffer, 100, file);

	FSCANF(file, "%*s %d %*s\n", numPoints);

	printf("Number of Points = %d\n", (*numPoints));

	points = (double *) malloc(sizeof(double)*(*numPoints)*3);

	while (count<(*numPoints)*3) {
		numLine = whiteSpaceCounter(file);
		for (i=0; i<numLine; i++) {
			FSCANF(file, "%lf ", &points[count]);
			count++;
		}
		FSCANF(file, "%lf\n", &points[count]);
		count++;
	}

#ifdef _MSC_VER
	FSCANF(file, "%s %d %*s\n", buffer, 6, &numCells);
#else
      FSCANF(file, "%s %d %*s\n", buffer, numCells);
#endif	

	printf("Number of Cells = %d\n", (*numCells));
	cells = (int *) malloc(sizeof(int)*(*numCells)*4);
	
	if (strncmp(buffer, cellString, 6)!=0) {
        fprintf(stderr, "Incorrectly formatted VTK file -\nCells description in the wrong place\n");
        exit(EXIT_FAILURE);
    }  

	for (i=0; i<(*numCells); i++) {
		FSCANF(file, "%*d %d %d %d %d\n", &cells[i*4], &cells[i*4+1], &cells[i*4+2], &cells[i*4+3]);
	}

	fclose(file);
}

int whiteSpaceCounter(FILE* file) {
	int count = 0;
	int c;
	int wasSpace = 0;

	fpos_t pos;

	fgetpos (file,&pos);

	c = fgetc(file);

	// Ignore whitespace at start of line
	if (c==32) {
		while (c==32 && c!=10) {
			c = fgetc(file);
		}
	}

	while (1) {
		c = fgetc(file);
		if (c==10) {
			count -= ( wasSpace==0 ? 0 : 1);
			break;
		} else {
			if (c==32 && wasSpace==0) {
				count++;
				wasSpace=1;
			} else if (c!=32 && wasSpace==1) {
				wasSpace=0;
			}
		}
	}

	fsetpos (file,&pos);

	return count;
}



