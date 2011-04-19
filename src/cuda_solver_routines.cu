#include <stdlib.h>
#include <stdio.h>

#include <cutil_inline.h>

//extern __shared__ double shared_array[];
//const int algebraicCount = 25;
//const int rateStateCount = 8;
//const int constantCount = 0;
//const double FLOPSPerFunction = 193.0f;
//const int DEFAULT_TESTING_THREADS = 2000000;
//const int sharedMemoryCellModel = 0;
//const char* cellModelName = "LR R3";
//
//////////////////////////////////////////////////////////////////////////////////
//// Cell Model Device Functions
//////////////////////////////////////////////////////////////////////////////////
//__device__ void computeRates(float VOI, double* DUMMY, double* STATES, double* ALGEBRAIC, double* RATES)
//{
//	ALGEBRAIC[0] = -25.5; // Add stimulus in proper
//	ALGEBRAIC[1] = (0.32f*STATES[0]+15.0816f)/(1.0f - (expf(- 0.1f*STATES[0]-4.713f))); // 7 ops
//
//	if (STATES[0] < -40.0f) {
//		ALGEBRAIC[2] = 0.135f*(expf(((80.0f+STATES[0])/- 6.8f))); // 4 ops
//		ALGEBRAIC[3] = (( - 127140*(expf(0.2444f*STATES[0])) - 3.47400e-05*(expf(-0.04391f*STATES[0])))*(STATES[0]+37.78f))/(1.0f+(expf(0.311f*STATES[0]+24.64053))); // 14 ops
//		ALGEBRAIC[9] = 3.56f*(expf(0.079f*STATES[0]))+ 310000*(expf(0.35f*STATES[0]));  // 7 ops
//		ALGEBRAIC[10] = 0.1212f*(expf(-0.01052f*STATES[0]))/(1.0f+(expf(-0.1378f*STATES[0]-5.531292f))); // 8 ops
//	} else {
//		ALGEBRAIC[2] = 0.00000;
//		ALGEBRAIC[3] = 0.00000;
//		ALGEBRAIC[9] = 1.0f/( 0.13f*(1.0f+(expf(((STATES[0]+10.66f)/- 11.1f)))));
//		ALGEBRAIC[10] = ( 0.3f*(expf(-2.53500e-07*STATES[0])))/(1.0f+(expf(-0.1f*STATES[0]-3.2f)));
//	}
//	if (STATES[0] < -100.0f) {
//		ALGEBRAIC[16] = 2.837f*(expf(0.04f*STATES[0]+3.08f) - 1.0f)/((STATES[0]+77.0f)*(expf(0.04f*STATES[0]+1.4f))); // 11 ops
//	} else {
//		ALGEBRAIC[16] = 1.0f;
//	}
//
//	ALGEBRAIC[4] = (0.095f*(expf(-0.01f*STATES[0] + 0.5f)))/(1.0f+(expf(-0.072*STATES[0] + 0.36f))); // 9 ops
//	ALGEBRAIC[5] = (0.012f*(expf(-0.008f*STATES[0]-0.224f)))/(1.0f+(expf(0.15f*STATES[0]+4.2f))); // 9 ops
//	ALGEBRAIC[6] = (0.0005f*(expf(0.083f*STATES[0]+4.15f)))/(1.0f+(expf(0.057f*STATES[0]+2.85f))); // 9 ops
//	ALGEBRAIC[7] =  23*(powf(STATES[1], 3.0f))*STATES[2]*STATES[3]*(STATES[0] - 54.794463f); // 6 ops
//	ALGEBRAIC[8] =  0.08f*(expf(-STATES[0]/11.0000)); // 3 ops
//	ALGEBRAIC[11] = (0.07f*(expf(-0.017f*STATES[0]-0.748f)))/(1.0f+(expf(0.05f*STATES[0]+2.2f))); // 9 ops
//	ALGEBRAIC[12] = (0.0065f*(expf(-0.02f*STATES[0]-0.6f)))/(1.0f+(expf(-0.2f*STATES[0]-6.0f))); // 9 ops
//	ALGEBRAIC[13] = (0.0013f*(expf(-0.06f*STATES[0]-1.2f)))/(1.0f+(expf(-0.04f*STATES[0]-0.8f))); // 9 ops
//	ALGEBRAIC[14] = 7.7f - 13.0287f*logf(STATES[4]); // 3 ops
//	ALGEBRAIC[15] =  0.09f*STATES[5]*STATES[6]*(STATES[0] - ALGEBRAIC[14]); // 4 ops
//	ALGEBRAIC[17] =  0.282f*STATES[7]*ALGEBRAIC[16]*(STATES[0] + 77.56758f); // 4 ops
//	ALGEBRAIC[18] = 1.02f/(1.0f+(expf(0.2385f*STATES[0] + 6.83967915f))); // 4 ops
//	ALGEBRAIC[19] = (0.49124f*(expf( 0.08032f *STATES[0] + 7.49939f) + expf(0.06175f*STATES[0] - 31.271255925f)))/(1.00000+expf(-0.514300*STATES[0] - 214.85137268791f)); // 13 ops
//	ALGEBRAIC[20] = ALGEBRAIC[18]/(ALGEBRAIC[18]+ALGEBRAIC[19]); // 2 ops
//	ALGEBRAIC[21] =  0.6047f*ALGEBRAIC[20]*(STATES[0] + 87.89290f); // 3 ops
//	ALGEBRAIC[22] = 1.0f/(1.0f+(exp(((7.488f - STATES[0])/5.98f)))); // 5 ops
//	ALGEBRAIC[23] =  0.0183f*ALGEBRAIC[22]*(STATES[0] + 87.89290f); // 3 ops
//	ALGEBRAIC[24] =  0.03921f*STATES[0] +2.3475027f; // 3 ops
//
//	RATES[0] =  -(ALGEBRAIC[0]+ALGEBRAIC[7]+ALGEBRAIC[15]+ALGEBRAIC[17]+ALGEBRAIC[21]+ALGEBRAIC[23]+ALGEBRAIC[24]); // 7 ops
//	RATES[1] =  ALGEBRAIC[1]*(1.00000 - STATES[1]) -  ALGEBRAIC[8]*STATES[1]; // 4 ops
//	RATES[2] =  ALGEBRAIC[2]*(1.00000 - STATES[2]) -  ALGEBRAIC[9]*STATES[2]; // 4 ops
//	RATES[3] =  ALGEBRAIC[3]*(1.00000 - STATES[3]) -  ALGEBRAIC[10]*STATES[3]; // 4 ops
//	RATES[4] =  - 0.0001f*ALGEBRAIC[15]+ 0.000007f - 0.07f*STATES[4]; // 4 ops
//	RATES[5] =  ALGEBRAIC[4]*(1.00000 - STATES[5]) -  ALGEBRAIC[11]*STATES[5]; // 4 ops
//	RATES[6] =  ALGEBRAIC[5]*(1.00000 - STATES[6]) -  ALGEBRAIC[12]*STATES[6]; // 4 ops
//	RATES[7] =  ALGEBRAIC[6]*(1.00000 - STATES[7]) -  ALGEBRAIC[13]*STATES[7]; // 4 ops
//}

extern __shared__ double shared_array[];
const int rateStateCount = 2;
const int constantCount = 0;
const int algebraicCount = 1;
const double FLOPSPerFunction = 9.0f;
const int DEFAULT_TESTING_THREADS = 20000000;
const int sharedMemoryCellModel = 0;
const char* cellModelName = "FN R1";

////////////////////////////////////////////////////////////////////////////////
// Cell Model Device Functions
////////////////////////////////////////////////////////////////////////////////
__device__ void computeRates(double time, double* constants, double* states, double* algebraic, double* rates)
{
	rates[1] =  0.005f*(states[0] - 3.0f*states[1]);
	rates[0] =  ((states[0]*(states[0] - -0.08f)*(1.0f - states[0]) - states[1])+ algebraic[0]);
}

void initProblem(int num_threads, double* STATES) {
	printf("test\n\n\n\n\n\n");
}

//////////////////////////////////////////////////////////////////////////////////
//// Cell Model Host Functions ////////// Should Not Be Needed Later /////////////
//////////////////////////////////////////////////////////////////////////////////
//void initProblem(int num_threads, double* STATES)
//{
//	int i;
//
//	STATES[0] = -84.3801107371;
//	STATES[1] = 0.00171338077730188;
//	STATES[2] = 0.982660523699656;
//	STATES[3] = 0.989108212766685;
//	STATES[4] = 0.00017948816388306;
//	STATES[5] = 0.00302126301779861;
//	STATES[6] = 0.999967936476325;
//	STATES[7] = 0.0417603108167287;
//
//	for (i=1; i<num_threads; i++)
//		memcpy(STATES + i*rateStateCount, STATES, rateStateCount*sizeof(double));
//}

const int DEFAULT_TESTING_TIMESTEPS = 1000;
const double FLOPSPerTimeStep = 2.0f;
const int FunctionEvals = 1;
const char* integratorName = "E R1";
const int sharedMemoryIntegrator = 0;

////////////////////////////////////////////////////////////////////////////////
// ODE Integrator Device Functions
////////////////////////////////////////////////////////////////////////////////
__device__ void integrator(int timeSteps, float stepSize, double* constants, double* states, double* algebraic)
{
	int i,j;

	double rates[rateStateCount];

#pragma unroll 40
	for (i=1; i<timeSteps+1; i++) {
		computeRates(i*stepSize, constants, states, algebraic, rates);

		for (j=0; j<rateStateCount; j++) {
			states[j] += stepSize*rates[j];
		}
	}
}

//const int DEFAULT_TESTING_TIMESTEPS = 1000;
//const double FLOPSPerTimeStep = 22.0f;
//const int FunctionEvals = 4;
//const char* integratorName = "RK R2";
//const int sharedMemoryIntegrator = 0;
//
//////////////////////////////////////////////////////////////////////////////////
//// ODE Integrator Device Functions
//////////////////////////////////////////////////////////////////////////////////
//__device__ void integrator(int timeSteps, float stepSize, double* constants, double* states, double* algebraic)
//{
//	int i,j;
//
//	double previousStates[rateStateCount];
//	double kutta[rateStateCount];
//	double offsets[rateStateCount];
//
//	for (j=0; j<rateStateCount; j++) {
//		previousStates[j] = states[j];
//	}
//
//
//#pragma unroll 40
//	for (i=1; i<timeSteps+1; i++) {
//		computeRates(i*stepSize, constants, previousStates, algebraic, kutta);
//
//		for (j=0; j<rateStateCount; j++) {
//			offsets[j] = previousStates[j] + 0.5f*kutta[j]*stepSize; // 3 ops
//			states[j] += (stepSize/6)*(kutta[j]); // 3 ops
//		}
//		computeRates(i*stepSize, constants, offsets, algebraic, kutta);
//
//		for (j=0; j<rateStateCount; j++) {
//			offsets[j] = previousStates[j] + 0.5f*kutta[j]*stepSize; // 3 ops
//			states[j] += (stepSize/6)*(2*kutta[j]); // 4 ops
//		}
//		computeRates(i*stepSize, constants, offsets, algebraic, kutta);
//
//		for (j=0; j<rateStateCount; j++) {
//			offsets[j] = previousStates[j] + kutta[j]*stepSize; // 2 ops
//			states[j] += (stepSize/6)*(2*kutta[j]); // 4 ops
//		}
//		computeRates(i*stepSize, constants, offsets, algebraic, kutta);
//
//		for (j=0; j<rateStateCount; j++) {
//			states[j] += (stepSize/6)*(kutta[j]); // 3 ops
//			previousStates[j] = states[j];
//		}
//	}
//}

const int sharedMemoryDevice = rateStateCount + algebraicCount;

////////////////////////////////////////////////////////////////////////////////
// Solver Kernel
////////////////////////////////////////////////////////////////////////////////
__global__ void solveSystem(int timeSteps, float stepSize, double* states)
{
	int i;

	double* threadAlgebraic = &shared_array[algebraicCount*threadIdx.x + sharedMemoryIntegrator*blockDim.x];
	double* threadStates = &shared_array[(algebraicCount + sharedMemoryIntegrator)*blockDim.x];

#if constantCount > 0
	double threadConstants[constantCount];

	intitialiseConstants(threadConstants);
#endif

	for (i=0; i<rateStateCount; i++) {
		threadStates[threadIdx.x + i*blockDim.x] = states[threadIdx.x + i*blockDim.x + blockIdx.x*blockDim.x*rateStateCount];
	}

#if constantCount > 0
	integrator(timeSteps, stepSize, threadConstants, threadStates + threadIdx.x*rateStateCount, threadAlgebraic);
#else
	integrator(timeSteps, stepSize, NULL, threadStates + threadIdx.x*rateStateCount, threadAlgebraic);
#endif

	for (i=0; i<rateStateCount; i++) {
		states[threadIdx.x + i*blockDim.x + blockIdx.x*blockDim.x*rateStateCount] = threadStates[threadIdx.x + i*blockDim.x];
	}
}

////////////////////////////////////////////////////////////////////////////////
// Check Error Checking
////////////////////////////////////////////////////////////////////////////////
void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg,
                                  cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }                         


}

void *safeMalloc(size_t size)  // Only called once - Pointless??
{
	void *ret = malloc(size);
	if (ret == NULL) {
		fprintf(stderr, "malloc of %zu bytes failed", size);
		exit(EXIT_FAILURE);
	}
	return ret;
}

void cleanup() { cudaThreadExit(); }

void domainDecomposition(unsigned int* threads_per_domain, unsigned int* spill, unsigned int num_threads, int unsigned num_streams, size_t* pagedMemorySize, unsigned int* threads_per_block,
		unsigned int* num_partitions, unsigned int* blocksPerDomain, unsigned int* num_blocks, size_t* sharedMem){

	size_t freeGPUMemory = 0;
	size_t test = 0;
	int cuda_device = 0;
	unsigned int remainingBlocks = 0;
	unsigned int numberRemainingDomains = 0;
	unsigned int count;

    cudaDeviceProp deviceProp;
    cutilSafeCall( cudaGetDeviceProperties(&deviceProp, cuda_device) );

    (*threads_per_block) = (*threads_per_block) > deviceProp.maxThreadsPerBlock ? deviceProp.maxThreadsPerBlock : (*threads_per_block);

	// Adjust threads per blocks so that states variables fit in shared memory
	count=0;
	(*sharedMem) = (size_t) (sharedMemoryIntegrator + sharedMemoryDevice + sharedMemoryCellModel) * (*threads_per_block) * sizeof(double);
	while ( (*sharedMem) > deviceProp.sharedMemPerBlock*0.5 && count < 200 ) {
		if ((*threads_per_block) > 32 ) {
			(*threads_per_block)-=32;
		} else {
			fprintf(stderr, "Cannot fit variables in shared memory");
			exit(EXIT_FAILURE);
		}
		count++;
		(*sharedMem) = (size_t) (sharedMemoryIntegrator + sharedMemoryDevice + sharedMemoryCellModel) * (*threads_per_block) * sizeof(double);
	}
	(*num_blocks)=num_threads/(*threads_per_block);
	spill[2]=num_threads % (*threads_per_block);
	//if (spill[2] !=0) (*num_blocks)++; // Round up calculation

	(*blocksPerDomain) = (*num_blocks)/(num_streams*(*num_partitions));
	if ((*blocksPerDomain) == 0) {
		fprintf(stderr, "Too many streams and partitions for your problem size\nReduce these and try again");
		exit(EXIT_FAILURE);
	} else if ((*blocksPerDomain) > deviceProp.maxGridSize[0]) {
		fprintf(stderr, "Too many blocks allocated to each domain for your problem size\nIncrease number of partitions and try again");
		exit(EXIT_FAILURE);
	}
	remainingBlocks = (*num_blocks)%(num_streams*(*num_partitions));
	if (remainingBlocks > 0 && remainingBlocks <= num_streams*(*blocksPerDomain)) {
		(*num_partitions)++;
	} else if (remainingBlocks > num_streams*(*blocksPerDomain)) {
		numberRemainingDomains = ceil((float)remainingBlocks/(*blocksPerDomain));
		(*num_partitions) += ceil((float)numberRemainingDomains/num_streams);
	}

	if ((*num_blocks)/(num_streams*(*num_partitions))==0) {
		(*num_partitions)--;
	}

	remainingBlocks = (*num_blocks)%(*num_partitions);

	// Adjust number of partitions so that data fits in GPU memory
	count = 0;
	cutilSafeCall( cudaMemGetInfo(&freeGPUMemory, &test) );
	*pagedMemorySize = sizeof(double)*rateStateCount*(*blocksPerDomain)*num_streams*(*threads_per_block);
	while ( *pagedMemorySize > freeGPUMemory*0.9 && count < 200 ) {
		if ((*blocksPerDomain) > 1 ) {
			(*num_partitions) = (*num_blocks)/num_streams/--(*blocksPerDomain);
			remainingBlocks = (*num_blocks)%(num_streams*(*num_partitions));
			if (remainingBlocks != 0) {
				(*num_partitions)++;
			} else if (remainingBlocks > num_streams*(*blocksPerDomain)) {
				numberRemainingDomains = ceil((float)remainingBlocks/(*blocksPerDomain));
				(*num_partitions) += ceil((float)numberRemainingDomains/num_streams);
			}
			remainingBlocks = (*num_blocks)%(*num_partitions);
		} else {
			fprintf(stderr, "Cannot fit variables in GPU device memory\nReduce threads per block and try again");
			exit(EXIT_FAILURE);
		}
		count++;

		*pagedMemorySize = sizeof(double)*rateStateCount*(*blocksPerDomain)*num_streams*(*threads_per_block);
	}

	if (*pagedMemorySize > freeGPUMemory*0.9) {
		fprintf(stderr, "Cannot fit variables in GPU device memory\nReduce threads per block and try again");
		exit(EXIT_FAILURE);
	}

	(*threads_per_domain) = (*threads_per_block)*(*blocksPerDomain);
	spill[0] = remainingBlocks/(*blocksPerDomain);
	spill[1] = remainingBlocks%(*blocksPerDomain);

//	if (renainingBlocks != 0) {
//		blocksLastStream =
//	}
//	remainingBlocks = num_threads%(num_blocks*(*num_partitions))
//	for (i=0; i<(*num_partitions); i++) {
//		if (i == (*num_partitions) && remainingBlocks != 0) {
//
//			for (j=0; j<remainingStreams; j++) {
//				domainIndex = j + i*num_streams;
//				(*threads_per_domain) = (*threads_per_block)*(*blocksPerDomain);
//				count += threads_per_domain[domainIndex];
//			}
//			if (blocksLastStream!=0) {
//				domainIndex ++;
//				threads_per_domain[domainIndex] = num_threads - count;
//			}
//		} else {
//			for (j=0; j<num_streams; j++) {
//				domainIndex = j + i*num_streams;
//				threads_per_domain[domainIndex] = (*threads_per_block)*(*blocksPerDomain);
//				count += threads_per_domain[domainIndex];
//			}
//		}
//	}
}

////////////////////////////////////////////////////////////////////////////////
// CPU side solution routine
////////////////////////////////////////////////////////////////////////////////
void solve(double* h_states, double startTime, double endTime, double stepSize,
		   unsigned int num_threads, unsigned int threads_per_block, unsigned int num_partitions, 
		   int unsigned num_streams, FILE *timing_file)
{
    unsigned int i,j;
	unsigned int num_blocks = 0;
	unsigned int local_offset = 0;
	unsigned int global_offset = 0;
	double nFLOPS = 0;
	double dSeconds = 0;
	double gflops = 0;
	double kernel_dSeconds = 0;
	double kernel_nFLOPS = 0;
	double kernel_gflops = 0;
	unsigned int timer = 0;
	unsigned int kernel_timer = 0;
	unsigned int domainIndex = 0;
	unsigned int blocksPerDomain = 0;
	unsigned int lastFullDomain = 0;
	unsigned int threads_per_domain = 0;
	unsigned int spill[3] = { 0 };

	size_t pagedMemorySize;
	size_t sharedMem = 0;
	size_t lastStreamMemorySize = 0;

	double *d_states = 0;
	double *h_paged_states = 0;

	cudaStream_t *streams;
	int cuda_device = 0;

	unsigned int timeSteps = (endTime-startTime)/stepSize;

    // Check for a CUDA compatible device
    int num_devices=0;
    cutilSafeCall( cudaGetDeviceCount(&num_devices) );
    if(0==num_devices)
    {
        printf("Your system does not have a CUDA capable device\n");
        return;
	}
	
    // Set appropriate device flags and print relevant information
	cutilSafeCall( cudaSetDevice( cuda_device ) );
	cutilSafeCall( cudaSetDeviceFlags(cudaDeviceBlockingSync) );
    cudaDeviceProp deviceProp;
    cutilSafeCall( cudaGetDeviceProperties(&deviceProp, cuda_device) );
    if( (1 == deviceProp.major) && (deviceProp.minor < 1))
        printf("%s does not have compute capability 1.1 or later\n", deviceProp.name);

    domainDecomposition(&threads_per_domain, spill, num_threads, num_streams, &pagedMemorySize, &threads_per_block,
    		&num_partitions, &blocksPerDomain, &num_blocks, &sharedMem);


    lastFullDomain = num_blocks/blocksPerDomain;

	// Allocate CUDA device and host pinned memory
	cutilSafeCall( cudaMallocHost((void**) &h_paged_states, pagedMemorySize) );
	cutilSafeCall( cudaMalloc((void **) &d_states, pagedMemorySize) );

	// Create streams
	streams = (cudaStream_t*) malloc(num_streams * sizeof(cudaStream_t));
    for(i = 0; i < num_streams; i++)
        cutilSafeCall( cudaStreamCreate(&(streams[i])) );

    // Setup execution parameters
    dim3  grid(blocksPerDomain, 1, 1);
    dim3  threads(threads_per_block, 1, 1);

    //#ifdef DEBUG
    if (timing_file) {
    	printf("> Device name : %s\n", deviceProp.name );
		printf("> CUDA Capable SM %d.%d hardware with %d multi-processors\n",
			deviceProp.major, deviceProp.minor, deviceProp.multiProcessorCount);
		printf("> Cell Model = %s, Integrator = %s\n", cellModelName, integratorName);
		printf("> num_threads = %d, num_blocks = %d, threads_per_block = %d, num_partitions = %d, timeSteps = %d, num_streams = %d\n",
			num_threads, num_blocks, threads_per_block, num_partitions, timeSteps, num_streams);
		//#endif
    	printf("grid.x %d threads.x %d sharedMem %d\n", grid.x, threads.x, sharedMem);
    	printf("Spills %d %d %d\n", spill[0], spill[1], spill[2]);

		// Setup and start global timer
		timer = 0;
		cutCreateTimer(&timer);
		cutilCheckError(cutStartTimer(timer));

		// Test kernel speed in default stream (timing is more accurate in default stream)
		memcpy(h_paged_states, h_states, pagedMemorySize/num_streams);
		cutilSafeCall( cudaMemcpy(d_states, h_paged_states, pagedMemorySize/num_streams, 
			cudaMemcpyHostToDevice) );
		// Start kernel timer
		kernel_timer = 0;
		cutCreateTimer(&kernel_timer);
		cutilCheckError(cutStartTimer(kernel_timer));
		// Start kernel
		solveSystem<<<grid, threads, sharedMem>>>(timeSteps, stepSize, d_states);
		checkCUDAError("Single Kernel Execution");
		cutilSafeCall( cudaThreadSynchronize() );
		// Stop kernel Timer 
		cutilCheckError(cutStopTimer(kernel_timer));
		cutilSafeCall( cudaMemcpy(h_paged_states, d_states, pagedMemorySize/num_streams, 
			cudaMemcpyDeviceToHost) );
		memcpy(h_states, h_paged_states, pagedMemorySize/num_streams);

		// Prefetch data for next partition in first stream
		if (num_partitions>1) {
			global_offset = rateStateCount * num_streams * grid.x * threads.x;
			memcpy(h_paged_states, h_states + global_offset, pagedMemorySize/num_streams);
		}
	} else {
		memcpy(h_paged_states, h_states, pagedMemorySize);
	}

	// Queue kernel calls into streams to hide memory transfers (num_partitions sets of kernel calls in each stream)
	for(i = 0; i < num_partitions+1; i++) {
		// Asynchronously launch num_streams memcopies
		for(j = 0; j < num_streams; j++) {
			domainIndex = j + i*num_streams;
			if (domainIndex <= lastFullDomain && (timing_file==NULL || !(i==0 && j==0))) {
				local_offset = j * rateStateCount * grid.x * threads.x ;
				if (domainIndex == lastFullDomain && (spill[1]!=0 || spill[2]!=0)) {
					//printf("last async in %d, size %d\n", domainIndex, lastStreamMemorySize);
					cutilSafeCall( cudaMemcpyAsync(d_states + local_offset, h_paged_states + local_offset,
							lastStreamMemorySize, cudaMemcpyHostToDevice, streams[j]) );
				} else {
					//printf("normal async in %d, size %d\n", domainIndex, pagedMemorySize/num_streams);
					cutilSafeCall( cudaMemcpyAsync(d_states + local_offset, h_paged_states + local_offset,
							pagedMemorySize/num_streams, cudaMemcpyHostToDevice, streams[j]) );
				}
			}
		}
		// Execute the kernel
		// Asynchronously launch num_streams kernels, each operating on its own portion of data
		for(j = 0; j < num_streams; j++) {
			domainIndex = j + i*num_streams;
			if (domainIndex <= lastFullDomain && (timing_file==NULL || !(i==0 && j==0))) {
				local_offset = j * grid.x * threads.x ;
				if (domainIndex == lastFullDomain && (spill[1]!=0 || spill[2]!=0)) {
				    grid.x = spill[1]+1;
					solveSystem<<<grid, threads, sharedMem, streams[j]>>>(timeSteps, stepSize,
							d_states + rateStateCount*local_offset);
				} else {
					solveSystem<<<grid, threads, sharedMem, streams[j]>>>(timeSteps, stepSize,
							d_states + rateStateCount*local_offset);
				}
			}
		}

		// Asynchronoously launch num_streams memcopies
		for(j = 0; j < num_streams; j++){
			domainIndex = j + i*num_streams;
			if (domainIndex <= lastFullDomain && (timing_file==NULL || !(i==0 && j==0))) {
				local_offset = j * rateStateCount * grid.x * threads.x ;
				if (domainIndex == lastFullDomain && (spill[1]!=0 || spill[2]!=0)) {
					//printf("last async out %d, size %d\n", domainIndex, lastStreamMemorySize);
					cutilSafeCall( cudaMemcpyAsync(h_paged_states + local_offset, d_states + local_offset,
						lastStreamMemorySize, cudaMemcpyDeviceToHost, streams[j]) );
				} else {
					//printf("normal async out %d, size %d\n", domainIndex, pagedMemorySize/num_streams);
					cutilSafeCall( cudaMemcpyAsync(h_paged_states + local_offset, d_states + local_offset,
						pagedMemorySize/num_streams, cudaMemcpyDeviceToHost, streams[j]) );
				}
			}
		}
		
		// Execute memcpys in and out of paged memory when CUDA calls in the streams have finished
		for(j = 0; j < num_streams; j++){
			domainIndex = j + i*num_streams;
			if (domainIndex <= lastFullDomain && (timing_file==NULL || !(i==0 && j==0))) {
				cudaStreamSynchronize(streams[j]);

				local_offset = j * rateStateCount * grid.x * threads.x ;
				global_offset = i * num_streams * grid.x * threads.x;
				
				if (domainIndex == lastFullDomain && (spill[1]!=0 || spill[2]!=0)) {
					//printf("last memcpy out %d\n", domainIndex);
					memcpy(h_states + rateStateCount * global_offset + local_offset, h_paged_states + local_offset,
							lastStreamMemorySize);
				} else {
					//printf("normal memcpy out %d\n", domainIndex);
					memcpy(h_states + rateStateCount * global_offset + local_offset, h_paged_states + local_offset,
						pagedMemorySize/num_streams);
				}

				global_offset = (i + 1) * num_streams * grid.x * threads.x;
				if (domainIndex == lastFullDomain - num_streams && (spill[1]!=0 || spill[2]!=0)) {
					//printf("last memcpy in %d\n", domainIndex);
					lastStreamMemorySize = sizeof(double)*rateStateCount*(spill[1]*threads_per_block+spill[2]);
					memcpy(h_paged_states + local_offset, h_states + rateStateCount * global_offset + local_offset,
							lastStreamMemorySize);
				} else if (domainIndex < lastFullDomain - num_streams) {
					//printf("normal memcpy in %d\n", domainIndex);
					memcpy(h_paged_states + local_offset, h_states + rateStateCount * global_offset + local_offset,
						pagedMemorySize/num_streams);
				}
			}
		}
	}

	if (timing_file) {
		// Stop global timer
		cutilCheckError(cutStopTimer(timer));

		// Calculate timing statistics
		dSeconds = cutGetTimerValue(timer)/1000.0;
		nFLOPS = (FLOPSPerTimeStep*rateStateCount + FLOPSPerFunction*FunctionEvals)*timeSteps*num_threads;
		gflops = 1.0e-9 * nFLOPS/dSeconds;

		kernel_dSeconds = cutGetTimerValue(kernel_timer)/1000.0;
		kernel_nFLOPS = (FLOPSPerTimeStep*rateStateCount + FLOPSPerFunction*FunctionEvals)*timeSteps*num_threads/num_streams/num_partitions;
		kernel_gflops = 1.0e-9 * kernel_nFLOPS/kernel_dSeconds;

		// Store Stats
		fprintf(timing_file,"%s\t%s\t%d\t%d\t%d\t%d\t%d\t%f\t%f\t%f\t%f\t%f\n", cellModelName, integratorName, num_threads, num_blocks, 
			threads_per_block, num_partitions, num_streams, dSeconds, gflops, kernel_dSeconds, 
			kernel_gflops, gflops/kernel_gflops*100);
	}

	// Deallocate Memory and Release Threads
	for(i = 0; i < num_streams; i++)
        	cutilSafeCall( cudaStreamDestroy(streams[i]) );
	cutilSafeCall( cudaFree(d_states) );
	cutilSafeCall( cudaFreeHost(h_paged_states) );
	cudaThreadExit();
}


////////////////////////////////////////////////////////////////////////////////
// Auxilliary Testing Functions
////////////////////////////////////////////////////////////////////////////////

void solveProblem(unsigned int timeSteps, unsigned int num_threads, unsigned int threads_per_block, 
				  unsigned int num_partitions, unsigned int num_streams, FILE *timing_file, FILE *results_file)
{
	unsigned int i, j; 
	float startTime = 0.0f;
	float endTime = 0.2f;

	double* h_states = NULL;

	h_states = (double *) safeMalloc(sizeof(double)*rateStateCount*num_threads);

	initProblem(num_threads, h_states);

	solve(h_states, startTime, endTime, timeSteps, num_threads, threads_per_block, num_partitions, num_streams, timing_file);

	if (results_file) {
		fprintf(results_file,"\n\n");

		for (i=0; i<num_threads; i++) {
			fprintf(results_file,"%d", i+1);
			for (j=0; j<rateStateCount; j++) {
				fprintf(results_file,"\t%f", h_states[i*rateStateCount+j]);
			}
			fprintf(results_file,"\n");
		}
	}

	free(h_states);
}

void startToFinish (unsigned int num_threads,unsigned int threads_per_block,unsigned int num_partitions,unsigned int num_streams,unsigned int timeSteps)
{
	//int i,j,k;
 
	FILE *file = NULL; 
	FILE *file1 = NULL;
 
	file = fopen("performance_data.txt", "rt");
	//file1 = fopen("results_data.txt", "wt");

	if (!file) {
		file = fopen("performance_data.txt", "wt");
		fprintf(file,"Cell Model\tIntegrator\tNumber of Threads\tNumber 0f Blocks\tThreads Per Block\tNumber of Partitions\tNumber of Streams\tTotal Computational Time(s)\tTotal GFLOPS\tSingle Kernel Computaional Time(s)\tKernel GFLOPS\tDevice Utilisation\n");
	} else {
		fclose(file);
		file = fopen("performance_data.txt", "at");
		if (!file) {
			fprintf(stderr, "Performance Data file could not be opened or created.");
			exit(EXIT_FAILURE);
		}
	}

	atexit( cleanup );

	//for(i=32; i<=1024; i++) if (num_threads%i==0) for(j=0; j<5; j++) solveProblem(timeSteps, num_threads, i, num_partitions, num_streams, file, file1);

	solveProblem(timeSteps, num_threads, threads_per_block, num_partitions, num_streams, file, file1);

	/*for(i=2; i<21; i+=4) {
		for(j=10; j<40; j+=4) {
			for (k=0; k<10; k++) {
				 solveProblem(timeSteps, num_threads, threads_per_block, j, i, file, file1);
			}
		}
	}*/

	if (file) fclose(file);
	if (file1) fclose(file1);
}
