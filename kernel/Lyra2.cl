/**
 * Implementation of the Lyra2 Password Hashing Scheme (PHS).
 *
 * Author: The Lyra PHC team (http://www.lyra-kdf.net/) -- 2014.
 *
 * This software is hereby placed in the public domain.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
 
typedef unsigned char byte;
typedef long int64_t;
typedef unsigned long uint64_t;
 
void * memset(void * ptr, int value, size_t num)
{
	void * ret = ptr;
	
	while(num--)
	{
		*(char *)ptr = value;
        ptr = (char *)ptr + 1;
	}
	
	return ret;
}
 
 /*Blake2b IV Array*/
static const uint64_t blake2b_IV[8] =
{
  0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL,
  0x3c6ef372fe94f82bULL, 0xa54ff53a5f1d36f1ULL,
  0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL,
  0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL
};

/*Blake2b's rotation*/
static inline uint64_t rotr64( const uint64_t w, const unsigned c ){
    return ( w >> c ) | ( w << ( 64 - c ) );
}

/*Blake2b's G function*/
#define G(r,i,a,b,c,d) \
  do { \
    a = a + b; \
    d = rotr64(d ^ a, 32); \
    c = c + d; \
    b = rotr64(b ^ c, 24); \
    a = a + b; \
    d = rotr64(d ^ a, 16); \
    c = c + d; \
    b = rotr64(b ^ c, 63); \
  } while(0)


/*One Round of the Blake2b's compression function*/
#define ROUND_LYRA(r)  \
    G(r,0,v[ 0],v[ 4],v[ 8],v[12]); \
    G(r,1,v[ 1],v[ 5],v[ 9],v[13]); \
    G(r,2,v[ 2],v[ 6],v[10],v[14]); \
    G(r,3,v[ 3],v[ 7],v[11],v[15]); \
    G(r,4,v[ 0],v[ 5],v[10],v[15]); \
    G(r,5,v[ 1],v[ 6],v[11],v[12]); \
    G(r,6,v[ 2],v[ 7],v[ 8],v[13]); \
    G(r,7,v[ 3],v[ 4],v[ 9],v[14]);



//Block length required so Blake2's Initialization Vector (IV) is not overwritten (THIS SHOULD NOT BE MODIFIED)
#define BLOCK_LEN_BLAKE2_SAFE_INT64 8                                   //512 bits (=64 bytes, =8 uint64_t)
#define BLOCK_LEN_BLAKE2_SAFE_BYTES (BLOCK_LEN_BLAKE2_SAFE_INT64 * 8)   //same as above, in bytes


#ifdef BLOCK_LEN_BITS
        #define BLOCK_LEN_INT64 (BLOCK_LEN_BITS/64)      //Block length: 768 bits (=96 bytes, =12 uint64_t)
        #define BLOCK_LEN_BYTES (BLOCK_LEN_BITS/8)       //Block length, in bytes
#else   //default block lenght: 768 bits
        #define BLOCK_LEN_INT64 12                       //Block length: 768 bits (=96 bytes, =12 uint64_t)
        #define BLOCK_LEN_BYTES (BLOCK_LEN_INT64 * 8)    //Block length, in bytes
#endif

#ifndef N_COLS
        #define N_COLS 8                              //Number of columns in the memory matrix: fixed to 64 by default
#endif

#define ROW_LEN_INT64 (BLOCK_LEN_INT64 * N_COLS) //Total length of a row: N_COLS blocks
#define ROW_LEN_BYTES (ROW_LEN_INT64 * 8)        //Number of bytes per row
 
 /**
 * Initializes the Sponge State. The first 512 bits are set to zeros and the remainder 
 * receive Blake2b's IV as per Blake2b's specification. <b>Note:</b> Even though sponges
 * typically have their internal state initialized with zeros, Blake2b's G function
 * has a fixed point: if the internal state and message are both filled with zeros. the 
 * resulting permutation will always be a block filled with zeros; this happens because 
 * Blake2b does not use the constants originally employed in Blake2 inside its G function, 
 * relying on the IV for avoiding possible fixed points.
 * 
 * @param state         The 1024-bit array to be initialized
 */
inline void initState(uint64_t state[/*16*/]) {
    //First 512 bis are zeros
    memset(state, 0, 64);
	/*for(int i = 0; i < 64; i++)
	{
		state[i] = 0;
	}*/
	
    //Remainder BLOCK_LEN_BLAKE2_SAFE_BYTES are reserved to the IV
    state[8] = blake2b_IV[0];
    state[9] = blake2b_IV[1];
    state[10] = blake2b_IV[2];
    state[11] = blake2b_IV[3];
    state[12] = blake2b_IV[4];
    state[13] = blake2b_IV[5];
    state[14] = blake2b_IV[6];
    state[15] = blake2b_IV[7];
}

/**
 * Execute Blake2b's G function, with all 12 rounds.
 * 
 * @param v     A 1024-bit (16 uint64_t) array to be processed by Blake2b's G function
 */
inline static void blake2bLyra(uint64_t *v) {
    ROUND_LYRA(0);
    ROUND_LYRA(1);
    ROUND_LYRA(2);
    ROUND_LYRA(3);
    ROUND_LYRA(4);
    ROUND_LYRA(5);
    ROUND_LYRA(6);
    ROUND_LYRA(7);
    ROUND_LYRA(8);
    ROUND_LYRA(9);
    ROUND_LYRA(10);
    ROUND_LYRA(11);
}

/**
 * Executes a reduced version of Blake2b's G function with only one round
 * @param v     A 1024-bit (16 uint64_t) array to be processed by Blake2b's G function
 */
inline static void reducedBlake2bLyra(uint64_t *v) {
    ROUND_LYRA(0);
}

/**
 * Performs a squeeze operation, using Blake2b's G function as the 
 * internal permutation
 * 
 * @param state      The current state of the sponge 
 * @param out        Array that will receive the data squeezed
 * @param len        The number of bytes to be squeezed into the "out" array
 */
inline void squeeze(uint64_t *state, uchar *out, unsigned int len) {
    int fullBlocks = len / BLOCK_LEN_BYTES;
    uchar *ptr = out;
    int i;
    //Squeezes full blocks
    for (i = 0; i < fullBlocks; i++) {
	//memcpy(ptr, state, BLOCK_LEN_BYTES);
	for(int i2 = 0; i2 < BLOCK_LEN_BYTES; i2++)
	{
		ptr[i2] = state[i2];
	}
	blake2bLyra(state);
	ptr += BLOCK_LEN_BYTES;
    }

    //Squeezes remaining bytes
    //memcpy(ptr, state, (len % BLOCK_LEN_BYTES));
	for(i = 0; i < (len % BLOCK_LEN_BYTES); i++)
	{
		ptr[i] = state[i];
	}
}

/**
 * Performs an absorb operation for a single block (BLOCK_LEN_INT64 words
 * of type uint64_t), using Blake2b's G function as the internal permutation
 * 
 * @param state The current state of the sponge 
 * @param in    The block to be absorbed (BLOCK_LEN_INT64 words)
 */
inline void absorbBlock(uint64_t *state, const uint64_t *in) {
    //XORs the first BLOCK_LEN_INT64 words of "in" with the current state
    state[0] ^= in[0];
    state[1] ^= in[1];
    state[2] ^= in[2];
    state[3] ^= in[3];
    state[4] ^= in[4];
    state[5] ^= in[5];
    state[6] ^= in[6];
    state[7] ^= in[7];
    state[8] ^= in[8];
    state[9] ^= in[9];
    state[10] ^= in[10];
    state[11] ^= in[11];

    //Applies the transformation f to the sponge's state
    blake2bLyra(state);
}

/**
 * Performs an absorb operation for a single block (BLOCK_LEN_BLAKE2_SAFE_INT64 
 * words of type uint64_t), using Blake2b's G function as the internal permutation
 * 
 * @param state The current state of the sponge 
 * @param in    The block to be absorbed (BLOCK_LEN_BLAKE2_SAFE_INT64 words)
 */
inline void absorbBlockBlake2Safe(uint64_t *state, const uint64_t *in) {
    //XORs the first BLOCK_LEN_BLAKE2_SAFE_INT64 words of "in" with the current state
    state[0] ^= in[0];
    state[1] ^= in[1];
    state[2] ^= in[2];
    state[3] ^= in[3];
    state[4] ^= in[4];
    state[5] ^= in[5];
    state[6] ^= in[6];
    state[7] ^= in[7];

    //Applies the transformation f to the sponge's state
    blake2bLyra(state);
}

/** 
 * Performs a reduced squeeze operation for a single row, from the highest to 
 * the lowest index, using the reduced-round Blake2b's G function as the 
 * internal permutation
 * 
 * @param state     The current state of the sponge 
 * @param rowOut    Row to receive the data squeezed
 */
inline void reducedSqueezeRow0(uint64_t* state, uint64_t* rowOut) {
    uint64_t* ptrWord = rowOut + (N_COLS-1)*BLOCK_LEN_INT64; //In Lyra2: pointer to M[0][C-1]
    int i;
    //M[row][C-1-col] = H.reduced_squeeze()    
    for (i = 0; i < N_COLS; i++) {
	ptrWord[0] = state[0];
	ptrWord[1] = state[1];
	ptrWord[2] = state[2];
	ptrWord[3] = state[3];
	ptrWord[4] = state[4];
	ptrWord[5] = state[5];
	ptrWord[6] = state[6];
	ptrWord[7] = state[7];
	ptrWord[8] = state[8];
	ptrWord[9] = state[9];
	ptrWord[10] = state[10];
	ptrWord[11] = state[11];

	//Goes to next block (column) that will receive the squeezed data
	ptrWord -= BLOCK_LEN_INT64;

	//Applies the reduced-round transformation f to the sponge's state
	reducedBlake2bLyra(state);
    }
}

/** 
 * Performs a reduced duplex operation for a single row, from the highest to 
 * the lowest index, using the reduced-round Blake2b's G function as the 
 * internal permutation
 * 
 * @param state		The current state of the sponge 
 * @param rowIn		Row to feed the sponge
 * @param rowOut	Row to receive the sponge's output
 */
inline void reducedDuplexRow1(uint64_t *state, uint64_t *rowIn, uint64_t *rowOut) {
    uint64_t* ptrWordIn = rowIn;				//In Lyra2: pointer to prev
    uint64_t* ptrWordOut = rowOut + (N_COLS-1)*BLOCK_LEN_INT64; //In Lyra2: pointer to row
    int i;

    for (i = 0; i < N_COLS; i++) {

	//Absorbing "M[prev][col]"
	state[0]  ^= (ptrWordIn[0]);
	state[1]  ^= (ptrWordIn[1]);
	state[2]  ^= (ptrWordIn[2]);
	state[3]  ^= (ptrWordIn[3]);
	state[4]  ^= (ptrWordIn[4]);
	state[5]  ^= (ptrWordIn[5]);
	state[6]  ^= (ptrWordIn[6]);
	state[7]  ^= (ptrWordIn[7]);
	state[8]  ^= (ptrWordIn[8]);
	state[9]  ^= (ptrWordIn[9]);
	state[10] ^= (ptrWordIn[10]);
	state[11] ^= (ptrWordIn[11]);

	//Applies the reduced-round transformation f to the sponge's state
	reducedBlake2bLyra(state);

	//M[row][C-1-col] = M[prev][col] XOR rand
	ptrWordOut[0] = ptrWordIn[0]  ^ state[0];
	ptrWordOut[1] = ptrWordIn[1]  ^ state[1];
	ptrWordOut[2] = ptrWordIn[2]  ^ state[2];
	ptrWordOut[3] = ptrWordIn[3]  ^ state[3];
	ptrWordOut[4] = ptrWordIn[4]  ^ state[4];
	ptrWordOut[5] = ptrWordIn[5]  ^ state[5];
	ptrWordOut[6] = ptrWordIn[6]  ^ state[6];
	ptrWordOut[7] = ptrWordIn[7]  ^ state[7];
	ptrWordOut[8] = ptrWordIn[8]  ^ state[8];
	ptrWordOut[9] = ptrWordIn[9]  ^ state[9];
	ptrWordOut[10] = ptrWordIn[10] ^ state[10];
	ptrWordOut[11] = ptrWordIn[11] ^ state[11];
	
	
	//Input: next column (i.e., next block in sequence)
	ptrWordIn += BLOCK_LEN_INT64;
	//Output: goes to previous column
	ptrWordOut -= BLOCK_LEN_INT64;
    }
}

/**
 * Performs a duplexing operation over "M[rowInOut][col] [+] M[rowIn][col]" (i.e., 
 * the wordwise addition of two columns, ignoring carries between words). The
 * output of this operation, "rand", is then used to make 
 * "M[rowOut][(N_COLS-1)-col] = M[rowIn][col] XOR rand" and 
 * "M[rowInOut][col] =  M[rowInOut][col] XOR rotW(rand)", where rotW is a 64-bit 
 * rotation to the left and N_COLS is a system parameter.
 *
 * @param state          The current state of the sponge 
 * @param rowIn          Row used only as input
 * @param rowInOut       Row used as input and to receive output after rotation
 * @param rowOut         Row receiving the output
 *
 */
inline void reducedDuplexRowSetup(uint64_t *state, uint64_t *rowIn, uint64_t *rowInOut, uint64_t *rowOut) {
    uint64_t* ptrWordIn = rowIn;				//In Lyra2: pointer to prev
    uint64_t* ptrWordInOut = rowInOut;				//In Lyra2: pointer to row*
    uint64_t* ptrWordOut = rowOut + (N_COLS-1)*BLOCK_LEN_INT64; //In Lyra2: pointer to row
    int i;

    for (i = 0; i < N_COLS; i++) {
	//Absorbing "M[prev] [+] M[row*]"
	state[0]  ^= (ptrWordIn[0]  + ptrWordInOut[0]);
	state[1]  ^= (ptrWordIn[1]  + ptrWordInOut[1]);
	state[2]  ^= (ptrWordIn[2]  + ptrWordInOut[2]);
	state[3]  ^= (ptrWordIn[3]  + ptrWordInOut[3]);
	state[4]  ^= (ptrWordIn[4]  + ptrWordInOut[4]);
	state[5]  ^= (ptrWordIn[5]  + ptrWordInOut[5]);
	state[6]  ^= (ptrWordIn[6]  + ptrWordInOut[6]);
	state[7]  ^= (ptrWordIn[7]  + ptrWordInOut[7]);
	state[8]  ^= (ptrWordIn[8]  + ptrWordInOut[8]);
	state[9]  ^= (ptrWordIn[9]  + ptrWordInOut[9]);
	state[10] ^= (ptrWordIn[10] + ptrWordInOut[10]);
	state[11] ^= (ptrWordIn[11] + ptrWordInOut[11]);

	//Applies the reduced-round transformation f to the sponge's state
	reducedBlake2bLyra(state);

	//M[row][col] = M[prev][col] XOR rand
	ptrWordOut[0] = ptrWordIn[0]  ^ state[0];
	ptrWordOut[1] = ptrWordIn[1]  ^ state[1];
	ptrWordOut[2] = ptrWordIn[2]  ^ state[2];
	ptrWordOut[3] = ptrWordIn[3]  ^ state[3];
	ptrWordOut[4] = ptrWordIn[4]  ^ state[4];
	ptrWordOut[5] = ptrWordIn[5]  ^ state[5];
	ptrWordOut[6] = ptrWordIn[6]  ^ state[6];
	ptrWordOut[7] = ptrWordIn[7]  ^ state[7];
	ptrWordOut[8] = ptrWordIn[8]  ^ state[8];
	ptrWordOut[9] = ptrWordIn[9]  ^ state[9];
	ptrWordOut[10] = ptrWordIn[10] ^ state[10];
	ptrWordOut[11] = ptrWordIn[11] ^ state[11];
	
	//M[row*][col] = M[row*][col] XOR rotW(rand)
	ptrWordInOut[0]  ^= state[11];
	ptrWordInOut[1]  ^= state[0];
	ptrWordInOut[2]  ^= state[1];
	ptrWordInOut[3]  ^= state[2];
	ptrWordInOut[4]  ^= state[3];
	ptrWordInOut[5]  ^= state[4];
	ptrWordInOut[6]  ^= state[5];
	ptrWordInOut[7]  ^= state[6];
	ptrWordInOut[8]  ^= state[7];
	ptrWordInOut[9]  ^= state[8];
	ptrWordInOut[10] ^= state[9];
	ptrWordInOut[11] ^= state[10];

	//Inputs: next column (i.e., next block in sequence)
	ptrWordInOut += BLOCK_LEN_INT64;
	ptrWordIn += BLOCK_LEN_INT64;
	//Output: goes to previous column
	ptrWordOut -= BLOCK_LEN_INT64;
    }
}

/**
 * Performs a duplexing operation over "M[rowInOut][col] [+] M[rowIn][col]" (i.e., 
 * the wordwise addition of two columns, ignoring carries between words). The
 * output of this operation, "rand", is then used to make 
 * "M[rowOut][col] = M[rowOut][col] XOR rand" and 
 * "M[rowInOut][col] =  M[rowInOut][col] XOR rotW(rand)", where rotW is a 64-bit 
 * rotation to the left.
 *
 * @param state          The current state of the sponge 
 * @param rowIn          Row used only as input
 * @param rowInOut       Row used as input and to receive output after rotation
 * @param rowOut         Row receiving the output
 *
 */
inline void reducedDuplexRow(uint64_t *state, uint64_t *rowIn, uint64_t *rowInOut, uint64_t *rowOut) {
    uint64_t* ptrWordInOut = rowInOut; //In Lyra2: pointer to row*
    uint64_t* ptrWordIn = rowIn; //In Lyra2: pointer to prev
    uint64_t* ptrWordOut = rowOut; //In Lyra2: pointer to row
    int i;

    for (i = 0; i < N_COLS; i++) {

	//Absorbing "M[prev] [+] M[row*]"
	state[0]  ^= (ptrWordIn[0]  + ptrWordInOut[0]);
	state[1]  ^= (ptrWordIn[1]  + ptrWordInOut[1]);
	state[2]  ^= (ptrWordIn[2]  + ptrWordInOut[2]);
	state[3]  ^= (ptrWordIn[3]  + ptrWordInOut[3]);
	state[4]  ^= (ptrWordIn[4]  + ptrWordInOut[4]);
	state[5]  ^= (ptrWordIn[5]  + ptrWordInOut[5]);
	state[6]  ^= (ptrWordIn[6]  + ptrWordInOut[6]);
	state[7]  ^= (ptrWordIn[7]  + ptrWordInOut[7]);
	state[8]  ^= (ptrWordIn[8]  + ptrWordInOut[8]);
	state[9]  ^= (ptrWordIn[9]  + ptrWordInOut[9]);
	state[10] ^= (ptrWordIn[10] + ptrWordInOut[10]);
	state[11] ^= (ptrWordIn[11] + ptrWordInOut[11]);

	//Applies the reduced-round transformation f to the sponge's state
	reducedBlake2bLyra(state);

	//M[rowOut][col] = M[rowOut][col] XOR rand
	ptrWordOut[0] ^= state[0];
	ptrWordOut[1] ^= state[1];
	ptrWordOut[2] ^= state[2];
	ptrWordOut[3] ^= state[3];
	ptrWordOut[4] ^= state[4];
	ptrWordOut[5] ^= state[5];
	ptrWordOut[6] ^= state[6];
	ptrWordOut[7] ^= state[7];
	ptrWordOut[8] ^= state[8];
	ptrWordOut[9] ^= state[9];
	ptrWordOut[10] ^= state[10];
	ptrWordOut[11] ^= state[11];

	//M[rowInOut][col] = M[rowInOut][col] XOR rotW(rand)
	ptrWordInOut[0] ^= state[11];
	ptrWordInOut[1] ^= state[0];
	ptrWordInOut[2] ^= state[1];
	ptrWordInOut[3] ^= state[2];
	ptrWordInOut[4] ^= state[3];
	ptrWordInOut[5] ^= state[4];
	ptrWordInOut[6] ^= state[5];
	ptrWordInOut[7] ^= state[6];
	ptrWordInOut[8] ^= state[7];
	ptrWordInOut[9] ^= state[8];
	ptrWordInOut[10] ^= state[9];
	ptrWordInOut[11] ^= state[10];

	//Goes to next block
	ptrWordOut += BLOCK_LEN_INT64;
	ptrWordInOut += BLOCK_LEN_INT64;
	ptrWordIn += BLOCK_LEN_INT64;
    }
}

void * memcpy(void * dst, const void * src, size_t count)
{
        void * ret = dst;

        while(count--) 
		{
            *(char *)dst = *(char *)src;
            dst = (char *)dst + 1;
            src = (char *)src + 1;
        }

        return ret;
}



/**
 * Executes Lyra2 based on the G function from Blake2b. This version supports salts and passwords
 * whose combined length is smaller than the size of the memory matrix, (i.e., (nRows x nCols x b) bits,
 * where "b" is the underlying sponge's bitrate). In this implementation, the "basil" is composed by all
 * integer parameters (treated as type "unsigned int") in the order they are provided, plus the value
 * of nCols, (i.e., basil = kLen || pwdlen || saltlen || timeCost || nRows || nCols).
 *
 * @param K The derived key to be output by the algorithm
 * @param kLen Desired key length
 * @param pwd User password
 * @param pwdlen Password length
 * @param salt Salt
 * @param saltlen Salt length
 * @param timeCost Parameter to determine the processing time (T)
 * @param nRows Number or rows of the memory matrix (R)
 * @param nCols Number of columns of the memory matrix (C)
 *
 * @return 0 if the key is generated correctly; -1 if there is an error (usually due to lack of memory for allocation)
 */
 
#define numRows 8
#define numCols 8
#define keyLength 64
#define passwordLength 64
#define saltLength 64

void LYRA2(uchar *K, const uchar *pwd, const uchar *salt, unsigned int timeCost) {

    //============================= Basic variables ============================//
    int64_t row = 2; //index of row to be processed
    int64_t prev = 1; //index of prev (last row ever computed/modified)
    int64_t rowa = 0; //index of row* (a previous row, deterministically picked during Setup and randomly picked while Wandering)
    int64_t tau; //Time Loop iterator
    int64_t step = 1; //Visitation step (used during Setup and Wandering phases)
    int64_t window = 2; //Visitation window (used to define which rows can be revisited during Setup)
    int64_t gap = 1; //Modifier to the step, assuming the values 1 or -1
    int64_t i; //auxiliary iteration counter
    //==========================================================================/
	
	const unsigned int nRows = numRows;
	const unsigned int nCols = numCols;
	const unsigned int kLen = keyLength;
	const unsigned int pwdlen = passwordLength;
	const unsigned int saltlen = saltLength;
		
	uint64_t wholeMatrixarr[(numRows * ROW_LEN_BYTES) / 8];
	uint64_t memMatrixarr[numRows];
	uint64_t statearr[16];
	
	uint64_t *uintpointers;
	uintpointers = memMatrixarr;
	
	uint64_t *wholeMatrix;// = wholeMatrixarr;
	uint64_t **memMatrix;// = memMatrixarr;
	uint64_t *state;// = statearr;
	
	state = statearr;
	wholeMatrix = wholeMatrixarr;
	memMatrix = &uintpointers;

    //========== Initializing the Memory Matrix and pointers to it =============//
    //Tries to allocate enough space for the whole memory matrix
    //i = (int64_t) ((int64_t) nRows * (int64_t) ROW_LEN_BYTES);
    //uint64_t *wholeMatrix = malloc(i);

    //Allocates pointers to each row of the matrix
    //uint64_t **memMatrix = malloc(nRows * sizeof (uint64_t*));
	
    //Places the pointers in the correct positions
    uint64_t *ptrWord = wholeMatrix;
    for (i = 0; i < nRows; i++) {
      memMatrix[i] = ptrWord;
      ptrWord += ROW_LEN_INT64;
    }
    //==========================================================================/

    //============= Getting the password + salt + basil padded with 10*1 ===============//
    //OBS.:The memory matrix will temporarily hold the password: not for saving memory,
    //but this ensures that the password copied locally will be overwritten as soon as possible

    //First, we clean enough blocks for the password, salt, basil and padding
    uint64_t nBlocksInput = ((saltlen + pwdlen + 6 * sizeof (int)) / BLOCK_LEN_BLAKE2_SAFE_BYTES) + 1;
    uchar *ptrByte = (uchar*) wholeMatrix;
	
	memset(ptrByte, 0, nBlocksInput * BLOCK_LEN_BLAKE2_SAFE_BYTES);
	/*for(i = 0; i < nBlocksInput * BLOCK_LEN_BLAKE2_SAFE_BYTES; i++)
	{
		ptrByte[i] = 0;
	}*/

    //Prepends the password
    memcpy(ptrByte, pwd, pwdlen);
	/*for(i = 0; i < pwdlen; i++)
	{
		ptrByte[i] = pwd[i];
	}*/
	//ptrByte = as_uchar(pwd);
    ptrByte += pwdlen;

    //Concatenates the salt
    memcpy(ptrByte, salt, saltlen);
	/*for(i = 0; i < saltlen; i++)
	{
		ptrByte[i] = salt[i];
	}*/
    ptrByte += saltlen;

    //Concatenates the basil: every integer passed as parameter, in the order they are provided by the interface
    memcpy(ptrByte, &kLen, sizeof (int));
	/*for(i = 0; i < sizeof(int); i++)
	{
		ptrByte[i] = &kLen[i];
	}*/
    ptrByte += sizeof (int);
    memcpy(ptrByte, &pwdlen, sizeof (int));
	/*for(i = 0; i < sizeof(int); i++)
	{
		ptrByte[i] = pwdlen[i];
	}*/
    ptrByte += sizeof (int);
    memcpy(ptrByte, &saltlen, sizeof (int));
	/*for(i = 0; i < sizeof(int); i++)
	{
		ptrByte[i] = saltlen[i];
	}*/
    ptrByte += sizeof (int);
    memcpy(ptrByte, &timeCost, sizeof (int));
	/*for(i = 0; i < sizeof(int); i++)
	{
		ptrByte[i] = timeCost[i];
	}*/
    ptrByte += sizeof (int);
    memcpy(ptrByte, &nRows, sizeof (int));
	/*for(i = 0; i < sizeof(int); i++)
	{
		ptrByte[i] = nRows[i];
	}*/
    ptrByte += sizeof (int);
    memcpy(ptrByte, &nCols, sizeof (int));
	/*for(i = 0; i < sizeof(int); i++)
	{
		ptrByte[i] = nCols[i];
	}*/
    ptrByte += sizeof (int);

    //Now comes the padding
    *ptrByte = 0x80; //first byte of padding: right after the password
    ptrByte = (uchar*) wholeMatrix; //resets the pointer to the start of the memory matrix
    ptrByte += nBlocksInput * BLOCK_LEN_BLAKE2_SAFE_BYTES - 1; //sets the pointer to the correct position: end of incomplete block
    *ptrByte ^= 0x01; //last byte of padding: at the end of the last incomplete block
    //==========================================================================/

    //======================= Initializing the Sponge State ====================//
    //Sponge state: 16 uint64_t, BLOCK_LEN_INT64 words of them for the bitrate (b) and the remainder for the capacity (c)
    //uint64_t *state = malloc(16 * sizeof (uint64_t));
    initState(state);
    //==========================================================================/

    //================================ Setup Phase =============================//
    //Absorbing salt, password and basil: this is the only place in which the block length is hard-coded to 512 bits
    ptrWord = wholeMatrix;
    for (i = 0; i < nBlocksInput; i++) {
      absorbBlockBlake2Safe(state, ptrWord); //absorbs each block of pad(pwd || salt || basil)
      ptrWord += BLOCK_LEN_BLAKE2_SAFE_BYTES; //goes to next block of pad(pwd || salt || basil)
    }

    //Initializes M[0] and M[1]
    reducedSqueezeRow0(state, memMatrix[0]); //The locally copied password is most likely overwritten here
    reducedDuplexRow1(state, memMatrix[0], memMatrix[1]);

    do {
      //M[row] = rand; //M[row*] = M[row*] XOR rotW(rand)
      reducedDuplexRowSetup(state, memMatrix[prev], memMatrix[rowa], memMatrix[row]);


      //updates the value of row* (deterministically picked during Setup))
      rowa = (rowa + step) & (window - 1);
      //update prev: it now points to the last row ever computed
      prev = row;
      //updates row: goes to the next row to be computed
      row++;

      //Checks if all rows in the window where visited.
      if (rowa == 0) {
      step = window + gap; //changes the step: approximately doubles its value
      window *= 2; //doubles the size of the re-visitation window
      gap = -gap; //inverts the modifier to the step
    }

    } while (row < nRows);
    //==========================================================================/

    //============================ Wandering Phase =============================//
    row = 0; //Resets the visitation to the first row of the memory matrix
    for (tau = 1; tau <= timeCost; tau++) {
    	//Step is approximately half the number of all rows of the memory matrix for an odd tau; otherwise, it is -1
    	step = (tau % 2 == 0) ? -1 : nRows / 2 - 1;
    	do {
  	    //Selects a pseudorandom index row*
  	    //------------------------------------------------------------------------------------------
  	    //rowa = ((unsigned int)state[0]) & (nRows-1);	//(USE THIS IF nRows IS A POWER OF 2)
  	    rowa = ((unsigned int) (state[0])) % nRows; //(USE THIS FOR THE "GENERIC" CASE)
  	    //------------------------------------------------------------------------------------------

  	    //Performs a reduced-round duplexing operation over M[row*] XOR M[prev], updating both M[row*] and M[row]
  	    reducedDuplexRow(state, memMatrix[prev], memMatrix[rowa], memMatrix[row]);

  	    //update prev: it now points to the last row ever computed
  	    prev = row;

  	    //updates row: goes to the next row to be computed
  	    //------------------------------------------------------------------------------------------
  	    //row = (row + step) & (nRows-1);	//(USE THIS IF nRows IS A POWER OF 2)
  	    row = (row + step) % nRows; //(USE THIS FOR THE "GENERIC" CASE)
  	    //------------------------------------------------------------------------------------------

      } while (row != 0);
    }
    //==========================================================================/

    //============================ Wrap-up Phase ===============================//
    //Absorbs the last block of the memory matrix
    absorbBlock(state, memMatrix[rowa]);

    //Squeezes the key
    squeeze(state, K, kLen);
    //==========================================================================/
}
