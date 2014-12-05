

/*Blake2b IV Array*/
__constant static const sph_u64 blake2b_IV[8] =
{
  0x6a09e667f3bcc908ULL, 0xbb67ae8584caa73bULL,
  0x3c6ef372fe94f82bULL, 0xa54ff53a5f1d36f1ULL,
  0x510e527fade682d1ULL, 0x9b05688c2b3e6c1fULL,
  0x1f83d9abfb41bd6bULL, 0x5be0cd19137e2179ULL
};

/*Blake2b's rotation*/
static inline sph_u64 rotr64( const sph_u64 w, const unsigned c ){
    return ( w >> c ) | ( w << ( 64 - c ) );
}

/*Blake2b's G function*/
#define G(a,b,c,d) \
  do { \
a += b; d ^= a; d = SPH_ROTR64(d, 32); \
c += d; b ^= c; b = SPH_ROTR64(b, 24); \
a += b; d ^= a; d = SPH_ROTR64(d, 16); \
c += d; b ^= c; b = SPH_ROTR64(b, 63); \
  } while(0)


/*One Round of the Blake2b's compression function*/
#define round_lyra(v)  \
 do { \
    G(v[ 0],v[ 4],v[ 8],v[12]); \
    G(v[ 1],v[ 5],v[ 9],v[13]); \
    G(v[ 2],v[ 6],v[10],v[14]); \
    G(v[ 3],v[ 7],v[11],v[15]); \
    G(v[ 0],v[ 5],v[10],v[15]); \
    G(v[ 1],v[ 6],v[11],v[12]); \
    G(v[ 2],v[ 7],v[ 8],v[13]); \
    G(v[ 3],v[ 4],v[ 9],v[14]); \
 } while(0)


#define reduceDuplexRowSetup(rowIn, rowInOut, rowOut) \
   { \
	for (int i = 0; i < 8; i++) \
				{ \
\
		for (int j = 0; j < 12; j++) {state[j] ^= Matrix[12 * i + j][rowIn] + Matrix[12 * i + j][rowInOut];} \
		round_lyra(state); \
		for (int j = 0; j < 12; j++) {Matrix[j + 84 - 12 * i][rowOut] = Matrix[12 * i + j][rowIn] ^ state[j];} \
\
		Matrix[0 + 12 * i][rowInOut] ^= state[11]; \
		Matrix[1 + 12 * i][rowInOut] ^= state[0]; \
		Matrix[2 + 12 * i][rowInOut] ^= state[1]; \
		Matrix[3 + 12 * i][rowInOut] ^= state[2]; \
		Matrix[4 + 12 * i][rowInOut] ^= state[3]; \
		Matrix[5 + 12 * i][rowInOut] ^= state[4]; \
		Matrix[6 + 12 * i][rowInOut] ^= state[5]; \
		Matrix[7 + 12 * i][rowInOut] ^= state[6]; \
		Matrix[8 + 12 * i][rowInOut] ^= state[7]; \
		Matrix[9 + 12 * i][rowInOut] ^= state[8]; \
		Matrix[10 + 12 * i][rowInOut] ^= state[9]; \
		Matrix[11 + 12 * i][rowInOut] ^= state[10]; \
				} \
 \
   } 

#define reduceDuplexRow(rowIn, rowInOut, rowOut) \
  { \
	 for (int i = 0; i < 8; i++) \
	 	 	 	 	 { \
		 for (int j = 0; j < 12; j++) \
			 state[j] ^= Matrix[12 * i + j][rowIn] + Matrix[12 * i + j][rowInOut]; \
 \
		 round_lyra(state); \
		 for (int j = 0; j < 12; j++) {Matrix[j + 12 * i][rowOut] ^= state[j];} \
\
		 Matrix[0 + 12 * i][rowInOut] ^= state[11]; \
		 Matrix[1 + 12 * i][rowInOut] ^= state[0]; \
		 Matrix[2 + 12 * i][rowInOut] ^= state[1]; \
		 Matrix[3 + 12 * i][rowInOut] ^= state[2]; \
		 Matrix[4 + 12 * i][rowInOut] ^= state[3]; \
		 Matrix[5 + 12 * i][rowInOut] ^= state[4]; \
		 Matrix[6 + 12 * i][rowInOut] ^= state[5]; \
		 Matrix[7 + 12 * i][rowInOut] ^= state[6]; \
		 Matrix[8 + 12 * i][rowInOut] ^= state[7]; \
		 Matrix[9 + 12 * i][rowInOut] ^= state[8]; \
		 Matrix[10 + 12 * i][rowInOut] ^= state[9]; \
		 Matrix[11 + 12 * i][rowInOut] ^= state[10]; \
	 	 	 	 	 } \
 \
  } 
#define absorbblock(in)  { \
	state[0] ^= Matrix[0][in]; \
	state[1] ^= Matrix[1][in]; \
	state[2] ^= Matrix[2][in]; \
	state[3] ^= Matrix[3][in]; \
	state[4] ^= Matrix[4][in]; \
	state[5] ^= Matrix[5][in]; \
	state[6] ^= Matrix[6][in]; \
	state[7] ^= Matrix[7][in]; \
	state[8] ^= Matrix[8][in]; \
	state[9] ^= Matrix[9][in]; \
	state[10] ^= Matrix[10][in]; \
	state[11] ^= Matrix[11][in]; \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
	round_lyra(state); \
  } 