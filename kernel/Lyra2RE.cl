/*
 * Lyra2RE kernel implementation.
 *
 * ==========================(LICENSE BEGIN)============================
 * Copyright (c) 2014 djm34
 * Copyright (c) 2014 James Lovejoy
 * Copyright (c) 2014 phm
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *
 * ===========================(LICENSE END)=============================
 *
 * @author   James Lovejoy <jameslovejoy1@gmail.com>
 */

#pragma OPENCL EXTENSION cl_amd_printf : enable

#ifndef LYRA2RE_CL
#define LYRA2RE_CL

#if __ENDIAN_LITTLE__
#define SPH_LITTLE_ENDIAN 1
#else
#define SPH_BIG_ENDIAN 1
#endif

#define SPH_UPTR sph_u64

typedef unsigned int sph_u32;
typedef int sph_s32;
#ifndef __OPENCL_VERSION__
typedef unsigned long long sph_u64;
typedef long long sph_s64;
#else
typedef unsigned long sph_u64;
typedef long sph_s64;
#endif


#define SPH_64 1
#define SPH_64_TRUE 1

#define SPH_C32(x)    ((sph_u32)(x ## U))
#define SPH_T32(x)    ((x) & SPH_C32(0xFFFFFFFF))
#define SPH_ROTL32(x, n)   SPH_T32(((x) << (n)) | ((x) >> (32 - (n))))
#define SPH_ROTR32(x, n)   SPH_ROTL32(x, (32 - (n)))

#define SPH_C64(x)    ((sph_u64)(x ## UL))
#define SPH_T64(x)    ((x) & SPH_C64(0xFFFFFFFFFFFFFFFF))
#define SPH_ROTL64(x, n)   SPH_T64(((x) << (n)) | ((x) >> (64 - (n))))
#define SPH_ROTR64(x, n)   SPH_ROTL64(x, (64 - (n)))

#define SPH_ECHO_64 1
#define SPH_KECCAK_64 1
#define SPH_SIMD_NOCOPY 0
#define SPH_KECCAK_NOCOPY 0
#define SPH_SMALL_FOOTPRINT_GROESTL 0
#define SPH_GROESTL_BIG_ENDIAN 0
#define SPH_CUBEHASH_UNROLL 0



#include "blake256.cl"
#include "groestl256.cl"
#include "Lyra2.cl"
#include "keccak.cl"
#include "skein.cl"

#define SWAP4(x) as_uint(as_uchar4(x).wzyx)
#define SWAP8(x) as_ulong(as_uchar8(x).s76543210)

#if SPH_BIG_ENDIAN
  #define DEC64E(x) (x)
  #define DEC64BE(x) (*(const __global sph_u64 *) (x));
  #define DEC64LE(x) SWAP8(*(const __global sph_u64 *) (x));
  #define DEC32LE(x) (*(const __global sph_u32 *) (x));
#else
  #define DEC64E(x) SWAP8(x)
  #define DEC64BE(x) SWAP8(*(const __global sph_u64 *) (x));
  #define DEC64LE(x) (*(const __global sph_u64 *) (x));
#define DEC32LE(x) SWAP4(*(const __global sph_u32 *) (x));
#endif

typedef union {
  unsigned char h1[64];
  uint h4[16];
  ulong h8[8];
} hash_t;

__attribute__((reqd_work_group_size(WORKSIZE, 1, 1)))
__kernel void search(__global char* block, __global hash_t* hashes)
{

 uint gid = get_global_id(0);
  __global hash_t *hash = &(hashes[gid-get_global_offset(0)]);

	sph_u32 h[8];
	for (int i = 0; i<8; i++) { h[i] = c_IV256[i]; }  

	sph_u32 m[16];
	sph_u32 v[16];

// Compress first round ==> temporary this should be done cpu
	for (int i = 0; i < 16; i++) {
		m[i] = DEC32LE(block +i*4);
	}

// {printf("The Message 0 and 1 %08x %08x\n",m[0],m[1]);}
	for (int i = 0; i < 8; i++) {v[i] = h[i];}

	v[8] =  c_u256[0];
	v[9] =  c_u256[1];
	v[10] = c_u256[2];
	v[11] = c_u256[3];

	v[12] = c_u256[4] ^ 512;
	v[13] = c_u256[5] ^ 512;
	v[14] = c_u256[6];
	v[15] = c_u256[7];


	for (int r = 0; r < 14; r++) {
		GS(0, 4, 0x8, 0xC, 0x0);
		GS(1, 5, 0x9, 0xD, 0x2);
		GS(2, 6, 0xA, 0xE, 0x4);
		GS(3, 7, 0xB, 0xF, 0x6);
		GS(0, 5, 0xA, 0xF, 0x8);
		GS(1, 6, 0xB, 0xC, 0xA);
		GS(2, 7, 0x8, 0xD, 0xC);
		GS(3, 4, 0x9, 0xE, 0xE);
	}

	for (int i = 0; i < 16; i++) {
		 int j = i & 7;
		 h[j] ^= v[i];
	}
//////////////////////

// compress 2nd round
 m[0] = DEC32LE(block + 64);
 m[1] = DEC32LE(block + 68);
 m[2] = DEC32LE(block + 72);
 m[3] = SWAP4(gid);

	for (int i = 4; i < 16; i++) {m[i] = c_Padding[i];}

	for (int i = 0; i < 8; i++) {v[i] = h[i];}

	v[8] =  c_u256[0];
	v[9] =  c_u256[1];
	v[10] = c_u256[2];
	v[11] = c_u256[3];
	v[12] = c_u256[4] ^ 640;
	v[13] = c_u256[5] ^ 640;
	v[14] = c_u256[6];
	v[15] = c_u256[7];

	for (int r = 0; r < 14; r++) {	
		GS(0, 4, 0x8, 0xC, 0x0);
		GS(1, 5, 0x9, 0xD, 0x2);
		GS(2, 6, 0xA, 0xE, 0x4);
		GS(3, 7, 0xB, 0xF, 0x6);
		GS(0, 5, 0xA, 0xF, 0x8);
		GS(1, 6, 0xB, 0xC, 0xA);
		GS(2, 7, 0x8, 0xD, 0xC);
		GS(3, 4, 0x9, 0xE, 0xE);
	}

	for (int i = 0; i < 16; i++) {
		 int j = i & 7;
		h[j] ^= v[i];}

for (int i=0;i<8;i++) {hash->h4[i]=SWAP4(h[i]);}

 barrier(CLK_GLOBAL_MEM_FENCE);

}

// keccak256


__attribute__((reqd_work_group_size(WORKSIZE, 1, 1)))
__kernel void search1(__global hash_t* hashes)
{
  uint gid = get_global_id(0);
  __global hash_t *hash = &(hashes[gid-get_global_offset(0)]);

 		sph_u64 keccak_gpu_state[25];

		for (int i = 0; i<25; i++) {
			if (i<4) { keccak_gpu_state[i] = hash->h8[i]; }
			else    { keccak_gpu_state[i] = 0; }
		}
		keccak_gpu_state[4] = 0x0000000000000001;
		keccak_gpu_state[16] = 0x8000000000000000;

		keccak_block(keccak_gpu_state);
		for (int i = 0; i<4; i++) { hash->h8[i] = keccak_gpu_state[i]; }

  barrier(CLK_GLOBAL_MEM_FENCE);


}

/// lyra2 algo 


__attribute__((reqd_work_group_size(WORKSIZE, 1, 1)))
__kernel void search2(__global hash_t* hashes)
{
 uint gid = get_global_id(0);
  __global hash_t *hash = &(hashes[gid-get_global_offset(0)]);



  sph_u64 state[16];
//  uint2 state[16];
  for (int i = 0; i<4; i++) { state[i] = hash->h8[i];} //password
  for (int i = 0; i<4; i++) { state[i + 4] = state[i]; } //salt 

  for (int i = 0; i<8; i++) { state[i + 8] = blake2b_IV[i]; }

  //     blake2blyra x2 

  for (int i = 0; i<24; i++) { round_lyra(state); } //because 12 is not enough

  sph_u64 Matrix[96][8]; // very uncool
  /// reducedSqueezeRow0

  for (int i = 0; i < 8; i++)
  {
	  for (int j = 0; j<12; j++) { Matrix[j + 84 - 12 * i][0] = state[j]; }
	  round_lyra(state);
  }

  /// reducedSqueezeRow1

  for (int i = 0; i < 8; i++)
  {
	  for (int j = 0; j<12; j++) { state[j] ^= Matrix[j + 12 * i][0]; }
	  round_lyra(state);
	  for (int j = 0; j<12; j++) { Matrix[j + 84 - 12 * i][1] = Matrix[j + 12 * i][0] ^ state[j]; }
  }

 
  reduceDuplexRowSetup(1, 0, 2);
  reduceDuplexRowSetup(2, 1, 3);
  reduceDuplexRowSetup(3, 0, 4);
  reduceDuplexRowSetup(4, 3, 5);
  reduceDuplexRowSetup(5, 2, 6);
  reduceDuplexRowSetup(6, 1, 7);

  sph_u64 rowa;
  rowa = state[0] & 7;

  reduceDuplexRow(7, rowa, 0);
  rowa = state[0] & 7;
  reduceDuplexRow(0, rowa, 3);
  rowa = state[0] & 7;
  reduceDuplexRow(3, rowa, 6);
  rowa = state[0] & 7;
  reduceDuplexRow(6, rowa, 1);
  rowa = state[0] & 7;
  reduceDuplexRow(1, rowa, 4);
  rowa = state[0] & 7;
  reduceDuplexRow(4, rowa, 7);
  rowa = state[0] & 7;
  reduceDuplexRow(7, rowa, 2);
  rowa = state[0] & 7;
  reduceDuplexRow(2, rowa, 5);

  absorbblock(rowa);

  for (int i = 0; i<4; i++) {hash->h8[i] = state[i];} 

  barrier(CLK_GLOBAL_MEM_FENCE);

}

//skein256

__attribute__((reqd_work_group_size(WORKSIZE, 1, 1)))
__kernel void search3(__global hash_t* hashes)
{
 uint gid = get_global_id(0);
  __global hash_t *hash = &(hashes[gid-get_global_offset(0)]);


		sph_u64 h[9];
		sph_u64 t[3];
        sph_u64 dt0,dt1,dt2,dt3;
		sph_u64 p0, p1, p2, p3, p4, p5, p6, p7;
        h[8] = skein_ks_parity;
		for (int i = 0; i<8; i++) {
			h[i] = SKEIN_IV512_256[i];
			h[8] ^= h[i];}
		    
			t[0]=t12[0];
			t[1]=t12[1];
			t[2]=t12[2];

        dt0=hash->h8[0];
        dt1=hash->h8[1];
        dt2=hash->h8[2];
        dt3=hash->h8[3];

		p0 = h[0] + dt0;
		p1 = h[1] + dt1;
		p2 = h[2] + dt2;
		p3 = h[3] + dt3;
		p4 = h[4];
		p5 = h[5] + t[0];
		p6 = h[6] + t[1];
		p7 = h[7];

        #pragma unroll 
		for (int i = 1; i<19; i+=2) {Round_8_512(p0,p1,p2,p3,p4,p5,p6,p7,i);}
        p0 ^= dt0;
        p1 ^= dt1;
        p2 ^= dt2;
        p3 ^= dt3;

		h[0] = p0;
		h[1] = p1;
		h[2] = p2;
		h[3] = p3;
		h[4] = p4;
		h[5] = p5;
		h[6] = p6;
		h[7] = p7;
		h[8] = skein_ks_parity;
        
		for (int i = 0; i<8; i++) { h[8] ^= h[i]; }
		
		t[0] = t12[3];
		t[1] = t12[4];
		t[2] = t12[5];
		p5 += t[0];  //p5 already equal h[5] 
		p6 += t[1];
       
        #pragma unroll
		for (int i = 1; i<19; i+=2) { Round_8_512(p0, p1, p2, p3, p4, p5, p6, p7, i); }

		hash->h8[0]      = p0;
		hash->h8[1]      = p1;
		hash->h8[2]      = p2;
		hash->h8[3]      = p3;
	


  barrier(CLK_GLOBAL_MEM_FENCE);

}

//groestl256

__attribute__((reqd_work_group_size(WORKSIZE, 1, 1)))
__kernel void search4(__global hash_t* hashes, __global uint* output, const uint target)
{
  uint gid = get_global_id(0);
  __global hash_t *hash = &(hashes[gid-get_global_offset(0)]);

  sph_u32 message[16],state[16];

  for (int k = 0; k<8; k++) {message[k]=hash->h4[k];}

  for (int k = 9; k<15; k++) { message[k] = 0;}

  message[8] = 0x80;
  message[15] = 0x01000000;

  for (int u = 0; u<16; u++) {state[u] = message[u];}
  state[15] ^= 0x10000;

  PERM_SMALL_P(state);
  state[15] ^= 0x10000;
  PERM_SMALL_Q(message);

  for (int u = 0; u<16; u++) {state[u] ^= message[u];}

  for (int u = 0; u<16; u++) {message[u] = state[u];}

  PERM_SMALL_P(message);

  state[15] ^= message[15];
  barrier(CLK_GLOBAL_MEM_FENCE);

  bool result = ( state[15] <= target);
  if (result) {
  output[atomic_inc(output+0xFF)] = SWAP4(gid);}

}
  
 


#endif // LYRA2RE_CL