#ifndef ALLIUM_H
#define ALLIUM_H

#include "miner.h"

extern int allium_test(unsigned char *pdata, const unsigned char *ptarget,
			uint32_t nonce);
extern void allium_regenhash(struct work *work);

#endif /* ALLIUM_H */
