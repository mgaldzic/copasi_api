/* Begin CVS Header
   $Source: /fs/turing/cvs/copasi_dev/copasi/randomGenerator/Cr250.cpp,v $
   $Revision: 1.9 $
   $Name: Build-33 $
   $Author: shoops $
   $Date: 2006/04/27 01:31:00 $
   End CVS Header */

// Copyright � 2005 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/******************************************************************************
 *  Module:  Cr250.cpp
 *  Description: implements R250 random number generator, from S.
 *  Kirkpatrick and E. Stoll, Journal of Computational Physics, 40,
 *  p. 517 (1981).
 *  Written by:    W. L. Maier
 *
 * METHOD...
 *  16 parallel copies of a linear shift register with
 *  period 2^250 - 1.  FAR longer period than the usual
 *  linear congruent generator, and commonly faster as
 *  well.  (For details see the above paper, and the
 *  article in DDJ referenced below.)
 *
 * HISTORY...
 *  Sep 02: Converted to c++ class by Stefan Hoops
 *  Sep 92: Number returned by dr250() is in range [0,1) instead
 *   of [0,1], so for example a random angle in the
 *   interval [0, 2*PI) can be calculated
 *   conveniently.  (J. R. Van Zandt <jrv@mbunix.mitre.org>)
 *  Aug 92: Initialization is optional.  Default condition is
 *   equivalent to initializing with the seed 12345,
 *   so that the sequence of random numbers begins:
 *   1173, 53403, 52352, 35341...  (9996 numbers
 *   skipped) ...57769, 14511, 46930, 11942, 7978,
 *   56163, 46506, 45768, 21162, 43113...  Using ^=
 *   operator rather than two separate statements.
 *   Initializing with own linear congruent
 *   pseudorandom number generator for portability.
 *   Function prototypes moved to a header file.
 *   Implemented r250n() to generate numbers
 *   uniformly distributed in a specific range
 *   [0,n), because the common expedient rand()%n is
 *   incorrect.  (J. R. Van Zandt <jrv@mbunix.mitre.org>)
 *  May 91: Published by W. L. Maier, "A Fast Pseudo Random Number
 *   Generator", Dr. Dobb's Journal #176.
 ******************************************************************************/

#include "copasi.h"
#include "CRandom.h"

Cr250::Cr250(unsigned C_INT32 seed):
    CRandom(),
    mIndex(0)
{
  setModulus(65535);
  initialize(seed);
}
Cr250::~Cr250(){}

void Cr250::initialize(unsigned C_INT32 seed)
{
  /*--------------------------------------------------------------------------*/
  int j, k;
  unsigned int mask;
  unsigned int msb;
  /*--------------------------------------------------------------------------*/

  mIndex = 0;
  mSeed = seed;

  /* Fill the r250 buffer with 15-bit values */
  for (j = 0; j < 250; j++)
    mBuffer[j] = myrand();

  /* Set some of the MS bits to 1 */
  for (j = 0; j < 250; j++)
    if (myrand() > 16384)
      mBuffer[j] |= 0x8000;

  msb = 0x8000;       /* To turn on the diagonal bit   */
  mask = 0xffff;      /* To turn off the leftmost bits */

  for (j = 0; j < 16; j++)
    {
      k = 11 * j + 3;         /* Select a word to operate on        */
      mBuffer[k] &= mask;     /* Turn off bits left of the diagonal */
      mBuffer[k] |= msb;      /* Turn on the diagonal bit           */
      mask >>= 1;
      msb >>= 1;
    }
  return;
}

unsigned C_INT32 Cr250::getRandomU()
{return r250();}

C_INT32 Cr250::getRandomS()
{return r250();}

C_FLOAT64 Cr250::getRandomCC()
{return r250() * mModulusInv;}

C_FLOAT64 Cr250::getRandomCO()
{return dr250();}

C_FLOAT64 Cr250::getRandomOO()
{return (r250() + .5) * mModulusInv1;}

unsigned C_INT32 Cr250::r250(void)
{
  register C_INT16 j;

  if (mIndex > 146)
    j = mIndex - 147;      /* Wrap pointer around */
  else
    j = mIndex + 103;

  mNumberU = mBuffer[mIndex] ^= mBuffer[j];

  if (mIndex > 248)      /* Increment pointer for next time */
    mIndex = 0;
  else
    mIndex++;

  return mNumberU;
}

unsigned C_INT32 Cr250::r250n(const unsigned C_INT16 & max)
{
  register unsigned C_INT16 limit;

  limit = (65535U / max) * max;

  do
    {
      r250();
      r250(); // Why a second call?
    }
  while (mNumberU >= limit);

  return mNumberU % max;
}

C_FLOAT64 Cr250::dr250()
{
  mNumberU = r250();
  return mNumberU / 65536.;   /* Return a number in [0.0 to 1.0) */
}

unsigned C_INT16 Cr250::myrand()
{
  mSeed = mSeed * 0x015a4e35L + 1;
  return (mSeed >> 16)&0x7fff;
}
