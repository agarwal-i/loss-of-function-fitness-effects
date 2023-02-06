#ifndef BRAND_HPP
#define BRAND_HPP

// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// All rights reserved.
// Some source from by Jasper Bedaux 2003/1/1 (see http://www.bedaux.net/mtrand/)
// and Isaku Wada, 2002/01/09
// Here we just create a pretty and easy to use class
// based on their work
// Mersenne Twister random number generator

/**
* @brief Random class for the brainable engine
*/
class BRand {
public:
	/** Singleton controller */
	static BRand Controller;

private:
	/**
	* @brief default constructor private to use singleton
	*/
	BRand() {
		seed(5489UL);
	}

	/**
	* @brief Destructor
	*/
	virtual ~BRand() {} // destructor

public:
	/**
	* @brief Set seed with 32 bit integer
	* @param seed
	*/
	void seed(unsigned long);

	/**
	* @brief overload operator() to make this a generator (functor)
	*/
	unsigned long operator()(){
		return rand_int32();
	}

	/**
	* @brief Get number in [0, 1]
	*/
	double nextClosed() {
		return static_cast<double>(rand_int32()) * (1. / 4294967295.); // divided by 2^32 - 1
	}

	/**
	* @brief Get number in (0, 1)
	*/
	double nextOpened() {
		return (static_cast<double>(rand_int32()) + .5) * (1. / 4294967296.);// divided by 2^32
	}


protected:
	/**
	* @brief To have the next generated number
	* @brief used by derived classes, otherwise not accessible; use the ()-operator
	*/
	inline unsigned long rand_int32() { // generate 32 bit random int
		if (P == N) gen_state(); // new state vector needed
		// gen_state() is split off to be non-inline, because it is only called once
		// in every 624 calls and otherwise irand() would become too big to get inlined
		unsigned long x = State[P++];
		x ^= (x >> 11);
		x ^= (x << 7) & 0x9D2C5680UL;
		x ^= (x << 15) & 0xEFC60000UL;
		return x ^ (x >> 18);
	}

private:
	static const int N = 624;
	static const int M = 397;

	// the variables below are static (no duplicates can exist)
	static unsigned long State[N]; // state vector array
	static int P; // position in state array

	/**
	* @biref used by gen_state()
	* @brief private functions used to generate the pseudo random numbers
	*/
	unsigned long twiddle(unsigned long u, unsigned long v) {
		return (((u & 0x80000000UL) | (v & 0x7FFFFFFFUL)) >> 1)
			^ ((v & 1UL) ? 0x9908B0DFUL : 0x0UL);
	}
	/**
	* @brief generate new state
	*/
	void gen_state();
};



#endif
