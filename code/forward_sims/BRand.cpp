#include "BRand.hpp"

// non-inline function definitions and static member definitions cannot
// reside in header file because of the risk of multiple declarations


/** Singleton */
BRand BRand::Controller;
/** Static variables */
unsigned long BRand::State[N] = {0x0UL};
int BRand::P = 0;

/** generate new state vector */
void BRand::gen_state() {
	for (int i = 0; i < (N - M); ++i){
		State[i] = State[i + M] ^ twiddle(State[i], State[i + 1]);
	}
	for (int i = N - M; i < (N - 1); ++i){
		State[i] = State[i + M - N] ^ twiddle(State[i], State[i + 1]);
	}
	State[N - 1] = State[M - 1] ^ twiddle(State[N - 1], State[0]);
	P = 0; // reset position
}

/** init by 32 bit seed */
void BRand::seed(unsigned long s) {
	State[0] = s & 0xFFFFFFFFUL; // for > 32 bit machines
	for (int i = 1; i < N; ++i) {
		State[i] = 1812433253UL * (State[i - 1] ^ (State[i - 1] >> 30)) + i;
		// see Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier
		// in the previous versions, MSBs of the seed affect only MSBs of the array state
		// 2002/01/09 modified by Makoto Matsumoto
		State[i] &= 0xFFFFFFFFUL; // for > 32 bit machines
	}
	P = N; // force gen_state() to be called for next random number
}
