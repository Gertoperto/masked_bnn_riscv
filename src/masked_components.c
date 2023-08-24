#include <stdbool.h>
#include <stdint.h>

#ifdef NO_ASM
// the compiler might potentially unmask the variables, this function is only really to test the correctness of our scheme
// this is not the case for the assembly components.
void masked_AND(bool s0_a, bool s1_a, bool s0_b, bool s1_b, bool* s0_out, bool* s1_out, bool r) {
	*s0_out = s1_a & s1_b ^ (s1_a & s0_b ^ (s0_a & s1_b ^ (s0_a & s0_b ^ r))); // @suppress("Suggested parenthesis around expression")
	*s1_out = r; // if the calling function holds on to r, we can get rid of s1_out.
}

int64_t m_and(uint32_t s0_a, uint32_t s1_a, uint32_t s0_b, uint32_t s1_b) {
#ifndef RNG_OFF
	uint32_t r = rand();
#else
	uint32_t r = 0;
#endif
	uint32_t s0_out = s1_a & s1_b ^ (s1_a & s0_b ^ (s0_a & s1_b ^ (s0_a & s0_b ^ r))); // @suppress("Suggested parenthesis around expression")
	uint64_t out = s0_out; out ^= ((uint64_t) r) << 32;
	return out;
}

// a regular one bit adder is as follows: it takes as input a carry bit, and two sum bits
// so a masked version takes masked shares of all inputs
// and outputs a masked carry_out, and masked sum bit.
// the three masked and operations in the adder each require one bit of randomness, meaning our function takes three rbits
void masked_1bit_ADD(bool s0_a, bool s1_a, bool s0_b, bool s1_b, bool s0_c, bool s1_c, bool* s0_cout, bool* s1_cout, bool* s0_sum, bool* s1_sum, bool r1, bool r2, bool r3)  {
	// sum bit is linear, sum = x xor y xor carry
	*s0_sum = s0_a ^ s0_b ^ s0_c;
	*s1_sum = s1_a ^ s1_b ^ s1_c;
	// carry bit is not linear (uses AND.)
	// cout = (a & b) xor (b & c) xor (c & a).
	//
	bool s0_ab, s1_ab, s0_bc, s1_bc, s0_ca, s1_ca;
	masked_AND(s0_a, s1_a, s0_b, s1_b, &s0_ab, &s1_ab, r1);
	masked_AND(s0_b, s1_b, s0_c, s1_c, &s0_bc, &s1_bc, r2);
	masked_AND(s0_c, s1_c, s0_a, s1_a, &s0_ca, &s1_ca, r3);

	*s0_cout = s0_ab ^ s0_bc ^ s0_ca;
	*s1_cout = s1_ab ^ s1_bc ^ s1_ca;
	// note for assembly implementation: do not need all three ab, bc, ca at once, can ab^bc before calculating c&a to save a temporary.
}

// n <= 32
// can we recycle r1,r2,r3 for each single bit addition? the shares are already masked :thinking:
// this n bit adder is essentially a ripple-carry adder.
// the use case of this adder is adding 8-bit integers to a larger integer, for all the bits beyond the 8th bit we could use a simpler adder that only accepst (a,0,c) (no b). ->
// sum_bit = a xor 0 xor c = a xor c
// carry_out = (a & 0) xor (0 & c) xor (c & a) = c & a
// only need a single random bit for that special case
// is that special case called a partial adder? ... nice.
//
// this special case might not work when adding negative numbers though, since they should be sign extended to whatever they get added to
// though then the special case becomes
// sum_bit = a xor 1 xor c = !a xor c
// carry_out = (a & 1) xor (1 & c) xor (c & a) = a xor c xor (c & a) = a or c
//
// all bits more significant than n are 0.
int64_t m_add_32b(uint32_t s0_a, uint32_t s1_a, uint32_t s0_b, uint32_t s1_b, bool s0_c, bool s1_c) {
	// carry out of one round is connected to carry in of the next round.
#ifndef RNG_OFF
	int r1 = rand(), r2 = rand(), r3 = rand();
#else
	int r1 = 0, r2 = 0, r3 = 0;
#endif

	bool s0_carry = s0_c, s1_carry = s1_c, s0_sum_out, s1_sum_out, s0_a_in, s1_a_in, s0_b_in, s1_b_in, r1_in, r2_in, r3_in;
	uint32_t n = 32, s0_sum_out_i = 0, s1_sum_out_i = 0;
	for (int i = 0; i < n; ++i) {
		// grab one bit from the input (LSB), and shift it to the right
		s0_a_in = s0_a & 1; s1_a_in = s1_a & 1; s0_b_in = s0_b & 1; s1_b_in = s1_b & 1;
		s0_a >>= 1; s1_a >>= 1; s0_b >>= 1; s1_b >>= 1;
		// grab one random bit from our random ints, and shift to the right to prep for next round
		r1_in = r1 & 1; r2_in = r2 & 1; r3_in = r3 & 1;
		r1 >>= 1; r2 >>= 1; r3 >>= 1;
		// feed all the inputs to our masked 1bit full adder, notably the carry gets copied from the previous round, and its address is reused for the next round.
		masked_1bit_ADD(s0_a_in, s1_a_in, s0_b_in, s1_b_in, s0_carry, s1_carry, &s0_carry, &s1_carry, &s0_sum_out, &s1_sum_out, r1_in, r2_in, r3_in);
		// shift the output one bit to the left, and XOR the answer bit into place
		s0_sum_out_i ^= s0_sum_out << i; s1_sum_out_i ^= s1_sum_out << i;
	}
	// add the carry as the MSB. (NO)
	// s0_sum_out_i ^= s0_carry << n_copy; s1_sum_out_i ^= s1_carry << n_copy;
	uint64_t out = s0_sum_out_i; out ^= ((uint64_t) s1_sum_out_i) << 32;
	return out;
}

int64_t m_popcount(uint32_t s0_in, uint32_t s1_in) {

	uint32_t s0_temp, s1_temp;
	uint64_t temp;
	// add single bits to get 16x 2-bit accumulators
	temp = m_add_32b(s0_in & 0x55555555, s1_in & 0x55555555, (s0_in >> 1) & 0x55555555, (s1_in >> 1) & 0x55555555, 0, 0);
	s0_temp = temp; s1_temp = temp >> 32;
	// add 2-bit accumulators to get 8x 4-bit accumulators
	temp = m_add_32b(s0_temp & 0x33333333, s1_temp & 0x33333333, (s0_temp >> 2) & 0x33333333, (s1_temp >> 2) & 0x33333333, 0, 0);
	s0_temp = temp; s1_temp = temp >> 32;
	// add 4-bit accumulators to get 4x 8-bit accumulators
	temp = m_add_32b(s0_temp & 0x07070707, s1_temp & 0x07070707, (s0_temp >> 4) & 0x07070707, (s1_temp >> 4) & 0x07070707, 0, 0);
	s0_temp = temp; s1_temp = temp >> 32;
	// add 8-bit accumulators to get 2x 16-bit accumulators
	temp = m_add_32b(s0_temp & 0x000F000F, s1_temp & 0x000F000F, (s0_temp >> 8) & 0x000F000F, (s1_temp >> 8) & 0x000F000F, 0, 0);
	s0_temp = temp; s1_temp = temp >> 32;
	// add 16-bit accumulators, our adder leaves the 16 most significant bits 0, no need for more masking after.
	temp = m_add_32b(s0_temp & 0x0000001F, s1_temp & 0x0000001F, (s0_temp >> 16) & 0x0000001F, (s1_temp >> 16) & 0x0000001F, 0, 0);
	return temp;
}
#endif
