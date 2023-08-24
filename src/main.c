/* Copyright 2019 SiFive, Inc */
/* SPDX-License-Identifier: Apache-2.0 */
#include <stdio.h>
#include <stdint.h>     /* for uinptr_t */
#include <stdlib.h>
#include <limits.h>
#include <stdbool.h>
#include <string.h>
#include <metal/hpm.h>
#include "static.h"

// #define RNG_OFF
// #define NO_ASM

#ifndef STATIC_H
#define CEILING(x,y) (((x) + (y) - 1) / (y))
#define INPUT_PIXELS 28*28
#define INPUT_PIXELS_WORD CEILING(INPUT_PIXELS, 32)
#define L1_NODES 1010
#define L1_NODES_WORD CEILING(L1_NODES, 32)
#define L2_NODES 1010
#define L2_NODES_WORD CEILING(L2_NODES, 32)
#define L3_NODES 1010
#define L3_NODES_WORD CEILING(L3_NODES, 32)
#define OUTPUT_NODES 10
#endif

#define VOLATILE_START_ADDR (uint32_t *) 0x80000000
#define VOLATILE_END_ADDR (uint32_t *)   0x80004000

#define POPCOUNT_RWORDS 3 * 5
#define ADD_RWORDS 3

#define MAX_ADDITIONS 100

#ifndef NO_ASM
extern int64_t m_add_32b(int32_t s0_a, int32_t s1_a, int32_t s0_b, int32_t s1_b, int32_t s0_c, int32_t s1_c);
extern int64_t m_popcount(int32_t s0_a, int32_t s1_a);
extern int64_t m_and(int32_t s0_a, int32_t s1_a, int32_t s0_b, int32_t s1_b);
#else
// this includes a C implmenentation of the above assembly components.
#include "masked_components.c"
#endif

// places random values on the heap for asm components
uint32_t* pre_rng(size_t n, uint32_t* addr) {
	for (int i = 0; i < n; ++i) {
		addr--;
		*addr = rand();
	}
	return addr;
}

#ifndef RNG_OFF
// mask the secrets before computation, 2 share xor masking
// bits are packed MSB first, so if we want to mask a single bit, 0b1(0..).
void mask_n_bits(int32_t n, int32_t* in, int32_t* s0_out, int32_t* s1_out) {
	int cur_rand;
	for (int i = 0; i < n/32; ++i) {
		cur_rand = rand();
		s0_out[i] = in[i] ^ cur_rand;
		s1_out[i] = cur_rand;
	}
	// n%32 bits remain to be masked
	if (n%32 != 0) {
		cur_rand = rand();
		// we have n%32 bits left to mask, so we shift out 32-n%32 bits from our random value
		// e.g. if 1 bit remains to be masked, we shift 0brrrr to get 0br000..0
		cur_rand = cur_rand << (32-(n%32));
		s0_out[n/32] = in[n/32] ^ cur_rand;
		s1_out[n/32] = cur_rand;
	}
}
#else
// unmasked version of mask_n_bits
void mask_n_bits(int32_t n, int32_t* in, int32_t* s0_out, int32_t* s1_out) {
	for (int i = 0; i < n/32; ++i) {
		s0_out[i] = in[i];
		s1_out[i] = 0;
	}
	if (n%32 != 0) s0_out[n/32] = in[n/32], s1_out[n/32] = 0;
}
#endif

// if there is only one index, the operation is not masked, but the result is also not interesting.
int masked_index_maxval(int n, int32_t *s0_vals, int32_t *s1_vals) {
	int32_t s0_max_index = 0, s1_max_index = 0, s0_max_val = s0_vals[0], s1_max_val = s1_vals[0], s0_sext_MSB, s1_sext_MSB, s0_diff, s1_diff, s0_and, s1_and;
	uint64_t diff, and64;
#if !(defined(RNG_OFF) || defined(NO_ASM)) // only pre-generate RNG when using the asm components with RNG on.
	uint32_t* addr = pre_rng((ADD_RWORDS + 4) * n, VOLATILE_END_ADDR);
	asm volatile ( // save to s6 the base addr where random numbers are stored.
			"mv s6, %[tAddr]\n"
			: : [tAddr] "r" (addr) : "s6");
#endif
	for (int i = 1; i < n; ++i) {
		// value is greater than the previous if their difference is negative (x1 - x2 < 0 -> x2 > x1).
		// so we calculate (masked) x1 + ~x2 + 1, note the +1 gets fed in as carry-in
		// also note ~x2 = ~s0_x2 ^ s1_x2; so we flip only one share.

		// this function is the masked version of these two operations:
		// max_index = (old_max_index & sext(MSB(max_conf - conf[i]))) ^ (i & sext(MSB(max_conf - conf(i)))
		// max_confidence = (old_max_confidence & sext(MSB(max_conf - conf[i]))) ^ (conf[i] & sext(MSB(max_conf - conf[i])))

		diff = m_add_32b(s0_max_val, s1_max_val, ~s0_vals[i], s1_vals[i], 0, 1);
		s0_diff = diff >> 32; s1_diff = diff;
		s0_sext_MSB = (s0_diff < 0) ? -1 : 0;
		s1_sext_MSB = (s1_diff < 0) ? -1 : 0;

		// note: reusing s0_diff, s1_diff as temporaries for the result of one masked AND
		diff = m_and(s0_max_index, s1_max_index, ~s0_sext_MSB, s1_sext_MSB);
		s0_diff = diff >> 32; s1_diff = diff;
		// we pass unmasked i ^ 0 as s0,s1, we don't care about leaking i, but after the the masked AND it will be randomised
		and64 = m_and(i, 0, s0_sext_MSB, s1_sext_MSB);
		s0_and = and64 >> 32; s1_and = and64;
		s0_max_index = s0_diff ^ s0_and;
		s1_max_index = s1_diff ^ s1_and;
		// same thing for the confidence score
		diff = m_and(s0_max_val, s1_max_val, ~s0_sext_MSB, s1_sext_MSB);
		s0_diff = diff >> 32; s1_diff = diff;
		and64 = m_and(s0_vals[i], s1_vals[i], s0_sext_MSB, s1_sext_MSB);
		s0_and = and64 >> 32; s1_and = and64;
		s0_max_val = s0_diff ^ s0_and;
		s1_max_val = s1_diff ^ s1_and;
	}
	return s0_max_index ^ s1_max_index;
	// we unmask the index and return it, it could easily be adapted to output both shares
	// however, in this case, this is the final result of our computation
	// and it is fine to unmask it here.
	// maybe more idiomatic to unmask it at the end in main?
}

// we calculate the dot product of vec(inputs) * vec(l1_weights), and activate if the dot prod >= 0.
void input_layer(uint8_t *input_pixels, uint32_t *weights, int32_t *s0_weights, int32_t *s1_weights, uint32_t *s0_l1_activations, uint32_t *s1_l1_activations, size_t nodes, size_t prev_nodes, size_t prev_nodes_word) { // @suppress("Type cannot be resolved")
	uint32_t* l1_weights_curr;
	for (int i = 0; i < nodes; ++i) { // loop over all nodes in current layer
		int32_t s0_sext_w, s1_sext_w, s0_curr_in, s1_curr_in, s0_curr_sum = 0, s1_curr_sum = 0;
		int64_t curr_sum;
		bool s0_w, s1_w;
		// the weights for every node is a long array, we want only the slice relevant to curr node
		l1_weights_curr = weights + i * prev_nodes_word;
		mask_n_bits(prev_nodes, l1_weights_curr, s0_weights, s1_weights);

		for (int j = 0; j < prev_nodes; ++j) { // loop over all nodes in prev layer
#if !(defined(RNG_OFF) || defined(NO_ASM)) // only pre-generate RNG when using the asm components with RNG on.
			if (j % MAX_ADDITIONS == 0) {
				uint32_t* addr = pre_rng(ADD_RWORDS * MAX_ADDITIONS, VOLATILE_END_ADDR);
				asm volatile ( // save to s6 the base addr where random numbers are stored, propagates across calls
						"mv s6, %[tAddr]\n"
						: : [tAddr] "r" (addr) : "s6");
			}
#endif
			// grab (MSB) masked weights out of the packed integer format, and shift to get ready for another round. consumes the masked weights!
			s0_w = s0_weights[j/32] < 0; s0_weights[j/32] <<= 1;
			s1_w = s1_weights[j/32] < 0; s1_weights[j/32] <<= 1;
			// sign extend the weights S.T. they can be XORd and ADDed. a positive weight (1) becomes -1 (0b111..11), otherwise 0.
			s0_sext_w = ~((int)s0_w - 1); // sign extend a boolean without branching
			s1_sext_w = ~((int)s1_w - 1);
			s0_curr_in = input_pixels[j]; // zero extend input char to 32bit
			s0_curr_in = ~s0_curr_in; // negate input I to get -I-1
 			s0_curr_in ^= s0_sext_w; // ^= with sign extended weight to either negate back to I if w=1.
			// feed the dot product (s0_curr_in, and s1_curr_in) to the adder.
			// notice that s1_curr_in is sext(s1_w).
 			// -I = ~I + 1, so ~(I ^ sext(w)) + !w -> if w = 0-> result=-I, if w = 1 -> result = I
 			// since +!w is only adding a single bit, we add it as carry-in to our adder
 			curr_sum = m_add_32b(s0_curr_sum, s1_curr_sum, s0_curr_in, s1_sext_w, !s0_w, s1_w);
			s0_curr_sum = curr_sum >> 32; s1_curr_sum = curr_sum;
		}
		// bits are packed MSB first
		// (negated most significant bit shifted into place : 1 or 0 from >=0) (aka activate if the sum is positive).
		s0_l1_activations[i/32] ^= (s0_curr_sum >= 0) << (31 - (i % 32));
		// (MSB (not flipped)) aka activate if sum is negative, part of the masking scheme, one share is flipped and the other is not.
		s1_l1_activations[i/32] ^= (s1_curr_sum < 0) << (31 - (i % 32));
		// the reason we flip one share: !a -> ~s0_a ^ s1_a where s0_a^s1_a = a.
	}
}

// The operation in this layer is optimised: (activation A) xnor (weight W) <=> ~A xor W; 32 multiplications happen at once in a single xnor
// The result of a single bit xnor means either +1 (1) or -1 (0) for our dot product, we count the set bits of all our xnors with a popcount.
// we calculate the dot product of vec(prev_layer) * vec(weights_curr_node), and activate if the dot prod >= 0.
// We store the activations packed in 32 bit format for easy xnoring on later layers.
void hidden_layer(uint32_t* weights, uint32_t* s0_weights, uint32_t* s1_weights, uint32_t* s0_activations_in, uint32_t* s1_activations_in, uint32_t* s0_activations_out, uint32_t* s1_activations_out, size_t node_count, size_t prev_node_count_word, size_t prev_node_count, bool output_layer) {
	asm volatile ( // these constants are propagated across every call to m_popcount.
			"li s1, 0x55555555\n"
			"li s2, 0x33333333\n"
			"li s3, 0x07070707\n"
			"li s4, 0x000F000F\n"
			"li s5, 0x0000001F\n"
			: : : "s1", "s2", "s3", "s4", "s5");
	uint32_t* curr_node_weights_ptr;
	for (int i = 0; i < node_count; ++i) { // loop over each node of the current layer
		int32_t s0_dot_prod = 0, s1_dot_prod = 0, s0_xnor, s1_xnor, s0_popcount, s1_popcount;
		int64_t popcount, dot_prod;
		// weights is a 2d arary
		curr_node_weights_ptr = weights + i * prev_node_count_word;
		mask_n_bits(prev_node_count, curr_node_weights_ptr, s0_weights, s1_weights);

#if !(defined(RNG_OFF) || defined(NO_ASM)) // only pre-generate RNG when using the asm components with RNG on.
		uint32_t* addr = pre_rng((POPCOUNT_RWORDS + ADD_RWORDS) * prev_node_count_word + ADD_RWORDS, VOLATILE_END_ADDR);
		asm volatile ( // save to s6 the base addr where random numbers are stored, propagated across calls
				"mv s6, %[tAddr]\n"
				: : [tAddr] "r" (addr) : "s6");
#endif
		for (int j = 0; j < prev_node_count_word; ++j) { // loop over the nodes of the previous layer, 32 bits at a time
			s0_xnor = ~s0_activations_in[j] ^ s0_weights[j];
			s1_xnor = s1_activations_in[j] ^ s1_weights[j];
			popcount = m_popcount(s0_xnor, s1_xnor);
			s0_popcount = popcount >> 32; s1_popcount = popcount;
			dot_prod = m_add_32b(s0_popcount, s1_popcount, s0_dot_prod, s1_dot_prod, 0, 0); // maybe can get away with adding fewer bits.
			s0_dot_prod = dot_prod >> 32; s1_dot_prod = dot_prod;
		}
		// final dot product = sum(popcounts) x 2 - prev_node_count
		dot_prod = m_add_32b(s0_dot_prod << 1, s1_dot_prod << 1, -prev_node_count, 0, 0, 0);
		s0_dot_prod = dot_prod >> 32; s1_dot_prod = dot_prod;
		if (!output_layer) {
			// (negated most significant bit shifted into place : 1 or 0 from >=0) (aka activate if the sum is positive).
			s0_activations_out[i/32] ^= (s0_dot_prod >= 0) << (31 - (i % 32));
			// (MSB (not flipped)) aka activate if sum is negative, part of the masking scheme, one share is flipped and the other is not.
			s1_activations_out[i/32] ^= (s1_dot_prod < 0) << (31 - (i % 32));
		} else {
			// if we are calculating the output layer, we do not calculate activations values, rather confidence scores.
			s0_activations_out[i] = s0_dot_prod;
			s1_activations_out[i] = s1_dot_prod;
		}
	}
}

uint64_t perf_testing() {
	struct metal_cpu *cpu = metal_cpu_get(metal_cpu_get_current_hartid());/* Get CPU device handle. */
	metal_hpm_init(cpu);

	uint64_t cycles = metal_hpm_read_counter(cpu, METAL_HPM_CYCLE);

	// your testable code here

	uint64_t cycles2 = metal_hpm_read_counter(cpu, METAL_HPM_CYCLE);
	int32_t total = cycles2 - cycles;
	int32_t total_upper = (cycles2 - cycles) >> 32;
	printf("clock cycles: %i, lower: %i, hex: 0x%08X%08X\n", total_upper, total, total_upper, total);
	fflush(stdout);
	return cycles2 - cycles;
}

// int metal_hpm_clear_counter(struct metal_cpu *gcpu, metal_hpm_counter counter);
#ifndef TESTING
int main () {
	struct metal_cpu *cpu = metal_cpu_get(metal_cpu_get_current_hartid());/* Get CPU device handle. */
	metal_hpm_init(cpu);
	srand(metal_hpm_read_counter(cpu, METAL_HPM_CYCLE)); // seed rng with cpu cycle count, fairly unpredictable
	// srand(0);

	// ## INPUT LAYER
	uint32_t s0_l1_activations[L1_NODES_WORD] = {0};
	uint32_t s1_l1_activations[L1_NODES_WORD] = {0};
	{
		uint8_t input_pixels[INPUT_PIXELS]; // @suppress("Type cannot be resolved")
		memcpy(input_pixels, input_pixel_bytes, INPUT_PIXELS);
		uint32_t* l1_weights = (uint32_t *) l1_weight_bytes; // each node of the input layer has a row with its weights to the input nodes (pixels).
		uint32_t s0_weights_placeholder[INPUT_PIXELS_WORD];
		uint32_t s1_weights_placeholder[INPUT_PIXELS_WORD];

		input_layer(input_pixels, l1_weights, s0_weights_placeholder, s1_weights_placeholder, s0_l1_activations, s1_l1_activations, L1_NODES, INPUT_PIXELS, INPUT_PIXELS_WORD);
	}
	 // ## LAYER 2.
	 uint32_t s0_l2_activations[L2_NODES_WORD] = {0};
	 uint32_t s1_l2_activations[L2_NODES_WORD] = {0};

	{
		uint32_t* l2_weights = (uint32_t *) l2_weight_bytes; // [L2_NODES][L1_NODES_WORD]; // each node of L2 has a row with weights to L1.
		uint32_t s0_l2_weights[L1_NODES_WORD];
		uint32_t s1_l2_weights[L1_NODES_WORD];

		hidden_layer(l2_weights, s0_l2_weights, s1_l2_weights, s0_l1_activations, s1_l1_activations, s0_l2_activations, s1_l2_activations, L2_NODES, L1_NODES_WORD, L1_NODES, false);
	}

	// ## LAYER 3.
	uint32_t s0_l3_activations[L3_NODES_WORD] = {0};
	uint32_t s1_l3_activations[L3_NODES_WORD] = {0};
	{
		uint32_t* l3_weights = (uint32_t *) l3_weight_bytes; //[L3_NODES][L2_NODES_WORD]; // each node of L2 has a row with weights to L1.
		uint32_t s0_l3_weights[L2_NODES_WORD];
		uint32_t s1_l3_weights[L2_NODES_WORD];

		hidden_layer(l3_weights, s0_l3_weights, s1_l3_weights, s0_l2_activations, s1_l2_activations, s0_l3_activations, s1_l3_activations, L3_NODES, L2_NODES_WORD, L2_NODES, false);
	}
	// ## LAYER 4 / OUTPUT LAYER.
  	int32_t s0_out_conf[OUTPUT_NODES] = {0};
	int32_t s1_out_conf[OUTPUT_NODES] = {0};
	{
		int32_t* out_weights = (uint32_t *) out_weight_bytes; //[OUTPUT_NODES][L3_NODES_WORD];
		int32_t s0_out_weights[L3_NODES_WORD];
		int32_t s1_out_weights[L3_NODES_WORD];

		hidden_layer(out_weights, s0_out_weights, s1_out_weights, s0_l3_activations, s1_l3_activations, s0_out_conf, s1_out_conf, OUTPUT_NODES, L3_NODES_WORD, L3_NODES, true);
	}
	// ## OUTPUT FUNCTION
	int res = masked_index_maxval(OUTPUT_NODES, s0_out_conf, s1_out_conf);

	printf("result %i\n", res); fflush(stdout);

	return 0;
}
#endif
