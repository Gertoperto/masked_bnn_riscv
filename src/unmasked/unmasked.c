#include <stdbool.h>
#include <stdint.h>
#include <limits.h>
#include "static.c"
#include <stdio.h>
#include <string.h>

int main() {
	// input layer : 784 nodes (28x28 8-bit pixels) : 
	// hidden layers 3x : 1010 nodes + activations :
	// the first hidden layer each node computes the sum of all input nodes times their respective weights 
	// so +1 / -1 * 8-bit int -> each node has 784 weights connected to the input layer
	// the other two hidden layers each compute the sum of weight * activation of the previous layer
	// activation -> sum is positve? activate. otherwise no -> return the sign bit

	// move input pixels from ROM to the stack (RAM)
	uint8_t input_pixels[INPUT_PIXELS]; // @suppress("Type cannot be resolved")
	memcpy(input_pixels, input_pixel_bytes, INPUT_PIXELS);

	uint32_t *weights_l1 = (uint32_t*) l1_weight_bytes; 
	int activations_l1[L1_NODES_WORD] = {0};
	for (int i = 0; i < L1_NODES; ++i) {
		int32_t curr_node_weights[INPUT_PIXELS_WORD]; // load weights for curr layer from ROM
		memcpy(curr_node_weights, weights_l1 + i * INPUT_PIXELS_WORD, INPUT_PIXELS_WORD * 4);
		int node_dot_product = 0; // calculate dot product
		for (int j = 0; j < INPUT_PIXELS; ++j) {
			bool w = curr_node_weights[j/32] < 0; // store most significant bit as `w`
			curr_node_weights[j/32] <<= 1; // shift the weights over for next round
			node_dot_product += (w) ? input_pixels[j] : -1 * input_pixels[j];
		}
		// shift activation into place
		activations_l1[i/32] ^= (node_dot_product >= 0) << (31-(i%32)); 
	}

	uint32_t *weights_l2 = (uint32_t*) l2_weight_bytes;
	int activations_l2[L2_NODES_WORD] = {0};
	for (int i = 0; i < L2_NODES; ++i) {
		int32_t *curr_node_weights = weights_l2 + i * L1_NODES_WORD;
		int32_t node_dot_product = 0;
		for (int j = 0; j < L1_NODES_WORD; ++j) {
			uint32_t xnor = ~activations_l1[j] ^ curr_node_weights[j];
			node_dot_product += __builtin_popcount(xnor); // popcount(xnor). 
		}
		node_dot_product *= 2;
		node_dot_product -= L1_NODES;
		activations_l2[i/32] ^= (node_dot_product >= 0) << (31-(i%32)); // maybe > 0 but probably 0 should be counted as positve/activate 
	}	

	uint32_t *weights_l3 = (uint32_t*) l3_weight_bytes;
	int activations_l3[L3_NODES_WORD] = {0};
	for (int i = 0; i < L3_NODES; ++i) {
		int32_t *curr_node_weights = weights_l3 + i * L2_NODES_WORD;
		int node_dot_product = 0;
		for (int j = 0; j < L2_NODES_WORD; ++j) {
			node_dot_product += __builtin_popcount(~activations_l2[j] ^ curr_node_weights[j]); // popcount(xnor). 
		}
		node_dot_product *= 2;
		node_dot_product -= L1_NODES;
		activations_l3[i/32] ^= (node_dot_product >= 0) << (31-(i%32)); // maybe > 0 but probably 0 should be counted as positve/activate 
	}

	uint32_t *weights_output = (uint32_t*) out_weight_bytes;
	int confidence_scores[OUTPUT_NODES] = {0}; // there are 10 output nodes, they calculate the same dot product but do not create an activation
	for (int i = 0; i < OUTPUT_NODES; ++i) {
		int32_t *curr_node_weights = weights_output + i * L3_NODES_WORD;
		int node_dot_product = 0;
		for (int j = 0; j < L3_NODES_WORD; ++j) {
			node_dot_product += __builtin_popcount(~activations_l3[j] ^ curr_node_weights[j]); // popcount(xnor). 
		}
		node_dot_product *= 2;
		node_dot_product -= L1_NODES;
		confidence_scores[i] = node_dot_product; 
	}

	int maxc = INT_MIN;
	int index = -1;
	for (int i = 0; i < OUTPUT_NODES; ++i) {
		printf("score i:%i : %i\n", i, confidence_scores[i]);
		if (confidence_scores[i] >= maxc) {
			maxc = confidence_scores[i]; index = i;
		}

	}
	printf("output: %i\n", index);
	fflush(stdout);

	return 0;
	// side node: can we precompute the entire sum of input pixels? then we only need to remove 2x the negative weight
}
