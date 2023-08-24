#ifndef STATIC_H
#define STATIC_H

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

#include <stdint.h>

extern const uint8_t input_pixel_bytes[]; // @suppress("Type cannot be resolved")
extern const uint8_t l1_weight_bytes[]; // @suppress("Type cannot be resolved")
extern const uint8_t l2_weight_bytes[]; // @suppress("Type cannot be resolved")
extern const uint8_t l3_weight_bytes[]; // @suppress("Type cannot be resolved")
extern const uint8_t out_weight_bytes[]; // @suppress("Type cannot be resolved")

#endif
