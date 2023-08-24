#include "static.h"

const uint8_t input_pixel_bytes[INPUT_PIXELS] = { // @suppress("Type cannot be resolved")
#include "static_data/input_pixels_hexstr.txt"
};
const uint8_t l1_weight_bytes[INPUT_PIXELS_WORD * 4 * L1_NODES] = { // @suppress("Type cannot be resolved")
#include "static_data/weights_l1_hexstr.txt"
};
const uint8_t l2_weight_bytes[L1_NODES_WORD * 4 * L2_NODES] = { // @suppress("Type cannot be resolved")
#include "static_data/weights_l2_hexstr.txt"
};
const uint8_t l3_weight_bytes[L2_NODES_WORD * 4 * L3_NODES] = { // @suppress("Type cannot be resolved")
#include "static_data/weights_l3_hexstr.txt"
};
const uint8_t out_weight_bytes[L3_NODES_WORD * 4 * OUTPUT_NODES] = { // @suppress("Type cannot be resolved")
#include "static_data/weights_out_hexstr.txt"
};
