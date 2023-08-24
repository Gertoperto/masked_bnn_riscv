# python 3.9
from random import randbytes, getrandbits

num_bytes = 28*28//8 * 1010
input_pixels = 28*28
l1_nodes = 1010
l2_nodes = 1010
l3_nodes = 1010
out_nodes = 10

# print(100)
# print(b",".join([x for x in randbytes(num_bytes)]))
# for byte in randbytes(num_bytes):
#     print(b"%b," % byte, end="")
# print(randbytes(num_bytes))
# print("\"", "".join('\\x{:02X}'.format(b) for b in randbytes(num_bytes)), "\"", end="", sep="")
# print(b2a_qp(randbytes(num_bytes)))

def gen(prev_nodes, curr_nodes):
    print("\"", end="", sep="")
    for _ in range(curr_nodes):
        for _ in range(prev_nodes // 32):
            print("\\x", "\\x".join('{:02X}'.format(b) for b in randbytes(4)), end="", sep="")
        if prev_nodes % 32:
            rand_int = getrandbits(prev_nodes % 32) << (32 - (prev_nodes % 32))
            print("\\x", "\\x".join("{:02X}".format(b) for b in rand_int.to_bytes(4, byteorder='little', signed=False)), end="", sep="")
    print("\"", end="", sep="")

# generate input bytes
# print("\"", "".join('\\x{:02X}'.format(b) for b in randbytes(input_pixels)), "\"", end="", sep="")

# gen(input_pixels, l1_nodes)

# gen(l1_nodes, l2_nodes)

# gen(l2_nodes, l3_nodes)

gen(l3_nodes, out_nodes)