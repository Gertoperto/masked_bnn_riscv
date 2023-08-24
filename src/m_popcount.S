.section .text
// #define NO_ASM
// #define RNG_OFF
#ifndef NO_ASM
.globl m_popcount
.type m_popcount, @function
#endif

#ifdef NO_CONST_PROP
// we do not use this, we load the consts from our main C file
load_consts: // only needs to happen once, other popcount calls can reuse it
.equ MASK0101, 0x55555555
.equ MASK0011, 0x33333333
.equ MASK0111, 0x07070707
.equ MASKBYTE, 0x000F000F
.equ MASKHALF, 0x0000001F
li s1, MASK0101
li s2, MASK0011
li s3, MASK0111
li s4, MASKBYTE
li s5, MASKHALF
#endif

// a0 a1 -> a
// shared input 32b
// TEMPORARIES:
// a2, a3 -> b
// a4, a5 -> c
// t0, t1 -> a&b
// t3, t5 -> random bits
// t2, t4 -> temporaries for calculations
// t6 -> bitmask
m_popcount:
mv tp, ra
sums_1bit:
// mask even and odd bits, shift to line up
srli a2, a0, 1
srli a3, a1, 1
and a0, a0, s1
and a1, a1, s1
and a2, a2, s1
and a3, a3, s1
// (a0,a1) (a2,a3) have odd bits, and even bits, stored in the even bit positions (respectively)
mv t6, s1 // use s1 as bitmask for precarry
jal pre_carry
jal sum_calc

// 2bit accumulators get added together to form 3 bit accumulators
// so we need to do two rounds of carry calculation
sums_2bit:
// mask even and odd bits, shift to line up
srli a2, a0, 2
srli a3, a1, 2
and a0, a0, s2
and a1, a1, s2
and a2, a2, s2
and a3, a3, s2
// now our registers look like this 0b00aa00aa -> a0,a1
//                                  0b00bb00bb -> a2,a3
and t6, s1, s2 // make new constant for 0b0001 (repeating) to mask multiple carry bits
jal pre_carry
jal calc_carry
jal sum_calc

// output looks like this: 0sss0sss0sss
// (3 sum bits and one zero, repeating)
sums_4bit:
// mask even and odd bits, shift to line up
srli a2, a0, 4
srli a3, a1, 4
and a0, a0, s3
and a1, a1, s3
and a2, a2, s3
and a3, a3, s3

and t6, s1, s2 // calculate bitmask -> 0b0000000100000001
and t6, t6, s3
jal pre_carry
jal calc_carry
jal calc_carry
jal sum_calc
// output looks like this: 0000ssss0000ssss

sums_8bit:
srli a2, a0, 8
srli a3, a1, 8
and a0, a0, s4
and a1, a1, s4
and a2, a2, s4
and a3, a3, s4

srli t6, t6, 3 // shift old bitmask back
and t6, t6, s4
jal pre_carry
jal calc_carry
jal calc_carry
jal calc_carry
jal sum_calc

sums_16bit:
srli a2, a0, 16
srli a3, a1, 16
and a0, a0, s5
and a1, a1, s5
and a2, a2, s5
and a3, a3, s5

li t6, 1
jal pre_carry
jal calc_carry
jal calc_carry
jal calc_carry
jal calc_carry
jal sum_calc

mv ra, tp
ret // return to caller with the sum in a0, a1 (6 bit sum).

// a0 a1 a3 a4 t6      outputs: a4 a5 - t0   t1
// a0 a1 b0 b1 bitmask          c0 c1 - a&b0 a&b1
// purpose is to calculate the first carry bit, and prepare a&b
pre_carry:
#ifndef RNG_OFF
lw t1, 0(s6)
#else
mv t1, x0
#endif

m_and_ab: // ab0 ->t0, ab1->t1 (r)
// s1_a & s1_b ^ (s1_a & s0_b ^ (s0_a & s1_b ^ (s0_a & s0_b ^ r)));
// s0_a & s0_b ^ r
and t0, a0, a2
xor t0, t0, t1 // x0 fake source of randomness, fixme
// s0_a & s1_b
and t2, a0, a3
xor t0, t0, t2
// s1_a & s0_b
and t2, a1, a2
xor t0, t0, t2
// s1_a & s1_b
and t2, a1, a3
xor t0, t0, t2

// the first bit has no carry in, so the first carry out is simply a&b
mv a4, x0 // clear temps to store calculated carry
mv a5, x0
// mask a&b to get only the carry bit of the first round
and t2, t0, t6
and t3, t1, t6
slli t2, t2, 1
slli t3, t3, 1
// xor the carry bit into place, a4,a5 keep all carries
xor a4, a4, t2
xor a5, a5, t3

slli t6, t6, 1

#ifndef RNG_OFF
lw t3, 4(s6)
lw t5, 8(s6)
addi s6, s6, 12
#else
mv t3, x0
mv t5, x0
#endif

xor t1, t1, t3 // r1 r2 and r3 do not change for every call of calc_carry
xor t1, t1, t5 // if there are enough temporaries, maybe move them to pre_carry

ret

// reg: a0 a1 a2 a3 a4 a5 t0   t1     t6      t3 t5
// var: a0 a1 b0 b1 c0 c1 a&b0 a&b1=r bitmask r  r
calc_carry:
m_and_ac: // result in ac0->t2, ac1->t3 (r)
// s1_a & s1_c ^ (s1_a & s0_c ^ (s0_a & s1_c ^ (s0_a & s0_c ^ r)));
// s0_a & s0_c ^ r
and t2, a0, a4
xor t2, t2, t3 // x0 fake source of randomness, fixme
// s0_a & s1_c
and t4, a0, a5
xor t2, t2, t4
// s1_a & s0_c
and t4, a1, a4
xor t2, t2, t4
// s1_a & s1_c
and t4, a1, a5
xor t2, t2, t4

m_and_bc: // result in bc0->t4, bc1->t5 (r)
// s1_b & s1_c ^ (s1_b & s0_c ^ (s0_b & s1_c ^ (s0_b & s0_c ^ r)));
// s0_b & s0_c ^ r
and t4, a2, a4
xor t4, t4, t5 // x0 fake source of randomness, fixme
xor t2, t2, t4
// s0_b & s1_c
and t4, a2, a5
xor t2, t2, t4
// s1_b & s0_c
and t4, a3, a4
xor t2, t2, t4
// s1_b & s1_c
and t4, a3, a5
xor t2, t2, t4

// ab ^ (bc ^ ac) -> t2,t3
xor t2, t2, t0
// xor t2, t2, t4 we xor'd all the components of b&c with a&c already to save a temporary
// t1 currently holds (ab^bc^ac)1 = (r^r^r).

// mask carry bit and shift into place
and t2, t2, t6
and t4, t1, t6
slli t2, t2, 1
slli t4, t4, 1
xor a4, a4, t2
xor a5, a5, t4

slli t6, t6, 1

ret

sum_calc:
// calculate the sum bits a^b^c -> a0, a1
xor a0, a0, a2
xor a0, a0, a4

xor a1, a1, a3
xor a1, a1, a5
ret
