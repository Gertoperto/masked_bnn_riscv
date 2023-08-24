.section .text
// #define NO_ASM
// #define RNG_OFF
#ifndef NO_ASM
.globl m_add_32b
.type m_add_32b, @function
#endif
// take 2 masked ints and 1 masked carry-in bit (6 args total) -> return 2 masked ints and idc abt c_out



// a0 a1 a2 a3 a4    a5
// a0 a1 b0 b1 c_in0 c_in1
// IMPORTANT only the LSB of c_in can be set, the rest should be zero

// temporaries:
// a6: (a&b)_0 (and r1 ^ r2 ^ r3, later).
// a7: (a&b)_1
// t0: mask bit
// t1-t2: calculations a&c, and shared carry_out temps
// t3-t4: calculations b&c
// t5-t6: random values (r)
m_add_32b:
// mv tp, ra // save return addr in frame pointer (idk if tp is important?) (no more sub calls)
// *s0_out = s1_a & s1_b ^ (s1_a & s0_b ^ (s0_a & s1_b ^ (s0_a & s0_b ^ r)));
// *s1_out = r; // if the calling function holds on to r, we can get rid of s1_out.
#ifdef RNG_OFF
mv a6, x0 // r = 0
mv t5, x0
mv t6, x0
#else
lw a6, 0(s6) // load r
lw t5, 4(s6) // load r for a&c, b&c carry calculations ahead of time
lw t6, 8(s6)
addi s6, s6, 12
#endif

m_and_ab: // (r) ab1 -> a6, ab0 -> a7
and t0, a0, a2 // s0_a & s0_b
and t1, a0, a3 // s0_a & s1_b
xor t0, t0, a6 // s0_a & s0_b ^ r
xor t0, t1, t0 // s0_a & s1_b ^ (s0_a & s0_b ^ r)
and t1, a1, a2 // s1_a & s0_b
xor t0, t1, t0 // s1_a & s0_b ^ (s0_a & s1_b ^ (s0_a & s0_b ^ r))
and t1, a1, a3 // s1_a & s1_b
xor a7, t1, t0 // s1_a & s1_b ^ (s1_a & s0_b ^ (s0_a & s1_b ^ (s0_a & s0_b ^ r)))

/*
and a6, a0, a2 // s0_a & s0_b as a6
xor a6, a6, x0 // a6 ^ r as a6 here x0 needs to be replaced with R, we have no source of randomness yet
and a7, a0, a3 // s0_a & s1_b : a7 next temporary
xor a6, a6, a7 // (s0_a & s1_b ^ (s0_a & s0_b ^ r))
and a7, a1, a2 // s1_a & s0_b
xor a6, a6, a7 // (s1_a & s0_b ^ (s0_a & s1_b ^ (s0_a & s0_b ^ r)))
and a7, a1, a3 // s1_a & s1_b
xor a6, a6, a7 // s1_a & s1_b ^ (s1_a & s0_b ^ (s0_a & s1_b ^ (s0_a & s0_b ^ r)))
mv a7, x0 // R -> a7 (FIXME) */

pre_carry_loop:
li t0, 1 // mask whichever bit of the carry we calculate
// a6 is only used again to calculate the second share of the carry:
// (a&b)1 ^ (b&c)1 ^ (c&a)1, so we can precalculate it and store it in a6
// if rng_off then a6 will just hold zero
xor a6, a6, t5
xor a6, a6, t6

carry_loop:
m_and_ac: // t1: ac0, t5: ac1 (r)
and t1, a0, a4 // s0_a & s0_c as t1
xor t1, t1, t5 // t1 ^ r as t1 here x0 needs to be replaced with R, we have no source of randomness yet
and t2, a0, a5 // s0_a & s1_c : t2 next temporary
xor t1, t1, t2 // (s0_a & s1_c ^ (s0_a & s0_c ^ r))
and t2, a1, a4 // s1_a & s0_c
xor t1, t1, t2 // (s1_a & s0_c ^ (s0_a & s1_c ^ (s0_a & s0_c ^ r)))
and t2, a1, a5 // s1_a & s1_c
xor t1, t1, t2 // s1_a & s1_c ^ (s1_a & s0_c ^ (s0_a & s1_c ^ (s0_a & s0_c ^ r)))

m_and_bc: // t3: bc0, t6: bc1 (r)
and t3, a2, a4 // s0_b & s0_c as t1
xor t3, t3, t6 // t1 ^ r as t1 here x0 needs to be replaced with R, we have no source of randomness yet
and t4, a2, a5 // s0_b & s1_c : t2 next temporary
xor t3, t3, t4 // (s0_b & s1_c ^ (s0_b & s0_c ^ r))
and t4, a3, a4 // s1_b & s0_c
xor t3, t3, t4 // (s1_b & s0_c ^ (s0_b & s1_c ^ (s0_b & s0_c ^ r)))
and t4, a3, a5 // s1_b & s1_c
xor t3, t3, t4 // s1_b & s1_c ^ (s1_b & s0_c ^ (s0_b & s1_c ^ (s0_b & s0_c ^ r)))

carry:
// t1: (a&b)0 ^ (b&c)0 ^ (c&a)0
xor t1, t1, t3
xor t1, t1, a7
// a6 holds: (a&b)1 ^ (b&c)1 ^ (c&a)1

// bitmask and shift carry_out into the carry_in bit of the next round
and t1, t1, t0
and t2, a6, t0 // no need to move a6 into a5, can just use a6 all the way thru
slli t1, t1, 1
slli t2, t2, 1
xor a4, a4, t1
xor a5, a5, t2
// (the last carryout is simply shifted off and lost, wont be needing it anyway).

check_done:
slli t0, t0, 1
bnez t0, carry_loop

done:
// calculate the sum bits and return
xor a0, a0, a2
xor a0, a0, a4
xor a1, a1, a3
xor a1, a1, a5
ret





