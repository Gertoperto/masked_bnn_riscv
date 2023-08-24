.section .text
// #define NO_ASM
// #define RNG_OFF
#ifndef NO_ASM
.globl m_and
.type m_and, @function
#endif



// *s0_out = s1_a & s1_b ^ (s1_a & s0_b ^ (s0_a & s1_b ^ (s0_a & s0_b ^ r)));
// *s1_out = r; // if the calling function holds on to r, we can get rid of s1_out.
// var: a0 a1 a2 a3
// reg: a0 a1 b0 b1
m_and:
and t0, a0, a2 // s0_a & s0_b
and t1, a0, a3 // s0_a & s1_b
#ifdef RNG_OFF
mv a0, x0
#else
lw a0, 0(s6)
addi s6, s6, 4
#endif
xor t0, t0, a0 // s0_a & s0_b ^ r
xor t0, t1, t0 // s0_a & s1_b ^ (s0_a & s0_b ^ r)
and t1, a1, a2 // s1_a & s0_b
xor t0, t1, t0 // s1_a & s0_b ^ (s0_a & s1_b ^ (s0_a & s0_b ^ r))
and t1, a1, a3 // s1_a & s1_b
xor a1, t1, t0 // s1_a & s1_b ^ (s1_a & s0_b ^ (s0_a & s1_b ^ (s0_a & s0_b ^ r)))
ret
