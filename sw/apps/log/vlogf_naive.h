// Copyright 2024 ETH Zurich and University of Bologna.
// Licensed under the Apache License, Version 2.0, see LICENSE for details.
// SPDX-License-Identifier: Apache-2.0
//
// Luca Colagrande <colluca@iis.ee.ethz.ch>

static inline void vexpf_naive(float *a, float *b) {

const struct logf_data __logf_data = {
  .tab = {
  { 0x1.661ec79f8f3bep+0, -0x1.57bf7808caadep-2 },
  { 0x1.571ed4aaf883dp+0, -0x1.2bef0a7c06ddbp-2 },
  { 0x1.49539f0f010bp+0, -0x1.01eae7f513a67p-2 },
  { 0x1.3c995b0b80385p+0, -0x1.b31d8a68224e9p-3 },
  { 0x1.30d190c8864a5p+0, -0x1.6574f0ac07758p-3 },
  { 0x1.25e227b0b8eap+0, -0x1.1aa2bc79c81p-3 },
  { 0x1.1bb4a4a1a343fp+0, -0x1.a4e76ce8c0e5ep-4 },
  { 0x1.12358f08ae5bap+0, -0x1.1973c5a611cccp-4 },
  { 0x1.0953f419900a7p+0, -0x1.252f438e10c1ep-5 },
  { 0x1p+0, 0x0p+0 },
  { 0x1.e608cfd9a47acp-1, 0x1.aa5aa5df25984p-5 },
  { 0x1.ca4b31f026aap-1, 0x1.c5e53aa362eb4p-4 },
  { 0x1.b2036576afce6p-1, 0x1.526e57720db08p-3 },
  { 0x1.9c2d163a1aa2dp-1, 0x1.bc2860d22477p-3 },
  { 0x1.886e6037841edp-1, 0x1.1058bc8a07ee1p-2 },
  { 0x1.767dcf5534862p-1, 0x1.4043057b6ee09p-2 },
  },
  .ln2 = 0x1.62e42fefa39efp-1,
  .poly = {
  -0x1.00ea348b88334p-2, 0x1.5575b0be00b6ap-2, -0x1.ffffef20a4123p-2,
  }
};

    #define T __logf_data.tab
    #define A __logf_data.poly
    #define Ln2 __logf_data.ln2
    #define N (1 << LOGF_TABLE_BITS)
    #define OFF 0x3f330000


    const uint32_t LOGF_TABLE_BITS = 4;
    const uint32_t OFF = 0x3f330000;
    const uint32_t N = 1 << LOGF_TABLE_BITS;
    const double Ln2 = 0x1.62e42fefa39efp-1;
    const double A[3] = {-0x1.00ea348b88334p-2, 0x1.5575b0be00b6ap-2, -0x1.ffffef20a4123p-2};

    uint64_t ki, t;

    // Loop over samples (unrolled by 4)
    for (int i = 0; i < LEN / 4; i++) {
        asm volatile(
fcvt.s.d fa0, ft0
fmv.x.w a0, fa0   // ix = asuint (x)
mv a2, a0
srai a6, a2, 23 // k += (ix>>23)-127
add a0, a1, a6 // k += (ix>>23)-127
            "fcvt.d.s fa1, %[input]            \n" // xd = (double_t) x
            "fmul.d   fa3, %[InvLn2N], fa1     \n" // z = InvLn2N * xd
            "fadd.d   fa1, fa3, %[SHIFT]       \n" // kd = (double) (z + SHIFT)
            "fsd      fa1, 0(%[ki])            \n" // ki = asuint64 (kd)
            "lw       a0, 0(%[ki])             \n" // ki = asuint64 (kd)
            "andi     a1, a0, 0x1f             \n" // ki % N
            "slli     a1, a1, 0x3              \n" // T[ki % N]
            "add      a1, %[T], a1             \n" // T[ki % N]
            "lw       a2, 0(a1)                \n" // t = T[ki % N]
            "lw       a1, 4(a1)                \n" // t = T[ki % N]
            "slli     a0, a0, 0xf              \n" // ki << (52 - EXP2F_TABLE_BITS)
            "sw       a2, 0(%[t])              \n" // store lower 32b of t (unaffected)
            "add      a0, a0, a1               \n" // t += ki << (52 - EXP2F_TABLE_BITS)
            "sw       a0, 4(%[t])              \n" // store upper 32b of t
            "fsub.d   fa2, fa1, %[SHIFT]       \n" // kd -= SHIFT
            "fsub.d   fa3, fa3, fa2            \n" // r = z - kd
            "fmadd.d  fa2, %[C0], fa3, %[C1]   \n" // z = C[0] * r + C[1]
            "fld      fa0, 0(%[t])             \n" // s = asdouble (t)
            "fmadd.d  fa4, %[C2], fa3, %[C3]   \n" // y = C[2] * r + C[3]
            "fmul.d   fa1, fa3, fa3            \n" // r2 = r * r
            "fmadd.d  fa4, fa2, fa1, fa4       \n" // y = z * r2 + y
            "fmul.d   fa4, fa4, fa0            \n" // y = y * s
            "fcvt.s.d %[output], fa4           \n" // (float) y
            : [output] "=f" (b[i])
            : [input] "f" (a[i]), [InvLn2N] "f" (InvLn2N), [SHIFT] "f" (SHIFT),
              [C0] "f" (C[0]), [C1] "f" (C[1]), [C2] "f" (C[2]), [C3] "f" (C[3]),
              [ki] "r" (&ki), [t] "r" (&t), [T] "r" (T)
            : "memory", "a0", "a1", "a2", "fa0", "fa1", "fa2", "fa3", "fa4"
        );
    }
}
