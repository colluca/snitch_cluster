#define FP0_ASM_BODY                   \
    "fmul.d   fa3, %[InvLn2N], ft0 \n" \
    "fmul.d   ft3, %[InvLn2N], ft0 \n" \
    "fmul.d   ft4, %[InvLn2N], ft0 \n" \
    "fmul.d   ft5, %[InvLn2N], ft0 \n" \
    "fadd.d   ft1, fa3, %[SHIFT]   \n" \
    "fadd.d   ft1, ft3, %[SHIFT]   \n" \
    "fadd.d   ft1, ft4, %[SHIFT]   \n" \
    "fadd.d   ft1, ft5, %[SHIFT]   \n" \
    "fmv.d    ft2, fa3             \n" \
    "fmv.d    ft2, ft3             \n" \
    "fmv.d    ft2, ft4             \n" \
    "fmv.d    ft2, ft5             \n"

#define INT_ASM_BODY        \
    "lw   a0,  0(%[ki]) \n" \
    "lw   a3,  8(%[ki]) \n" \
    "lw   a4, 16(%[ki]) \n" \
    "lw   a5, 24(%[ki]) \n" \
    "andi a1, a0, 0x1f  \n" \
    "andi a6, a3, 0x1f  \n" \
    "andi a7, a4, 0x1f  \n" \
    "andi t0, a5, 0x1f  \n" \
    "slli a1, a1, 0x3   \n" \
    "slli a6, a6, 0x3   \n" \
    "slli a7, a7, 0x3   \n" \
    "slli t0, t0, 0x3   \n" \
    "add  a1, %[T], a1  \n" \
    "add  a6, %[T], a6  \n" \
    "add  a7, %[T], a7  \n" \
    "add  t0, %[T], t0  \n" \
    "lw   a2, 0(a1)     \n" \
    "lw   t1, 0(a6)     \n" \
    "lw   t2, 0(a7)     \n" \
    "lw   t3, 0(t0)     \n" \
    "lw   a1, 4(a1)     \n" \
    "lw   a6, 4(a6)     \n" \
    "lw   a7, 4(a7)     \n" \
    "lw   t0, 4(t0)     \n" \
    "slli a0, a0, 0xf   \n" \
    "slli a3, a3, 0xf   \n" \
    "slli a4, a4, 0xf   \n" \
    "slli a5, a5, 0xf   \n" \
    "sw   a2,  0(%[t])  \n" \
    "sw   t1,  8(%[t])  \n" \
    "sw   t2, 16(%[t])  \n" \
    "sw   t3, 24(%[t])  \n" \
    "add  a0, a0, a1    \n" \
    "add  a3, a3, a6    \n" \
    "add  a4, a4, a7    \n" \
    "add  a5, a5, t0    \n" \
    "sw   a0,  4(%[t])  \n" \
    "sw   a3, 12(%[t])  \n" \
    "sw   a4, 20(%[t])  \n" \
    "sw   a5, 28(%[t])  \n"

#define FP1_ASM_BODY                     \
    "fsub.d   fa2, ft0, %[SHIFT]     \n" \
    "fsub.d   ft6, ft0, %[SHIFT]     \n" \
    "fsub.d   ft7, ft0, %[SHIFT]     \n" \
    "fsub.d   ft8, ft0, %[SHIFT]     \n" \
    "fsub.d   fa3, ft0, fa2          \n" \
    "fsub.d   ft3, ft0, ft6          \n" \
    "fsub.d   ft4, ft0, ft7          \n" \
    "fsub.d   ft5, ft0, ft8          \n" \
    "fmadd.d  fa2, %[C0], fa3, %[C1] \n" \
    "fmadd.d  ft6, %[C0], ft3, %[C1] \n" \
    "fmadd.d  ft7, %[C0], ft4, %[C1] \n" \
    "fmadd.d  ft8, %[C0], ft5, %[C1] \n" \
    "fmadd.d  fa4, %[C2], fa3, %[C3] \n" \
    "fmadd.d  fs0, %[C2], ft3, %[C3] \n" \
    "fmadd.d  fs1, %[C2], ft4, %[C3] \n" \
    "fmadd.d  fs2, %[C2], ft5, %[C3] \n" \
    "fmul.d   fa1, fa3, fa3          \n" \
    "fmul.d   fa5, ft3, ft3          \n" \
    "fmul.d   fa6, ft4, ft4          \n" \
    "fmul.d   fa7, ft5, ft5          \n" \
    "fmadd.d  fa4, fa2, fa1, fa4     \n" \
    "fmadd.d  fs0, ft6, fa5, fs0     \n" \
    "fmadd.d  fs1, ft7, fa6, fs1     \n" \
    "fmadd.d  fs2, ft8, fa7, fs2     \n" \
    "fmul.d   ft2, fa4, ft1          \n" \
    "fmul.d   ft2, fs0, ft1          \n" \
    "fmul.d   ft2, fs1, ft1          \n" \
    "fmul.d   ft2, fs2, ft1          \n"

#define FP0_FP1_ASM_BODY                 \
    "fmul.d   fa3, %[InvLn2N], ft0   \n" \
    "fmul.d   ft3, %[InvLn2N], ft0   \n" \
    "fmul.d   ft4, %[InvLn2N], ft0   \n" \
    "fmul.d   ft5, %[InvLn2N], ft0   \n" \
    "fadd.d   ft2, fa3, %[SHIFT]     \n" \
    "fadd.d   ft2, ft3, %[SHIFT]     \n" \
    "fadd.d   ft2, ft4, %[SHIFT]     \n" \
    "fadd.d   ft2, ft5, %[SHIFT]     \n" \
    "fsub.d   fa2, ft1, %[SHIFT]     \n" \
    "fsub.d   ft6, ft1, %[SHIFT]     \n" \
    "fsub.d   ft7, ft1, %[SHIFT]     \n" \
    "fsub.d   ft8, ft1, %[SHIFT]     \n" \
    "fmv.d    ft2, fa3               \n" \
    "fmv.d    ft2, ft3               \n" \
    "fmv.d    ft2, ft4               \n" \
    "fmv.d    ft2, ft5               \n" \
    "fsub.d   fa3, ft1, fa2          \n" \
    "fsub.d   ft3, ft1, ft6          \n" \
    "fsub.d   ft4, ft1, ft7          \n" \
    "fsub.d   ft5, ft1, ft8          \n" \
    "fmadd.d  fa2, %[C0], fa3, %[C1] \n" \
    "fmadd.d  ft6, %[C0], ft3, %[C1] \n" \
    "fmadd.d  ft7, %[C0], ft4, %[C1] \n" \
    "fmadd.d  ft8, %[C0], ft5, %[C1] \n" \
    "fmadd.d  fa4, %[C2], fa3, %[C3] \n" \
    "fmadd.d  fs0, %[C2], ft3, %[C3] \n" \
    "fmadd.d  fs1, %[C2], ft4, %[C3] \n" \
    "fmadd.d  fs2, %[C2], ft5, %[C3] \n" \
    "fmul.d   fa1, fa3, fa3          \n" \
    "fmul.d   fa5, ft3, ft3          \n" \
    "fmul.d   fa6, ft4, ft4          \n" \
    "fmul.d   fa7, ft5, ft5          \n" \
    "fmadd.d  fa4, fa2, fa1, fa4     \n" \
    "fmadd.d  fs0, ft6, fa5, fs0     \n" \
    "fmadd.d  fs1, ft7, fa6, fs1     \n" \
    "fmadd.d  fs2, ft8, fa7, fs2     \n" \
    "fmul.d   ft2, fa4, ft0          \n" \
    "fmul.d   ft2, fs0, ft0          \n" \
    "fmul.d   ft2, fs1, ft0          \n" \
    "fmul.d   ft2, fs2, ft0          \n"