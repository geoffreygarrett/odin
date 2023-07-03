# Long Double Data Type

## Introduction

`long double` is a floating-point data type in C and its related languages, providing more precision than the `double`
type in most cases. However, the language standard only necessitates it to be at least as precise as `double`, and its
actual implementation may not adhere strictly to the IEEE formats.

## Historical Background

The `long double` type was introduced in the original 1989 C standard, and further support was added with the C99
standard revision in 1999, extending the standard library to include `long double` operations such as `sinl()`
and `strtold()`[^1^]. `long double` constants are floating-point constants appended with an "L" or "l" suffix. Absence
of a suffix leaves the evaluation dependent on `FLT_EVAL_METHOD`.

## Implementations

The implementation of `long double` differs across architectures and compilers.

### x86 Architecture

On the x86 architecture, most C compilers map `long double` to the 80-bit extended precision type, as outlined in the
C99/C11 standards and supported by x86 hardware[^16^]. However, there are exceptions, such as Microsoft Visual C++ for
x86, where `long double` is equivalent to `double`[^2^].

Intel's C++ compiler on Microsoft Windows supports extended precision, although it requires the `/Qlong-double` switch
for `long double` to match the hardware's extended precision format[^3^].

Several compilers utilize `long double` for the IEEE 754 quadruple-precision binary floating-point format (binary128).
This holds true for operating systems such as HP-UX[^4^], Solaris/SPARC[^5^], 64-bit ARM (AArch64)[^7^], and z/OS
with `FLOAT(IEEE)`[^8^][^9^][^10^].

### PowerPC Systems

On some PowerPC systems, `long double` is implemented as double-double arithmetic[^11^], considering a `long double`
value as the exact sum of two double-precision values and providing at least 106-bit precision.

### GNU C Compiler

With the GNU C Compiler, `long double` is 80-bit extended precision on x86 processors[^16^], and on certain
architectures, `long double` can be double-double or 128-bit quadruple precision[^20^]. As of gcc 4.3, quadruple
precision is supported on x86, but it's implemented as the nonstandard type `__float128` rather
than `long double`[^21^].

## Processor Configurations

It's worth noting that while x86 architecture supports 80-bit extended-precision operations, processor configuration can
automatically round operations to double or single precision. Extended precision may be employed for intermediate
compiler-generated calculations even if the final results are stored at a lower precision. This is dependent on
the `FLT_EVAL_METHOD` and can be overridden within individual programs[^22^][^23^][^24^][^25^][^26^].

## Other Specifications

In CORBA (from the specification of 3.0), the `long double` data type represents an IEEE double-extended floating-point
number, having an exponent of at least 15 bits in length and a signed fraction of at least 64 bits[^27^].

## References

[^1^]: ANSI/ISO 9899-1990 American National Standard for Programming Languages - C, section 6.1.2.5.
[^2^]: "Long Double". learn.microsoft.com. Retrieved 2022-10-06.
[^3^]: Intel Developer Site
[^4^]: Hewlett Packard (1992). "Porting

C Programs". HP-UX Portability Guide - HP 9000 Computers (PDF) (2nd ed.). pp. 5-3 and 5-37.
[^5^]: "IEEE Arithmetic". docs.oracle.com. Retrieved 2022-10-06.
[^7^]: "Procedure Call Standard for the Arm® 64-bit Architecture (AArch64)". GitHub. 2020-10-01. Archived (PDF) from the
original on 2020-10-02.
[^8^]: "Floating-point types". IBM. 2020-10-09. Retrieved 2020-10-09.
[^9^]: Schwarz, Eric (June 22, 2015). "The IBM z13 SIMD Accelerators for Integer, String, and Floating-Point" (PDF).
Retrieved July 13, 2015.
[^10^]: Schwarz, E. M.; Krygowski, C. A. (September 1999). "The S/390 G5 floating-point unit". IBM Journal of Research
and Development. 43 (5/6): 707–721. CiteSeerX 10.1.1.117.6711. doi:10.1147/rd.435.0707.
[^11^]: "The saga of the Power ISA 128-bit long double". 2018-12-22. Retrieved 2021-12-26.
[^16^]: "x86 Options (Using the GNU Compiler Collection (GCC))". gcc.gnu.org. Retrieved 2022-10-06.
[^20^]: "SPARC Options (Using the GNU Compiler Collection (GCC))". gcc.gnu.org. Retrieved 2022-10-06.
[^21^]: "GCC 4.3 Release Series — Changes, New Features, and Fixes - GNU Project". gcc.gnu.org. Retrieved 2022-10-06.
[^22^]: Brian J. Gough and Richard M. Stallman, An Introduction to GCC, section 8.6 Floating-point issues (Network
Theory Ltd., 2004).
[^23^]: "Significant changes from NetBSD 6.0 to 7.0".
[^24^]: "Visual Studio 2005 Retired documentation". Microsoft Download Center. Retrieved 2022-10-06.
[^25^]: Intel C++ Compiler Documentation, Using the -fp-model (/fp) Option.
[^26^]: "IA-32 Function Calling Conventions".
[^27^]: "CORBA 3.0 specification"
