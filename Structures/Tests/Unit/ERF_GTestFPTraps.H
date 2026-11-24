#ifndef ERF_GTEST_FP_TRAPS_H
#define ERF_GTEST_FP_TRAPS_H

#include <cfenv>

#include <AMReX.H>

namespace erf_gtest {

/**
 * amrex::setFPExcept(excepts) with the floating-point flags cleared first; returns the earlier traps. A test
 * that ran with the traps off may leave a flag raised. On x86-64 glibc that flag can sit in the x87 status
 * word as well as in the SSE register: fesetexceptflag writes both, and erf_moordyn::FPTrapMask's destructor
 * calls it to give back the flags raised before the mask. feenableexcept then unmasks it, and the next x87
 * instruction that waits (the fstcw in fegetexcept, that is the next amrex::getFPExcept) stops the test
 * binary with SIGFPE for an operation of an earlier test. Every scope below that sets traps goes through here.
 */
inline amrex::FPExcept set_fp_traps (amrex::FPExcept excepts)
{
    std::feclearexcept(FE_ALL_EXCEPT);
    return amrex::setFPExcept(excepts);
}

/**
 * AMReX's three floating-point traps (invalid, zero and overflow, as amrex.fpe_trap_* sets them) on for a
 * scope, the earlier traps restored after it, a failed assertion included. active() is false where AMReX
 * cannot set traps (not Linux glibc or macOS): a test then skips. The flags raised before the scope are
 * cleared before the traps come on (set_fp_traps), and those raised in it before the earlier traps come back,
 * so that none fires later from a sticky flag.
 */
class FPTrapsOn
{
public:
    FPTrapsOn () : m_prev(set_fp_traps(amrex::FPExcept::all)) {}
    ~FPTrapsOn () { std::feclearexcept(FE_ALL_EXCEPT); amrex::setFPExcept(m_prev); }
    FPTrapsOn (const FPTrapsOn&) = delete;
    FPTrapsOn& operator= (const FPTrapsOn&) = delete;
    static bool active () { return amrex::getFPExcept() == amrex::FPExcept::all; }

private:
    amrex::FPExcept m_prev;
};

/**
 * The floating-point flags among invalid, divide-by-zero and overflow that f() raises, with the traps off for the
 * call (the earlier traps restored after it). A guarded division that an optimising compiler evaluates
 * speculatively raises its flag although the result is discarded; this reports it as a failed assertion instead of
 * the SIGILL that amrex.fpe_trap_* gives, which would stop the whole test binary. Only an optimised build
 * speculates, so a Debug build passes either way.
 */
template <class F>
int raised_fp_flags (F&& f)
{
    const amrex::FPExcept prev = amrex::setFPExcept(amrex::FPExcept::none);
    std::feclearexcept(FE_ALL_EXCEPT);
    f();
    const int raised = std::fetestexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
    std::feclearexcept(FE_ALL_EXCEPT);
    amrex::setFPExcept(prev);
    return raised;
}

} // namespace erf_gtest

#endif
