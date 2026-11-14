"""Attribute build-log warnings to FLARE or to Submodules/ the way FLARE CI's Report step does."""

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

import build_warnings

ROOT = "/home/runner/work/FLARE/FLARE"
ERF_TU = f"{ROOT}/Submodules/ERF/Source/IO/ERF_WriteScalarProfiles.cpp"
FLARE_TU = f"{ROOT}/Structures/Source/FLARE_Example.cpp"
HEADER_WARNING = ("/usr/include/c++/11/bits/stl_algobase.h:431:30: warning: "
                  "'void* __builtin_memcpy(void*, const void*, long unsigned int)' writing 32 bytes "
                  "into a region of size 16 overflows the destination [-Wstringop-overflow=]")


def gcc_block(tu, frame):
    """The GCC 11 diagnostic from FLARE CI run 43603472994, with its TU and last frame swapped in."""

    return f"""\
[ 38%] Building CXX object Submodules/ERF/Exec/CMakeFiles/erf_srclib.dir/__/Source/IO/ERF_Plotfile2DMetadata.cpp.o
In file included from /usr/include/c++/11/bits/char_traits.h:39,
                 from /usr/include/c++/11/string:40,
                 from /usr/include/c++/11/iomanip:40,
                 from {tu}:4:
In static member function 'static constexpr _Tp* std::__copy_move<_IsMove, true, std::random_access_iterator_tag>::__copy_m(const _Tp*, const _Tp*, _Tp*) [with _Tp = double; bool _IsMove = false]',
    inlined from 'constexpr _OI std::copy(_II, _II, _OI) [with _II = const double*; _OI = double*]' at /usr/include/c++/11/bits/stl_algobase.h:620:7,
    inlined from 'std::vector<_Tp, _Alloc>::iterator std::vector<_Tp, _Alloc>::insert(std::vector<_Tp, _Alloc>::const_iterator, std::initializer_list<_Tp>) [with _Tp = double; _Alloc = std::allocator<double>]' at /usr/include/c++/11/bits/stl_vector.h:1313:17,
    inlined from 'void ERF::sum_integrated_quantities(double)' at {frame}:172:17:
{HEADER_WARNING}
  431 |             __builtin_memmove(__result, __first, sizeof(_Tp) * _Num);
      |             ~~~~~~~~~~~~~~~~~^~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/usr/include/c++/11/ext/new_allocator.h: In member function 'void ERF::sum_integrated_quantities(double)':
/usr/include/c++/11/ext/new_allocator.h:127:48: note: at offset 48 into destination object of size 64 allocated by 'operator new'
""".splitlines()


def old_grep(lines):
    """The filter the Report step used before: drop lines naming Submodules/ or lto-wrapper."""

    return sorted({l for l in lines if ("warning:" in l.lower() or "error:" in l.lower())
                   and "Submodules/" not in l and "lto-wrapper: warning:" not in l})


class AttributionTests(unittest.TestCase):
    def kept(self, lines):
        return build_warnings.flare_warnings(lines, ROOT)

    def test_header_warning_inlined_from_erf_is_dropped(self):
        lines = gcc_block(ERF_TU, ERF_TU)
        self.assertEqual(old_grep(lines), [HEADER_WARNING])  # what turned the job red
        self.assertEqual(self.kept(lines), [])

    def test_header_warning_inlined_from_flare_is_kept(self):
        self.assertEqual(self.kept(gcc_block(FLARE_TU, FLARE_TU)), [HEADER_WARNING])

    def test_innermost_frame_decides_over_translation_unit(self):
        # FLARE code inlined into an ERF translation unit is FLARE's; the reverse is ERF's.
        self.assertEqual(self.kept(gcc_block(ERF_TU, FLARE_TU)), [HEADER_WARNING])
        self.assertEqual(self.kept(gcc_block(FLARE_TU, ERF_TU)), [])

    def test_first_root_frame_is_the_innermost(self):
        # GCC lists inline frames innermost first: a FLARE helper inlined into an ERF caller
        # is FLARE's, an ERF helper inlined into a FLARE caller is ERF's.
        for helper, caller, expected in ((FLARE_TU, ERF_TU, [HEADER_WARNING]), (ERF_TU, FLARE_TU, [])):
            with self.subTest(helper=helper):
                lines = gcc_block(caller, caller)
                last = next(i for i, l in enumerate(lines) if "inlined from 'void ERF" in l)
                lines.insert(last, f"    inlined from 'void helper()' at {helper}:9:5,")
                self.assertEqual(self.kept(lines), expected)

    def test_constructor_and_operator_context_lines(self):
        for header in ("In copy constructor 'X::X(const X&)',", "In move constructor 'X::X(X&&)',",
                       "In member function 'X& X::operator=(const X&)',", "In function template 'f',"):
            with self.subTest(header=header):
                lines = gcc_block(ERF_TU, ERF_TU)
                lines[lines.index(next(l for l in lines if l.startswith("In static")))] = header
                self.assertEqual(self.kept(lines), [])

    def test_translation_unit_decides_without_inline_frames(self):
        for tu, expected in ((ERF_TU, []), (FLARE_TU, [HEADER_WARNING])):
            with self.subTest(tu=tu):
                lines = [l for l in gcc_block(tu, tu) if "inlined from" not in l and "In static" not in l]
                self.assertEqual(self.kept(lines), expected)

    def test_template_required_from_frame(self):
        for tu, expected in ((ERF_TU, []), (FLARE_TU, [HEADER_WARNING])):
            with self.subTest(tu=tu):
                lines = [f"/usr/include/c++/11/bits/stl_vector.h: In instantiation of 'void f() [with T = double]':",
                         f"{tu}:50:9:   required from here", HEADER_WARNING]
                self.assertEqual(self.kept(lines), expected)

    def test_header_warning_without_frames_in_root_is_kept(self):
        self.assertEqual(self.kept(gcc_block("/usr/include/c++/11/vector", "/usr/include/c++/11/vector")),
                         [HEADER_WARNING])

    def test_context_does_not_leak_across_unrelated_lines(self):
        # A build line between the ERF context and the header warning breaks the block.
        lines = gcc_block(ERF_TU, ERF_TU)
        lines.insert(lines.index(HEADER_WARNING), "[ 39%] Building CXX object Structures/x.cpp.o")
        self.assertEqual(self.kept(lines), [HEADER_WARNING])
        # A later header warning with no context of its own is not attributed to ERF.
        self.assertEqual(self.kept(gcc_block(ERF_TU, ERF_TU) + [HEADER_WARNING]), [HEADER_WARNING])

    def test_flare_and_submodule_locations(self):
        flare = f"{FLARE_TU}:12:5: warning: unused variable 'x' [-Wunused-variable]"
        erf = f"{ERF_TU}:12:5: warning: unused variable 'x' [-Wunused-variable]"
        relative = "Structures/Source/FLARE_Example.cpp:3:1: error: expected ';'"
        lto = "lto-wrapper: warning: using serial compilation of 4 LTRANS jobs"
        self.assertEqual(self.kept([flare, erf, relative, lto, flare]), sorted([flare, relative]))

    def test_unparsed_lines_follow_the_old_filter(self):
        lines = ["CMake Warning: something in FLARE", "CMake Warning: Submodules/ERF thing"]
        self.assertEqual(self.kept(lines), old_grep(lines))


class MainTests(unittest.TestCase):
    def run_main(self, argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()):
            status = build_warnings.main(argv)
        return status, out.getvalue().splitlines()

    def test_count_line_and_trailing_slash_root(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "build-output.txt"
            log.write_text("\n".join(gcc_block(ERF_TU, ERF_TU) + gcc_block(FLARE_TU, FLARE_TU)) + "\n")
            status, out = self.run_main(["--root", ROOT + "/", str(log)])
        self.assertEqual(status, 0)
        self.assertEqual(out, [HEADER_WARNING, "Warnings: 1"])

    def test_missing_log_is_an_error(self):
        status, out = self.run_main(["--root", ROOT, "/nonexistent/build-output.txt"])
        self.assertEqual((status, out), (2, []))


if __name__ == "__main__":
    unittest.main()
