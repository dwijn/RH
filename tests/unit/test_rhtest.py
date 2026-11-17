#!/usr/bin/env python3
"""Self-test of the readers and comparisons in tests/rhtest on synthetic files
(standard library only, no RH program needed): pops files with an NLTE and an
LTE block, spectrum.out of both layouts with the vacuum-to-air record that
tells them apart (Q-21), solveray output, rhtest.compare, the command line of
compare_spectra.py, and the small helpers run_tests.py relies on. The formats
are the ones of writespect_xdr.c etc. (see rhtest/output.py), written here
independently with struct."""
import os
import struct
import subprocess
import sys
import tempfile
import unittest

TESTS = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, TESTS)
from rhtest import compare, keywords, output, run   # noqa: E402


def be(values):
    return struct.pack(f">{len(values)}d", *values)


class Files(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def write(self, name, data):
        p = os.path.join(self.tmp.name, name)
        with open(p, "wb") as f:
            f.write(data)
        return p

    # spectrum.out: Nspect, lambda, I[Nspect][n], bool + 200.0, then Q, U, V
    def spectrum(self, name, ns=3, n=4, polarized=False, scale=1.0, bump=None):
        lam = [400.0 + k for k in range(ns)]
        I = [scale * (1.0 + 0.1 * k + 0.01 * j) for k in range(ns) for j in range(n)]
        data = struct.pack(">i", ns) + be(lam) + be(I) + struct.pack(">id", 1, 200.0)
        if polarized:
            for s, base in enumerate((0.01, -0.02, 0.03)):
                arr = [base * (1.0 + k + 0.5 * j) for k in range(ns) for j in range(n)]
                if bump and bump[0] == "QUV"[s]:
                    arr[bump[1]] += bump[2]
                data += be(arr)
        return self.write(name, data), lam, I


class Readers(Files):
    def test_pops_nlte_and_lte(self):
        nl, ns = 3, 5
        nlte = [1.0 + k for k in range(nl * ns)]
        lte = [100.0 + k for k in range(nl * ns)]
        p = self.write("pops.X.out", b"\0" * 12 + struct.pack(">2i", nl, ns) + be(nlte) + be(lte))
        self.assertEqual(output.read_pops_array(p)[:2], (nl, ns))
        self.assertEqual(list(output.read_pops_array(p)[2]), nlte)
        self.assertEqual(list(output.read_pops_array(p, lte=True)[2]), lte)
        self.assertEqual(list(output.read_pops(p)), nlte)
        bad = self.write("pops.bad.out", b"\0" * 12 + struct.pack(">2i", nl, ns) + be(nlte))
        with self.assertRaises(ValueError):
            output.read_pops_array(bad)

    def test_spectrum_layouts(self):
        p, lam, I = self.spectrum("u.out")
        raw, l, arrs, ns = output.read_spectrum_file(p)
        self.assertEqual((ns, list(l), list(arrs["I"]), sorted(arrs)), (3, lam, I, ["I"]))
        q, lam, I = self.spectrum("p.out", polarized=True)
        raw, l, arrs, ns = output.read_spectrum_file(q)
        self.assertEqual(sorted(arrs), ["I", "Q", "U", "V"])
        self.assertEqual(list(arrs["I"]), I)
        self.assertAlmostEqual(arrs["V"][5], 0.03 * (1.0 + 1 + 0.5 * 1))
        # Q-21: the size of the polarized file also fits an unpolarized layout
        # with four times the data; only the record after I tells them apart
        with self.assertRaises(ValueError):
            output.read_spectrum_file(q, stokes=False)
        with self.assertRaises(ValueError):
            output.read_spectrum_file(p, stokes=True)
        self.assertEqual(sorted(output.read_spectrum_file(q, stokes=True)[2]), ["I", "Q", "U", "V"])
        short = self.write("s.out", b"\0\0")
        with self.assertRaises(Exception):
            output.read_spectrum_file(short)
        junk = self.write("j.out", struct.pack(">i", 3) + b"\0" * 100)
        with self.assertRaises(ValueError):
            output.read_spectrum_file(junk)

    def test_old_spectrum_readers(self):
        p, lam, I = self.spectrum("u.out", ns=3, n=4)
        l, rows = output.read_spectrum(p, 4)
        self.assertEqual((list(l), [list(r) for r in rows]), (lam, [I[0:4], I[4:8], I[8:12]]))
        q, lam, I = self.spectrum("p.out", ns=3, n=4, polarized=True)
        l, st = output.read_stokes_spectrum(q, 4)
        self.assertEqual([list(r) for r in st[0]], [I[0:4], I[4:8], I[8:12]])
        self.assertEqual(len(st), 4)

    def test_solveray(self):
        M = 3
        I, Q, U, V = ([s * (k + 1.0) for k in range(M)] for s in (1, 10, 100, 1000))
        data = be([0.5, 0.8]) + be(I) + struct.pack(">i", 0) + be(Q) + be(U) + be(V)
        p = self.write("spectrum_1", data)
        self.assertEqual([list(a) for a in output.read_solveray(p)], [I, Q, U, V])
        u = self.write("spectrum_2", be([0.5, 0.8]) + be(I) + struct.pack(">i", 0))
        self.assertEqual(list(output.read_solveray_intensity(u, M)), I)
        with self.assertRaises(ValueError):
            output.read_solveray(u)

    def test_doubles(self):
        p = self.write("J.dat", struct.pack("=3d", 1.0, 2.0, 3.0))
        self.assertEqual(list(output.read_doubles(p)), [1.0, 2.0, 3.0])
        with self.assertRaises(struct.error):
            output.doubles(b"\0" * 12, 0, 2)


class Compare(Files):
    def cmp(self, a, b, **kw):
        return compare.compare_spectra(a, b, **kw)

    def test_unpolarized(self):
        a, _, _ = self.spectrum("a.out")
        b, _, _ = self.spectrum("b.out", scale=1.0 + 1e-12)
        c, _, _ = self.spectrum("c.out", scale=1.0 + 1e-6)
        self.assertEqual(self.cmp(a, a, tol=1e-9), (0, "match (tol 1e-09): I=0.0e+00"))
        self.assertEqual(self.cmp(a, a, exact=True), (0, "byte-identical"))
        self.assertEqual(self.cmp(a, b, tol=1e-9)[0], 0)
        s, msg = self.cmp(a, c, tol=1e-9)
        self.assertEqual(s, 1)
        self.assertTrue(msg.startswith("MISMATCH (tol 1e-09): I=1.0e-06  worst I at lambda=40"), msg)
        s, msg = self.cmp(a, b, exact=True)
        self.assertEqual(s, 1)
        self.assertTrue(msg.startswith("NOT BYTE-IDENTICAL: ") and msg.endswith(" differing bytes"), msg)

    def test_polarized(self):
        a, _, _ = self.spectrum("a.out", polarized=True)
        b, _, _ = self.spectrum("b.out", polarized=True, bump=("V", 3, 3e-5))
        s, msg = self.cmp(a, b, tol=1e-9)
        self.assertEqual(s, 1)
        self.assertIn("I=0.0e+00  Q=0.0e+00  U=0.0e+00  V=", msg)
        self.assertEqual(self.cmp(a, a, tol=1e-9, stokes=True)[0], 0)
        s, msg = self.cmp(a, a, tol=1e-9, stokes=False)
        self.assertEqual(s, 2)
        self.assertTrue(msg.startswith("UNREADABLE: ") and "unpolarized layout" in msg, msg)

    def test_problems(self):
        a, _, _ = self.spectrum("a.out")
        b, _, _ = self.spectrum("b.out", ns=4)
        self.assertEqual(self.cmp(a, b)[0], 1)
        self.assertTrue(self.cmp(a, b)[1].startswith("WAVELENGTH GRID DIFFERS"))
        c, _, _ = self.spectrum("c.out", n=5)
        self.assertEqual(self.cmp(a, c)[1].split(":")[0], "SIZE DIFFERS")
        s, msg = self.cmp(a, os.path.join(self.tmp.name, "missing.out"))
        self.assertEqual(s, 2)
        self.assertTrue(msg.startswith("UNREADABLE: [Errno 2]"), msg)
        nan = self.write("n.out", struct.pack(">i", 3) + be([400.0, 401.0, 402.0]) +
                         be([float("nan")] + [1.0] * 11) + struct.pack(">id", 1, 200.0))
        self.assertEqual(self.cmp(a, nan), (1, "NON-FINITE VALUES in I"))

    def test_command_line(self):
        a, _, _ = self.spectrum("a.out", polarized=True)
        b, _, _ = self.spectrum("b.out", polarized=True, bump=("Q", 0, 1e-4))
        for args in (["--tol", "1e-9"], ["--tol", "1e-9", "--stokes"], ["--exact"],
                     ["--no-stokes"], ["--tol", "1"]):
            r = subprocess.run([sys.executable, "-I", os.path.join(TESTS, "compare_spectra.py"),
                                a, b, *args], stdout=subprocess.PIPE, universal_newlines=True)
            tol = float(args[args.index("--tol") + 1]) if "--tol" in args else 1e-8
            st = True if "--stokes" in args else (False if "--no-stokes" in args else None)
            status, msg = self.cmp(a, b, tol=tol, exact="--exact" in args, stokes=st)
            self.assertEqual((r.returncode, r.stdout), (status, msg + "\n"), args)


class Helpers(Files):
    def test_has_stokes_input(self):
        for text, expect in (("  STOKES_INPUT = field.B\n", True), ("STOKES_INPUT = none\n", False),
                             ("STOKES_INPUT=none # no field\n", False),
                             ("  STOKES_INPUT = nonexistent.B\n", True),
                             ("# STOKES_INPUT = field.B\n", False), ("N_THREADS = 1\n", False)):
            p = self.write("keyword.input", text.encode())
            self.assertEqual(keywords.has_stokes_input(p), expect, text)

    def test_log_problem(self):
        for text, expect in (("all fine\n", ""), ("x\nERROR TERMINATING_ERROR: boom\nnext\nmore\n",
                                                  "ERROR TERMINATING_ERROR: boom next "),
                             ("x NAN y\n", "NaN in log"), ("nan values\n", ""),
                             ("Segmentation fault\n", "Segmentation fault "),
                             ("seg\nSegmentation\nrt: runtime error: y\n", "rt: runtime error: y ")):
            p = self.write("log.log", text.encode())
            self.assertEqual(run.log_problem(p), expect, text)
        p = self.write("log.log", ("A" * 200 + " runtime error\n").encode())
        self.assertEqual(len(run.log_problem(p)), 150)

    def test_peak_rss(self):
        d = self.tmp.name
        self.assertIsNone(run.peak_rss_mb(d))
        if sys.platform != "darwin":
            with open(os.path.join(d, "time.txt"), "w") as f:
                f.write("\tMaximum resident set size (kbytes): 102400\n")
            self.assertEqual(run.peak_rss_mb(d), 100)
            with open(os.path.join(d, "time.txt"), "w") as f:
                f.write("\tExit status: 0\n")
            self.assertIsNone(run.peak_rss_mb(d))


if __name__ == "__main__":
    unittest.main()
