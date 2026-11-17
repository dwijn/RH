"""keyword.input and atoms.input."""
import os
import re


def setkw(text, key, value):
    """Set key = value in a keyword.input text (replacing the first
    uncommented setting, else appending it); value None removes the key by
    commenting it out, so that it is not set at all."""
    out, done = [], False
    for l in text.splitlines():
        s = l.strip()
        if not s.startswith("#") and "=" in s and s.split("=")[0].strip() == key:
            out.append(f"# {s}" if value is None else f"  {key} = {value}")
            done = True
        else:
            out.append(l)
    if not done and value is not None:
        out.append(f"  {key} = {value}")
    return "\n".join(out) + "\n"


def apply(text, pairs):
    """setkw for each (key, value) of pairs, in order."""
    for k, v in pairs:
        text = setkw(text, k, v)
    return text


def edit_file(path, pairs):
    """apply(), in place."""
    with open(path) as f:
        text = f.read()
    with open(path, "w") as f:
        f.write(apply(text, pairs))


def atoms_input(atoms):
    """Text of an atoms.input: atoms is a list of tuples (file, status,
    populations[, populations output file]), e.g.
    ("../../Atoms/CaII.atom", "ACTIVE", "ZERO_RADIATION", "pops.CA.out")."""
    return f"{len(atoms)}\n" + "".join("  " + "  ".join(a) + "\n" for a in atoms)


def write_atoms_input(d, atoms):
    with open(os.path.join(d, "atoms.input"), "w") as f:
        f.write(atoms_input(atoms))


# H passive in LTE and Ca II active from zero radiation (PRD lines), the
# atom set of the solver and PRD tests; run directories sit two levels down
H6_CAII = [("../../Atoms/H_6.atom", "PASSIVE", "LTE_POPULATIONS"),
           ("../../Atoms/CaII.atom", "ACTIVE", "ZERO_RADIATION", "pops.CA.out")]


def has_stokes_input(path):
    """The keyword.input at path sets STOKES_INPUT to something other than
    none (the run reads a magnetic field, so spectrum.out has Q, U, V).
    "none" must be a word (followed by white space, # or the end of the line);
    the grep of the old run_tests.sh took "nonexistent.B" for it too."""
    with open(path, errors="replace") as f:
        text = f.read()
    return (re.search(r"^[ \t\r\f\v]*STOKES_INPUT[ \t\r\f\v]*=", text, re.M) is not None and
            re.search(r"^[ \t\r\f\v]*STOKES_INPUT[ \t\r\f\v]*=[ \t\r\f\v]*none([ \t\r\f\v#]|$)",
                      text, re.M) is None)
