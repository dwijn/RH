"""The work directory of a test: <path>/runs/<name> run directories, and
<path>/{Atoms,Atmos,Molecules} links to the source tree (the run directories
sit two levels down, so their input files say ../../Atoms/...)."""
import atexit
import os
import shutil
import sys
import tempfile


class WorkDir:
    """work given (--work): used and always kept. Otherwise a temporary
    directory under $TMPDIR, removed by finish(True); kept and named when the
    script ends any other way (finish(False), sys.exit, an uncaught
    exception)."""

    def __init__(self, src, work=None):
        self.src = os.path.abspath(src)
        self.explicit = bool(work)
        self.path = os.path.abspath(work) if work else tempfile.mkdtemp(prefix="rh_e2e_")
        self.done = False
        atexit.register(self._name_if_kept)
        os.makedirs(os.path.join(self.path, "runs"), exist_ok=True)
        for d in ("Atoms", "Atmos", "Molecules"):
            link = os.path.join(self.path, d)
            if not os.path.lexists(link):
                os.symlink(os.path.join(self.src, d), link)

    def join(self, *parts):
        return os.path.join(self.path, *parts)

    def _name_if_kept(self):
        if not self.done:
            self.done = True
            sys.stdout.flush()
            print(f"runs kept in {self.path}", file=sys.stderr)

    def finish(self, ok):
        """End of the test: remove a temporary directory when ok, otherwise
        keep and name it. Call once, last."""
        self.done = True
        if ok and not self.explicit:
            shutil.rmtree(self.path)
        elif not ok:
            sys.stdout.flush()
            print(f"runs kept in {self.path}", file=sys.stderr)
