ids:
date: 2026-11-27
title: Errors that return, the mechanism: rh_status, Fail, RH_TRY, ExitOnError, RH_ASSERT (no error site converted yet)

The first part of Q-51 (errors that return, phase 3 of the library track). A function that can fail will return an `rh_status`, with the codes of the API in `docs/LIBRARY.md` §4, and a value it computes through an out-parameter. This change adds the pieces, in `error.h` and `error.c`; no error site uses them yet, so nothing else changes:

- `Fail(rh, code, routine, format, ...)` records the error in `rh->error` (code, routine and the whole formatted message) and returns the code. The first error is kept: worker threads that fail at the same moment race on an atomic flag, and the winner writes the record.
- `RH_TRY(call)` passes a failure on to the caller.
- `ExitOnError(rh, status)` is the bridge for code not yet converted: on a failure it prints the recorded error exactly as `Error()` prints a terminating error and exits with status 70, so a run behaves as before while the conversion moves through the library module by module.
- `RH_ASSERT(rh, condition, routine, format, ...)` is for invariants, which only a bug can break: it prints the message, the condition and where it is, and aborts, in Release builds too.
- `Failed(rh)` tells workers to stop early; `ClearError(rh)` forgets the record.

`Error()` and `ExitOnError()` share the code that prints the message and exits. The new unit test `test_status` checks each piece, among them that `ExitOnError` prints byte for byte what `Error` prints, and that eight threads failing at once leave one whole record.

Phase 3 of the library track (`docs/LIBRARY.md` §6.3); Q-51 is closed by the next part, which converts the runtime and numerics.
