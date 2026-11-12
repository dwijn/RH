"""
Minimal XDR (RFC 4506) packing and unpacking for rhanalyze.

Replaces the standard-library module xdrlib, deprecated in Python 3.11 and
removed in 3.13 (B-48), with the subset rhanalyze uses and the same
behaviour: big-endian 4-byte integers and floats, 8-byte doubles, opaque
data and strings padded to a multiple of 4 bytes, EOFError when the buffer
runs out, Error from done() when data is left over.
"""
import struct
from io import BytesIO


class Error(Exception):
    """Exception class for this module (as xdrlib.Error)."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        return repr(self.msg)

    def __str__(self):
        return str(self.msg)


class ConversionError(Error):
    """A value cannot be packed (as xdrlib.ConversionError)."""


def _packer(fmt):
    s = struct.Struct(fmt)

    def pack(self, x):
        try:
            self._buf.write(s.pack(x))
        except struct.error as e:
            raise ConversionError(e.args[0]) from None
    return pack


def _unpacker(fmt):
    s = struct.Struct(fmt)

    def unpack(self):
        i = self._pos
        self._pos = j = i + s.size
        data = self._buf[i:j]
        if len(data) < s.size:
            raise EOFError
        return s.unpack(data)[0]
    return unpack


class Packer:
    """Pack values into an XDR buffer."""

    def __init__(self):
        self.reset()

    def reset(self):
        self._buf = BytesIO()

    def get_buffer(self):
        return self._buf.getvalue()

    pack_int = _packer('>l')
    pack_uint = _packer('>L')
    pack_float = _packer('>f')
    pack_double = _packer('>d')

    def pack_fstring(self, n, s):
        if n < 0:
            raise ValueError('fstring size must be nonnegative')
        data = s[:n]
        n = ((n + 3) // 4) * 4
        self._buf.write(data + (n - len(data)) * b'\0')

    def pack_string(self, s):
        n = len(s)
        self.pack_uint(n)
        self.pack_fstring(n, s)

    def pack_farray(self, n, list, pack_item):
        if len(list) != n:
            raise ValueError('wrong array size')
        for item in list:
            pack_item(item)


class Unpacker:
    """Unpack values from an XDR buffer."""

    def __init__(self, data):
        self.reset(data)

    def reset(self, data):
        self._buf = data
        self._pos = 0

    def get_position(self):
        return self._pos

    def set_position(self, position):
        self._pos = position

    def get_buffer(self):
        return self._buf

    def done(self):
        if self._pos < len(self._buf):
            raise Error('unextracted data remains')

    unpack_int = _unpacker('>l')
    unpack_uint = _unpacker('>L')
    unpack_float = _unpacker('>f')
    unpack_double = _unpacker('>d')

    def unpack_fstring(self, n):
        if n < 0:
            raise ValueError('fstring size must be nonnegative')
        i = self._pos
        j = i + ((n + 3) // 4) * 4
        if j > len(self._buf):
            raise EOFError
        self._pos = j
        return self._buf[i:i + n]

    def unpack_string(self):
        return self.unpack_fstring(self.unpack_uint())

    def unpack_farray(self, n, unpack_item):
        return [unpack_item() for _ in range(n)]
