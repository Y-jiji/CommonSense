import struct
import pathlib
import numpy as np
import os

""" Utility functions for reading and writing binary files. """
# --------------------------------------------------------------
# odat file format:
# first byte is data type:
# the first 4 bits is log(length), the last four bits is type x1 - int, x2 - float, x3 unsigned
# for example 32bit float is 0x52.
# after this, the bulk data is [vec length] (int32) + data * vec length (data type defined above).


class odat:
    int32 = 0x51
    float32 = 0x52
    int64 = 0x61
    float64 = 0x62
    int16 = 0x41
    int8 = 0x31
    uint8 = 0x33
    formats = {
        int32: "i",
        float32: "f",
        int64: "q",
        float64: "d",
        int16: "h",
        int8: "b",
        uint8: "B",
    }
    extensions = {"fvecs": float32, "bvecs": int8, "ivecs": int32, "dvecs": float64}
    nptypes = {
        int32: np.int32,
        float32: np.float32,
        int64: np.int64,
        float64: np.float64,
        int16: np.int16,
        int8: np.int8,
        uint8: np.uint8,
    }
    revtypes = {
        np.dtype("int32"): int32,
        np.dtype("float32"): float32,
        np.dtype("int64"): int64,
        np.dtype("float64"): float64,
        np.dtype("int16"): int16,
        np.dtype("int8"): int8,
        np.dtype("uint8"): uint8,
    }

    def dsize(dtype):
        return 1 << (dtype // 16 - 3)


# Input type: 2D numpy array (uniform). or list of 1D numpy arrays (non-uniform)
# In these cases, no need to specify data type
# Or any 2D iterable data, but need to specify odat.dtype
def write_file(filename, data, dtype=None, create_parent=True, offset=0):
    if dtype is None:
        # if data is ndarray
        if isinstance(data, np.ndarray):
            dtype = odat.revtypes[data.dtype]
        # if data is a list of ndarrays, only the first item will be checked
        # if the input list is inconsistent, this code will be buggy
        elif hasattr(data, "__getitem__") and isinstance(data[0], np.ndarray):
            dtype = odat.revtypes[data[0].dtype]
        else:
            print("Error: Cannot deduce data type!")
            return

    if isinstance(filename, pathlib.Path):
        filename = str(filename)
    ext = filename.split(".")[-1]
    if ext != "odat":
        print("Error: Please change file extension to .odat!")
        return
    if create_parent:
        par_path = pathlib.Path(filename).parent
        par_path.mkdir(parents=True, exist_ok=True)
    open_code = "wb" if offset == 0 else "ab"

    with open(filename, open_code) as fp:
        if offset == 0:
            dt = struct.pack("b", dtype)
            fp.write(dt)
        else:
            fp.seek(offset)
        if isinstance(data, np.ndarray) and len(data.shape) == 1:
            data = data.reshape((1, -1))

        for y in data:
            d = struct.pack("I", y.size)
            fp.write(d)
            for x in y:
                a = struct.pack(odat.formats[dtype], x)
                fp.write(a)
        return fp.tell()


def see_int(array, pos):
    return array[pos : pos + 4].view(np.int32)[0]


def read_non_uniform(file, count, dtype, buffer, filesz, getoffset):
    block_size = 1048576
    pos = 0
    step = 0
    result = []
    while count > 0:
        if pos < buffer.size:
            linesz = see_int(buffer, pos)
            linesz_bytes = 4 + linesz * odat.dsize(dtype)
            if pos + linesz_bytes > buffer.size:
                pass
            else:
                vec = buffer[pos + 4 : pos + linesz_bytes].view(odat.nptypes[dtype])
                result.append(vec)
                count -= 1
                pos += linesz_bytes
                continue
        if file.tell() < filesz:
            buffer = np.fromfile(
                file, dtype=np.int8, count=block_size, offset=pos - buffer.size
            )
            pos = 0
        else:
            break
    offset = file.tell() - buffer.size + pos
    if getoffset:
        return result, offset
    else:
        return result


def read_impl(file, count, offset, dtype, filesz, getoffset):
    # returns a 2D numpy array if rows are uniform (same length)
    # otherwise returns a list of numpy arrays
    file.seek(offset)
    byte = file.read(4)
    sz = struct.unpack("I", byte)[0]
    file.seek(-4, 1)
    filesz = filesz - offset
    linesz_bytes = 4 + sz * odat.dsize(dtype)
    count = 0xFFFFFFFF if count < 0 else count
    count_non_uniform = count
    # infer count
    if count * linesz_bytes > filesz:
        if filesz % linesz_bytes == 0:
            count = filesz // linesz_bytes
        else:  # nonuniform
            return read_non_uniform(
                file, count_non_uniform, dtype, np.array([]), filesz, getoffset
            )

    buffer = np.fromfile(file, dtype=np.int8, count=count * linesz_bytes)
    buffer.resize((count, linesz_bytes))
    linesz = np.ascontiguousarray(buffer[:, :4]).view(np.int32)
    if not np.all(linesz == sz):
        buffer.resize(buffer.size)
        return read_non_uniform(
            file, count_non_uniform, dtype, buffer, filesz, getoffset
        )
    offset += count * linesz_bytes
    result = np.ascontiguousarray(buffer[:, 4:])
    result = result.view(odat.nptypes[dtype])
    if getoffset:
        return result, offset
    else:
        return result


# offset is number of bits since file start
# count is number of lines to read
# automatically infer file type from file extension
def read_file(filename, count=-1, offset=0, getoffset=False):
    file_ext = filename.split(".")[-1]
    with open(filename, "rb") as f:
        if file_ext == "odat":
            byte = f.read(1)
            dtype = struct.unpack("b", byte)[0]
            offset = 1 if offset == 0 else offset
        else:
            dtype = odat.extensions[file_ext]
        return read_impl(f, count, offset, dtype, os.path.getsize(filename), getoffset)


# read legacy file
def read_legacy_file(filename, dtype, count=-1, offset=0, getoffset=False):
    with open(filename, "rb") as f:
        return read_impl(f, count, offset, dtype, os.path.getsize(filename), getoffset)


# both reader and write do not keep file open
class ONIAKReader:
    def __init__(self, filename):
        self.filename = filename
        if not os.path.exists(filename):
            raise ValueError(filename, "file does not exist")
        self.pointer = 0

    def readline(self, line=1):
        result, offset = read_file(self.filename, line, self.pointer, getoffset=True)
        self.pointer = offset
        return result

    def reset(self):
        self.pointer = 0

    # only works if the file is uniform
    def get_size(self):
        file_ext = self.filename.split(".")[-1]
        filesize = os.path.getsize(self.filename)

        with open(self.filename, "rb") as f:
            if file_ext == "odat":
                byte = f.read(1)
                dtype = struct.unpack("b", byte)[0]
                filesize -= 1  # odat type byte
            else:
                dtype = odat.extensions[file_ext]

            byte = f.read(4)
            sz = struct.unpack("I", byte)[0]
        linesz_bytes = 4 + sz * odat.dsize(dtype)

        if filesize % linesz_bytes != 0:
            # this file is non-uniform
            return -1

        return filesize // linesz_bytes, sz


# file parents will be created if not existent yet
class ONIAKWriter:
    def __init__(self, filename, dtype = odat.float32, overwrite=False):
        self.filename = filename
        if not overwrite and os.path.exists(filename):
            raise ValueError("file cannot be overwritten")
        self.pointer = 0
        self.dtype = dtype

    def writelines(self, data):
        self.pointer = write_file(self.filename, data, dtype=self.dtype, offset=self.pointer)
