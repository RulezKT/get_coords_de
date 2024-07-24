import array
import math

from pydantic import BaseModel

# constants for file de440s.bsp
SEGMENT_START_TIME: int = -4_734_072_000
SEGMENT_LAST_TIME: int = 4_735_368_000


class FileRecords(BaseModel):
    rec_start_addr: int
    int_len: float
    rsize: float


DE440S_FILE_RECORDS: tuple[FileRecords] = (
    {
        # SSB
        "rec_start_addr": 0,
        "int_len": 0,
        "rsize": 0,
    },
    {
        # MERCURY_BARYCENTER
        "rec_start_addr": 8065,
        "int_len": 691200.0,
        "rsize": 44.0,
    },
    {
        # VENUS_BARYCENTER
        "rec_start_addr": 610869,
        "int_len": 1382400.0,
        "rsize": 32.0,
    },
    {
        # EARTH_BARYCENTER
        "rec_start_addr": 830073,
        "int_len": 1382400.0,
        "rsize": 41.0,
    },
    {
        # MARS_BARYCENTER
        "rec_start_addr": 1110927,
        "int_len": 2764800.0,
        "rsize": 35.0,
    },
    {
        # JUPITER_BARYCENTER
        "rec_start_addr": 1230806,
        "int_len": 2764800.0,
        "rsize": 26.0,
    },
    {
        # SATURN_BARYCENTER
        "rec_start_addr": 1319860,
        "int_len": 2764800.0,
        "rsize": 23.0,
    },
    {
        # URANUS_BARYCENTER
        "rec_start_addr": 1398639,
        "int_len": 2764800.0,
        "rsize": 20.0,
    },
    {
        # NEPTUNE_BARYCENTER
        "rec_start_addr": 1467143,
        "int_len": 2764800.0,
        "rsize": 20.0,
    },
    {
        # PLUTO_BARYCENTER
        "rec_start_addr": 1535647,
        "int_len": 2764800.0,
        "rsize": 20.0,
    },
    {
        # SUN
        "rec_start_addr": 1604151,
        "int_len": 1382400.0,
        "rsize": 35.0,
    },
    {
        # MOON
        "rec_start_addr": 1843905,
        "int_len": 345600.0,
        "rsize": 41.0,
    },
    {
        # EARTH
        "rec_start_addr": 2967309,
        "int_len": 345600.0,
        "rsize": 41.0,
    },
)


class Position(BaseModel):
    x: float
    y: float
    z: float


# v2024
# get coords x , y , z from de440s.bsp file
# target is 0 for all planets except for the moon
# moon is calculated relative to earth_bary so target is 3 , so 301/3 or 11/3
# earth(399 or 12) is calculated relative to earth_bary so target is 3, 399/3 or 12/3
# earth_bary(3) is calculated as usual so target is 0
def get_coords(date_in_seconds: int, target_code: int, file_in_mem: bytes) -> Position:
    if SEGMENT_START_TIME >= date_in_seconds or SEGMENT_LAST_TIME <= date_in_seconds:
        print("get_coords: Date is out of range")
        return {"x": 0, "y": 0, "z": 0}

    i = target_code
    if target_code == 399:
        i = 12
    elif target_code == 301:
        i = 11
    if i > 12:
        print("get_coords: Target code is out of range")
        return {"x": 0, "y": 0, "z": 0}

    intlen = DE440S_FILE_RECORDS[i]["int_len"]
    rsize = DE440S_FILE_RECORDS[i]["rsize"]
    start_addr = DE440S_FILE_RECORDS[i]["rec_start_addr"]

    internal_offset = (
        math.floor((date_in_seconds - SEGMENT_START_TIME) / intlen) * rsize
    )
    record = 8 * int(start_addr + internal_offset)

    order = (int(rsize) - 2) / 3 - 1

    # the only place, we read data from bsp file, loaded in memory
    data = array.array("d", file_in_mem[record - 8 : record + int(rsize) * 8])

    tau = (date_in_seconds - data[0]) / data[1]
    order = int(order)
    deg = order + 1

    return {
        "x": _chebyshev(order, tau, data[2 : 2 + deg]),
        "y": _chebyshev(order, tau, data[2 + deg : 2 + 2 * deg]),
        "z": _chebyshev(order, tau, data[2 + 2 * deg : 2 + 3 * deg]),
    }


def _chebyshev(order: int, x: float, data: list[float]) -> float:
    # Evaluate a Chebyshev polynomial
    two_x = 2 * x
    bkp2 = data[order]
    bkp1 = two_x * bkp2 + data[order - 1]
    for n in range(order - 2, 0, -1):
        bk = data[n] + two_x * bkp1 - bkp2
        bkp2 = bkp1
        bkp1 = bk

    return data[0] + x * bkp1 - bkp2
