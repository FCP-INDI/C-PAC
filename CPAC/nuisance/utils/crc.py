import base64

POLYNOMIAL = 0x1021
PRESET = 0xFFFF


def _initial(c):
    crc = 0
    c = c << 8
    for _ in range(8):
        if (crc ^ c) & 0x8000:
            crc = (crc << 1) ^ POLYNOMIAL
        else:
            crc = crc << 1
        c = c << 1
    return crc

_tab = [_initial(i) for i in range(256)]


def _update_crc(crc, c):
    cc = 0xff & c

    tmp = (crc >> 8) ^ cc
    crc = (crc << 8) ^ _tab[tmp & 0xff]
    crc = crc & 0xffff

    return crc


def crc(str):
    crc = PRESET
    for c in str:
        crc = _update_crc(crc, ord(c))
    return crc


def encode(string):
    scrc = str(crc(string))
    bcrc = scrc.encode('ascii')
    base64crc = base64.urlsafe_b64encode(bcrc).strip(b'=')
    return base64crc.decode()
