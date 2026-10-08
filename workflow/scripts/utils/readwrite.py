"""Technology-agnostic access to multiplexed proteomics images.

Every technology-specific detail of reading an image (qptiff vs OME-TIFF channel names and pixel
size) lives here, so that downstream scripts only see channel names, one channel plane at a time and
a pixel size in µm. Images are never loaded whole: CODEX samples reach 70 x 24k x 24k px (~78 GB).
"""

import re
import xml.etree.ElementTree as ET

import numpy as np
import tifffile

# channel names of the two-channel segmentation input written by prepare_image and read by Mesmer
SEGMENTATION_NUCLEAR = "nuclear"
SEGMENTATION_MEMBRANE = "membrane"


class MultiplexImage:
    """Level-0 channels of a multiplexed image (qptiff or OME-TIFF), read one plane at a time.

    Use as a context manager:
        with MultiplexImage(path) as img:
            dapi = img.read_channel("DAPI")
    """

    def __init__(self, path):
        self.path = path
        self._tif = tifffile.TiffFile(path)

        # level 0 of the first series: the full-resolution image (qptiff "Baseline", OME image 0)
        series = self._tif.series[0]
        self._level = series.levels[0]
        axes = self._level.axes
        if axes not in ("CYX", "YX"):
            raise ValueError(f"{path}: unsupported axes {axes!r}, expected CYX or YX")
        # single-channel YX images are exposed as CYX with one channel
        self.shape = self._level.shape if axes == "CYX" else (1,) + tuple(self._level.shape)
        self.dtype = self._level.dtype

        # channel names and pixel size, from the format's own metadata
        if self._tif.is_ome:
            self.channel_names = _ome_channel_names(self._tif.ome_metadata)
            self.pixel_size_um = _ome_pixel_size_um(self._tif.ome_metadata)
            self.channel_offsets = _ome_channel_offsets(self._tif.ome_metadata)
        else:
            self.channel_names = [_qptiff_channel_name(p.description) for p in self._level.pages[: self.shape[0]]]
            self.pixel_size_um = None
            self.channel_offsets = [0] * len(self.channel_names)
        # fall back to the generic TIFF resolution tags (qptiff, or OME without PhysicalSizeX)
        if self.pixel_size_um is None:
            self.pixel_size_um = _tiff_resolution_um(self._level.pages[0])

        # every channel must have a unique name, since channels are selected by name downstream
        if len(self.channel_names) != self.shape[0] or any(n is None for n in self.channel_names):
            raise ValueError(f"{path}: could not read a name for every channel: {self.channel_names}")
        if len(set(self.channel_names)) != len(self.channel_names):
            raise ValueError(f"{path}: duplicated channel names: {self.channel_names}")

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def close(self):
        self._tif.close()

    @property
    def n_channels(self):
        return self.shape[0]

    def channel_index(self, name):
        if name not in self.channel_names:
            raise KeyError(f"{self.path}: no channel {name!r}; channels are {self.channel_names}")
        return self.channel_names.index(name)

    def read_channel(self, channel, maxworkers=None, subtract_offset=False):
        """One full-resolution channel plane (Y, X), by index or name.

        subtract_offset: remove the channel's pixel value offset (channel_offsets), clipping at 0.
        """
        c = self.channel_index(channel) if isinstance(channel, str) else channel
        # each channel is one TIFF page, so only that plane is read from disk
        plane = self._level.pages[c].asarray(maxworkers=maxworkers)
        # in place, so no second plane is allocated: raise values below the offset to it, then subtract
        offset = self.channel_offsets[c]
        if subtract_offset and offset:
            np.maximum(plane, offset, out=plane)
            plane -= np.asarray(offset, dtype=plane.dtype)
        return plane


def _qptiff_channel_name(description):
    """Marker name of one qptiff page: <Biomarker>, else the filter <Name>."""
    # the page description is XML; take the first non-empty tag, Biomarker preferred
    for tag in ("Biomarker", "Name"):
        m = re.search(fr"<{tag}>(.*?)</{tag}>", description or "")
        if m and m.group(1).strip():
            return m.group(1).strip()
    return None


def _ome_root(ome_xml):
    # parse the OME-XML string into an element tree
    return ET.fromstring(ome_xml)


def _local(tag):
    # strip the XML namespace: "{http://...}Channel" -> "Channel"
    return tag.rsplit("}", 1)[-1]


def _ome_channel_names(ome_xml):
    """Channel names of the first OME image (namespace-agnostic)."""
    image = next(e for e in _ome_root(ome_xml).iter() if _local(e.tag) == "Image")
    # one <Channel Name=...> element per channel, in channel order
    return [e.get("Name") for e in image.iter() if _local(e.tag) == "Channel"]


def _ome_channel_offsets(ome_xml):
    """Pixel value offset of every channel of the first OME image, 0 where none is given.

    Xenium writes it in each channel's MapAnnotation (<M K="Pixel value offset">): 100 for the
    protein channels, 0 for DAPI and the morphology stains.
    """
    root = _ome_root(ome_xml)
    # offset of every MapAnnotation that declares one, by annotation ID
    offsets = {}
    for ann in root.iter():
        if _local(ann.tag) != "MapAnnotation":
            continue
        for m in ann.iter():
            if _local(m.tag) == "M" and m.get("K") == "Pixel value offset":
                offsets[ann.get("ID")] = int(float(m.text))
    # each channel points to its annotation(s) through <AnnotationRef ID=...>
    image = next(e for e in root.iter() if _local(e.tag) == "Image")
    channels = [e for e in image.iter() if _local(e.tag) == "Channel"]
    refs = [[r.get("ID") for r in ch if _local(r.tag) == "AnnotationRef"] for ch in channels]
    return [next((offsets[i] for i in ids if i in offsets), 0) for ids in refs]


# conversion of OME length units to µm
_UNIT_TO_UM = {"µm": 1.0, "um": 1.0, "nm": 1e-3, "mm": 1e3, "cm": 1e4, "m": 1e6}


def _ome_pixel_size_um(ome_xml):
    """PhysicalSizeX of the first OME image in µm, or None if absent."""
    pixels = next((e for e in _ome_root(ome_xml).iter() if _local(e.tag) == "Pixels"), None)
    if pixels is None or pixels.get("PhysicalSizeX") is None:
        return None
    # convert to µm (OME's default unit when none is given)
    unit = pixels.get("PhysicalSizeXUnit", "µm")
    return float(pixels.get("PhysicalSizeX")) * _UNIT_TO_UM.get(unit, 1.0)


def _tiff_resolution_um(page):
    """Pixel size in µm from the TIFF XResolution/ResolutionUnit tags, or None if absent."""
    res, unit = page.tags.get("XResolution"), page.tags.get("ResolutionUnit")
    if res is None or unit is None:
        return None
    # XResolution is a rational (numerator, denominator) in pixels per ResolutionUnit
    num, den = res.value
    if not num or not den:
        return None
    pixels_per_unit = num / den
    # ResolutionUnit: 2 = inch, 3 = centimetre
    um_per_unit = {2: 25400.0, 3: 10000.0}.get(int(unit.value))
    return um_per_unit / pixels_per_unit if um_per_unit else None


def write_ome_tiff(path, planes, shape, dtype, channel_names, pixel_size_um, tile=(512, 512)):
    """Write a CYX (or YX) OME-TIFF plane by plane, tiled and zlib-compressed.

    planes is an iterable of (Y, X) arrays, one per channel, so only one plane is in memory.
    """
    axes = "CYX" if len(shape) == 3 else "YX"

    # split each plane into tiles in row-major order, as tifffile expects from an iterator
    def tiles():
        for plane in planes:
            for y in range(0, plane.shape[0], tile[0]):
                for x in range(0, plane.shape[1], tile[1]):
                    block = plane[y : y + tile[0], x : x + tile[1]]
                    # edge tiles are zero-padded to the full tile size
                    if block.shape != tile:
                        block = np.pad(block, ((0, tile[0] - block.shape[0]), (0, tile[1] - block.shape[1])))
                    yield block

    # stream the tiles into a BigTIFF with pixel size both as TIFF resolution (px/cm) and in OME metadata
    with tifffile.TiffWriter(path, bigtiff=True, ome=True) as tif:
        tif.write(
            tiles(),
            shape=shape,
            dtype=dtype,
            tile=tile,
            compression="zlib",
            photometric="minisblack",
            resolution=(1e4 / pixel_size_um, 1e4 / pixel_size_um),
            resolutionunit="CENTIMETER",
            metadata={
                "axes": axes,
                "PhysicalSizeX": pixel_size_um,
                "PhysicalSizeXUnit": "µm",
                "PhysicalSizeY": pixel_size_um,
                "PhysicalSizeYUnit": "µm",
                "Channel": {"Name": list(channel_names)},
            },
        )
