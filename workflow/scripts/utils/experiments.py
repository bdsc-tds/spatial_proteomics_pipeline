"""Read, validate and scan the experiments.yml sample sheet.

Layout (as in xenium_analysis_pipeline): tissue -> panel -> donor -> [samples]. Keys starting with
"_" are metadata, not layers:
  - top level:   _base_path (raw proteomics images), _xenium_base_path (optional, raw Xenium runs)
  - panel level: technology-specific settings, see PANEL_TEMPLATE (required) and PANEL_OPTIONAL
  - donor:       list of sample names

Nothing per sample is stored; both paths are derived when the sheet is read:
  - image:  the single *.{_image_ext} file matching {_image_glob} (default "*", i.e. directly) in
            {_base_path}/{tissue}/{panel}/{donor}/{sample}/
  - xenium: the Xenium run {_xenium_base_path}/{tissue}/{panel}/{donor}/{sample}/ (layout of
            xenium_analysis_pipeline) whose donor name contains the proteomics donor and whose
            donor or sample name contains the proteomics sample; None if there is no such run or
            several (see match_xenium_runs)

This module only depends on the standard library (plus PyYAML in load_experiments) so that it can be
imported both by the Snakefile and by update_experiments.py.
"""

import glob
import os
import sys

# panel-level metadata written for every new panel; None values must be filled in by hand
PANEL_TEMPLATE = {
    "_technology": None,  # e.g. codex, comet
    "_image_ext": None,  # image extension inside each sample folder, e.g. qptiff, ome.tiff
    "_pixel_size_um": None,  # fallback when the image metadata carries no pixel size
    "_nuclear_channels": None,  # channels combined (max) into the nuclear signal, usually [DAPI]; the first is used for alignment
    "_membrane_channels": None,  # channels combined (max) into the membrane/cytoplasm signal for Mesmer
    "_background_subtraction": False,  # run background subtraction for this panel (COMET only)
}

# panel-level metadata that may be left out, with its default
PANEL_OPTIONAL = {
    # image location relative to the sample folder, as a glob; e.g. morphology_focus/ch0000_*.ome.tif
    # for a Xenium output folder, whose first morphology_focus file opens the whole multi-file image
    "_image_glob": "*",
}

# panel metadata holding channel lists (a single channel may be written as a plain string)
CHANNEL_LIST_KEYS = ("nuclear_channels", "membrane_channels")

# extensions tried for panels that are not in experiments.yml yet
DEFAULT_IMAGE_EXTS = ("qptiff", "ome.tiff", "ome.tif")


def is_meta(key):
    """Whether an experiments.yml key is metadata rather than a layer.

    key: any dict key of the experiments tree.
    Returns True when the key starts with "_".
    """
    return str(key).startswith("_")


def layer_items(node):
    """Child layers of one node of the experiments tree, metadata keys left out.

    node: dict of one layer (whole sheet, tissue or panel), or None.
    Returns a list of (name, child node) pairs.
    """
    return [(k, v) for k, v in (node or {}).items() if not is_meta(k)]


def donor_samples(donor_node):
    """Sample names of one donor.

    donor_node: value under a donor key: a list of names, a single name, or None.
    Returns a list of sample names (empty for None).
    """
    # normalize the three accepted forms (None, single string, list) into a list
    if donor_node is None:
        return []
    if isinstance(donor_node, str):
        return [donor_node]
    return list(donor_node)


def iter_samples(experiments):
    """Walk every sample of the sheet.

    experiments: parsed experiments.yml dict.
    Yields (tissue, panel, donor, sample, panel_node), panel_node being the panel's dict (with its
    "_" metadata).
    """
    for tissue, tissue_node in layer_items(experiments):
        for panel, panel_node in layer_items(tissue_node):
            for donor, donor_node in layer_items(panel_node):
                for sample in donor_samples(donor_node):
                    yield tissue, panel, donor, sample, panel_node


def resolve_path(path, experiments_path=None):
    """Make a path from experiments.yml absolute.

    path: path as written in the sheet (absolute, relative, or None).
    experiments_path: path of experiments.yml; relative paths are resolved against its folder.
    Returns the absolute path, or path unchanged when it is None, already absolute, or
    experiments_path is None.
    """
    if path is None or os.path.isabs(path) or experiments_path is None:
        return path
    # realpath so that a symlinked experiments.yml resolves against its real folder
    return os.path.normpath(os.path.join(os.path.dirname(os.path.realpath(experiments_path)), path))


def load_experiments(experiments_path):
    """Read experiments.yml.

    experiments_path: path of experiments.yml.
    Returns the parsed dict ({} for an empty file).
    """
    # imported here so the rest of the module stays standard-library only
    import yaml

    with open(experiments_path) as f:
        return yaml.safe_load(f) or {}


def find_images(sample_dir, exts, pattern=PANEL_OPTIONAL["_image_glob"]):
    """Image files of one sample folder.

    sample_dir: sample folder path.
    exts: accepted extensions without the dot, e.g. ["qptiff"].
    pattern: glob relative to sample_dir the files must match (default "*": directly in the folder).
    Returns the sorted paths, relative to sample_dir, ending in any of exts; [] if the folder is missing.
    """
    if not os.path.isdir(sample_dir):
        return []
    # regular files only (no subfolders) matching the pattern whose name ends in one of the extensions
    return sorted(
        os.path.relpath(f, sample_dir)
        for f in glob.glob(os.path.join(glob.escape(sample_dir), pattern))
        if os.path.isfile(f) and any(f.endswith("." + e) for e in exts)
    )


def get_samples(experiments, experiments_path=None):
    """Validate experiments.yml and flatten it into one entry per sample.

    experiments: parsed experiments.yml dict.
    experiments_path: path of experiments.yml, to resolve relative _base_path/_xenium_base_path.
    Returns {"tissue/panel/donor/sample": info}, info holding the four layer names, "image" (absolute
    image path), "xenium" (absolute Xenium run dir, None without a unique match) and "panel_meta"
    (panel metadata without the leading "_"). Raises ValueError listing every incomplete panel or
    sample folder without exactly one image; warns about ambiguous Xenium matches.
    """
    # raw data roots; the proteomics one is mandatory, the Xenium one optional
    base_path = resolve_path(experiments.get("_base_path"), experiments_path)
    xenium_base_path = resolve_path(experiments.get("_xenium_base_path"), experiments_path)
    if base_path is None:
        raise ValueError("experiments.yml: _base_path is not set")

    # index the Xenium runs once
    xenium_runs = index_xenium_runs(xenium_base_path) if xenium_base_path else []

    samples, errors, ambiguous = {}, [], []
    for tissue, panel, donor, sample, panel_node in iter_samples(experiments):
        sample_id = f"{tissue}/{panel}/{donor}/{sample}"

        # panel metadata: every template key must be present and non-null
        missing = [k for k in PANEL_TEMPLATE if panel_node.get(k) is None]
        if missing:
            errors.append(f"panel {tissue}/{panel}: missing or null {', '.join(missing)}")
            continue
        # metadata keys without the leading "_", e.g. panel_meta["image_ext"], optional keys defaulted
        panel_meta = {k.lstrip("_"): v for k, v in {**PANEL_OPTIONAL, **panel_node}.items() if is_meta(k)}
        # channel lists given as a single name become one-element lists
        for key in CHANNEL_LIST_KEYS:
            if isinstance(panel_meta[key], str):
                panel_meta[key] = [panel_meta[key]]

        # image: exactly one file with the panel's extension (and glob) in the sample folder
        sample_dir = os.path.join(base_path, tissue, panel, donor, sample)
        images = find_images(sample_dir, [panel_meta["image_ext"]], panel_meta["image_glob"])
        if len(images) != 1:
            errors.append(f"sample {sample_id}: expected one {panel_meta['image_glob']} *.{panel_meta['image_ext']} "
                          f"in {sample_dir}, found {len(images)}")
            continue

        # Xenium run: donor, then sample name match; ambiguous matches are left out
        matches = match_xenium_runs(donor, sample, xenium_runs)
        if len(matches) > 1:
            ambiguous.append(f"{sample_id}: {', '.join(matches)}")

        # flat per-sample entry with absolute paths
        samples[sample_id] = {
            "tissue": tissue,
            "panel": panel,
            "donor": donor,
            "sample": sample,
            "image": os.path.join(sample_dir, images[0]),
            "xenium": os.path.join(xenium_base_path, matches[0]) if len(matches) == 1 else None,
            "panel_meta": panel_meta,
        }

    # report every problem at once rather than one per run
    if errors:
        raise ValueError("experiments.yml is incomplete:\n  " + "\n  ".join(sorted(set(errors))))
    if ambiguous:
        print("WARNING: several Xenium runs match these samples, alignment skipped:\n  "
              + "\n  ".join(ambiguous), file=sys.stderr)
    return samples


def get_panel_exts(experiments):
    """Image extension of every panel that already defines one.

    experiments: parsed experiments.yml dict.
    Returns {panel: _image_ext}, across all tissues.
    """
    exts = {}
    # a panel name used in several tissues keeps the extension of the last one read
    for _, tissue_node in layer_items(experiments):
        for panel, panel_node in layer_items(tissue_node):
            if panel_node and panel_node.get("_image_ext"):
                exts[panel] = panel_node["_image_ext"]
    return exts


def get_panel_globs(experiments):
    """Image glob of every panel that defines one (_image_glob).

    experiments: parsed experiments.yml dict.
    Returns {panel: _image_glob}, across all tissues.
    """
    globs = {}
    # a panel name used in several tissues keeps the glob of the last one read
    for _, tissue_node in layer_items(experiments):
        for panel, panel_node in layer_items(tissue_node):
            if panel_node and panel_node.get("_image_glob"):
                globs[panel] = panel_node["_image_glob"]
    return globs


def scan_raw_samples(base_path, panel_exts=None, default_exts=DEFAULT_IMAGE_EXTS, panel_globs=None):
    """Find the sample folders on disk, {base_path}/{tissue}/{panel}/{donor}/{sample}/.

    base_path: raw proteomics directory (_base_path).
    panel_exts: {panel: extension} from get_panel_exts; panels not in it use default_exts.
    default_exts: extensions tried for panels not yet in experiments.yml.
    panel_globs: {panel: image glob} from get_panel_globs; panels not in it look directly in the folder.
    Returns {(tissue, panel, donor, sample): [image paths relative to the sample folder]}; the list is
    empty for sample folders that hold no image yet (e.g. CODEX samples that only have their crop geojson).
    """
    panel_exts = panel_exts or {}
    panel_globs = panel_globs or {}
    found = {}
    # every directory exactly four levels below base_path is a sample folder
    for sample_dir in sorted(glob.glob(os.path.join(base_path, "*", "*", "*", "*", ""))):
        key = tuple(os.path.relpath(sample_dir, base_path).split(os.sep))
        # skip folders whose path has a "_"-prefixed level, which would clash with metadata keys
        if any(is_meta(k) for k in key):
            continue
        # images: files with the panel's extension, or any default extension for unknown panels
        exts = [panel_exts[key[1]]] if key[1] in panel_exts else list(default_exts)
        found[key] = find_images(sample_dir, exts, panel_globs.get(key[1], PANEL_OPTIONAL["_image_glob"]))
    return found


def contains_name(name, query):
    """Whether one name contains another, compared across modalities.

    name: name searched in, e.g. a Xenium donor or run folder name.
    query: name searched for, e.g. a proteomics donor or sample name.
    Returns True when query appears in name as whole "_"-separated tokens, ignoring case and
    reading "-" as "_": donor_01 is in donor_01 and in
    output-XETG..__donor_01__.., but not in donor_011.
    """
    # normalize case and separators
    name, query = (str(s).lower().replace("-", "_") for s in (name, query))
    # pad both with "_" so only whole tokens match
    return f"_{query}_" in f"_{name}_"


def index_xenium_runs(xenium_base_path):
    """List the Xenium runs on disk, {xenium_base_path}/{tissue}/{panel}/{donor}/{sample}/.

    xenium_base_path: raw Xenium directory (_xenium_base_path).
    Returns a sorted list of (tissue, panel, donor, sample) tuples, one per run.
    """
    # every directory exactly four levels below the base path is a run, as laid out by xenium_analysis_pipeline
    run_dirs = glob.glob(os.path.join(xenium_base_path, "*", "*", "*", "*", ""))
    return sorted(tuple(os.path.relpath(d, xenium_base_path).split(os.sep)) for d in run_dirs)


def match_xenium_runs(donor, sample, xenium_runs):
    """Find the Xenium run(s) of one proteomics sample.

    donor, sample: proteomics donor and sample names.
    xenium_runs: (tissue, panel, donor, sample) tuples from index_xenium_runs.
    Returns the sorted run dirs, relative to _xenium_base_path, whose donor name contains donor and
    whose donor or sample (run folder) name contains sample (contains_name). Tissue and panel are
    ignored; more than one result means the sample is ambiguous.
    """
    # runs of the same donor
    donor_runs = [run for run in xenium_runs if contains_name(run[2], donor)]
    # among them, runs of the same sample
    return sorted(
        os.path.join(*run) for run in donor_runs if contains_name(run[2], sample) or contains_name(run[3], sample)
    )


def unlisted_samples(experiments, experiments_path=None):
    """Samples on disk that experiments.yml does not list yet.

    experiments: parsed experiments.yml dict.
    experiments_path: path of experiments.yml, to resolve a relative _base_path.
    Returns "tissue/panel/donor/sample" IDs whose raw folder holds exactly one image but that are
    missing from the sheet; [] when _base_path is unset or missing.
    """
    base_path = resolve_path(experiments.get("_base_path"), experiments_path)
    if base_path is None or not os.path.isdir(base_path):
        return []
    # compare the samples in the sheet with the sample folders on disk
    listed = {(t, p, d, s) for t, p, d, s, _ in iter_samples(experiments)}
    scanned = scan_raw_samples(base_path, get_panel_exts(experiments), panel_globs=get_panel_globs(experiments))
    return ["/".join(k) for k, images in scanned.items() if len(images) == 1 and k not in listed]
