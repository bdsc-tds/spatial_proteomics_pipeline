"""Create or update experiments.yml from the raw data directory.

Run before the first pipeline execution, and again whenever new samples arrive:

    python workflow/scripts/utils/update_experiments.py --config config/configs_codex/config.yml --dry-run

- Scans {_base_path}/{tissue}/{panel}/{donor}/{sample}/ for sample folders holding one image and
  adds the ones missing from experiments.yml, as `donor: [sample, ...]`.
- Never overwrites values that are already set: hand-edited panel metadata and comments survive
  (ruamel.yaml round-trip).
- New panels get the PANEL_TEMPLATE metadata block (or a copy of --panel_template's); null values
  must be filled in by hand before the pipeline accepts the panel.
- Image and Xenium paths are not stored: the pipeline derives them when it reads the sheet. When
  _xenium_base_path is set, this script previews which Xenium run each sample will be matched to
  (same donor, then same sample name, see experiments.match_xenium_runs).
- Samples listed in experiments.yml but gone from disk are only reported, unless --prune.
"""

import argparse
import copy
import difflib
import io
import os
import sys

from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from experiments import (
    DEFAULT_IMAGE_EXTS,
    donor_samples,
    PANEL_TEMPLATE,
    get_panel_exts,
    index_xenium_runs,
    iter_samples,
    layer_items,
    match_xenium_runs,
    resolve_path,
    scan_raw_samples,
)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config", required=True, help="pipeline config.yml; its `experiments` key locates experiments.yml")
    parser.add_argument("--base_path", help="raw image directory; only used when experiments.yml has no _base_path yet")
    parser.add_argument("--xenium_base_path", help="raw Xenium directory; only used when experiments.yml has no _xenium_base_path yet")
    parser.add_argument("--image_ext", help=f"image extension for panels not yet in experiments.yml (default: any of {', '.join(DEFAULT_IMAGE_EXTS)})")
    parser.add_argument("--panel_template", help="existing panel whose metadata is copied into new panels")
    parser.add_argument("--prune", action="store_true", help="remove samples whose raw folder no longer holds an image")
    parser.add_argument("--dry-run", action="store_true", help="print the resulting changes without writing")
    return parser.parse_args()


def new_yaml():
    # round-trip loader that keeps quotes and comments, with the indentation used in experiments.yml
    yaml = YAML()
    yaml.preserve_quotes = True
    yaml.indent(mapping=2, sequence=4, offset=2)
    # wide enough that long lines are never wrapped
    yaml.width = 4096
    return yaml


def dump(yaml, data):
    # serialize to a string (used both for writing and for diffing old vs. new)
    buf = io.StringIO()
    yaml.dump(data, buf)
    return buf.getvalue()


def experiments_path_from_config(yaml, config_path):
    # read the pipeline config and check it points to an experiments sheet
    with open(config_path) as f:
        cfg = yaml.load(f)
    if not cfg or not cfg.get("experiments"):
        sys.exit(f"{config_path}: no `experiments` entry")
    # relative to the real config file, so it also works through the config/config.yml symlink
    return resolve_path(cfg["experiments"], os.path.realpath(config_path))


def find_panel_template(experiments, panel_name):
    # copy the _-prefixed metadata keys (not the donors) of the first panel with this name, in any tissue
    for _, tissue_node in layer_items(experiments):
        for panel, panel_node in layer_items(tissue_node):
            if panel == panel_name:
                return CommentedMap((k, copy.deepcopy(v)) for k, v in panel_node.items() if str(k).startswith("_"))
    sys.exit(f"--panel_template {panel_name}: panel not found in experiments.yml")


def flow_list(items):
    """List dumped on one line: [a, b]."""
    seq = CommentedSeq(items)
    seq.fa.set_flow_style()
    return seq


def child(node, key, factory=CommentedMap):
    # get node[key], creating an empty entry first if it is missing or null
    if node.get(key) is None:
        node[key] = factory()
    return node[key]


def main():
    args = parse_args()
    yaml = new_yaml()
    experiments_path = experiments_path_from_config(yaml, args.config)

    # load the existing sheet, or start a new one from the CLI paths
    if os.path.exists(experiments_path):
        with open(experiments_path) as f:
            experiments = yaml.load(f) or CommentedMap()
        old_text = dump(yaml, experiments)
    else:
        experiments, old_text = CommentedMap(), ""
    # CLI paths only fill in raw-data roots the sheet does not define yet
    for key, value in (("_base_path", args.base_path), ("_xenium_base_path", args.xenium_base_path)):
        if experiments.get(key) is None and value:
            experiments[key] = value
    # resolve the raw image root (relative paths are relative to experiments.yml)
    base_path = resolve_path(experiments.get("_base_path"), experiments_path)
    if not base_path or not os.path.isdir(base_path):
        sys.exit(f"raw directory not found: {base_path!r} (set _base_path or pass --base_path)")

    # scan raw sample folders
    default_exts = [args.image_ext] if args.image_ext else DEFAULT_IMAGE_EXTS
    panel_exts = get_panel_exts(experiments)
    scanned = scan_raw_samples(base_path, panel_exts, default_exts)
    # split folders by number of images found: exactly one (usable), none, several (ambiguous, skipped)
    with_image = {k: v[0] for k, v in scanned.items() if len(v) == 1}
    no_image = [k for k, v in scanned.items() if not v]
    multi_image = {k: v for k, v in scanned.items() if len(v) > 1}

    # new panels start from a template metadata block
    template = find_panel_template(experiments, args.panel_template) if args.panel_template else None

    # add samples missing from the sheet, as a one-line [sample, ...] list per donor
    added = []
    for (tissue, panel, donor, sample), image in with_image.items():
        tissue_node = child(experiments, tissue)
        # first sample of a panel in this tissue: create the panel with template metadata
        if tissue_node.get(panel) is None:
            panel_node = copy.deepcopy(template) if template else CommentedMap(PANEL_TEMPLATE)
            # the extension always comes from the image actually found, not from the template panel,
            # out of the extensions the scan used for this panel name (it may already be declared by
            # the same panel in another tissue)
            scan_exts = [panel_exts[panel]] if panel in panel_exts else default_exts
            panel_node["_image_ext"] = max((e for e in scan_exts if image.endswith("." + e)), key=len)
            tissue_node[panel] = panel_node
        panel_node = tissue_node[panel]
        # normalize the donor entry (missing, or a single sample string) into a one-line list
        if not isinstance(panel_node.get(donor), list):
            panel_node[donor] = flow_list(donor_samples(panel_node.get(donor)))
        # append the sample only if not already listed
        if sample not in panel_node[donor]:
            panel_node[donor].append(sample)
            added.append((tissue, panel, donor, sample))

    # samples listed in the sheet whose raw image is gone
    gone = [(t, p, d, s) for t, p, d, s, _ in iter_samples(experiments) if (t, p, d, s) not in with_image]
    # with --prune, drop them (and the donor entry once its list is empty)
    if args.prune:
        for tissue, panel, donor, sample in gone:
            donor_list = experiments[tissue][panel][donor]
            donor_list.remove(sample)
            if not donor_list:
                del experiments[tissue][panel][donor]

    # preview the Xenium runs the pipeline will match on donor and sample name (nothing is stored)
    xenium_filled, xenium_ambiguous, xenium_unmatched = [], {}, []
    xenium_base_path = resolve_path(experiments.get("_xenium_base_path"), experiments_path)
    if xenium_base_path:
        if not os.path.isdir(xenium_base_path):
            sys.exit(f"Xenium directory not found: {xenium_base_path!r}")
        runs = index_xenium_runs(xenium_base_path)
        # classify each sample as matched to one run, ambiguous (several runs) or unmatched
        for tissue, panel, donor, sample, _ in iter_samples(experiments):
            matches = match_xenium_runs(donor, sample, runs)
            if len(matches) == 1:
                xenium_filled.append(((tissue, panel, donor, sample), matches[0]))
            elif matches:
                xenium_ambiguous[(tissue, panel, donor, sample)] = matches
            else:
                xenium_unmatched.append((tissue, panel, donor, sample))

    # summary
    fmt = lambda k: "/".join(k)
    print(f"experiments.yml: {experiments_path}")
    print(f"raw sample folders: {len(scanned)} ({len(with_image)} with one image, {len(no_image)} without image, {len(multi_image)} with several)")
    print(f"added samples: {len(added)}")
    for k in added:
        print(f"  + {fmt(k)}")
    for k, images in multi_image.items():
        print(f"  several images, skipped: {fmt(k)}: {', '.join(images)}")
    if gone:
        print(f"listed but no image on disk: {len(gone)}" + (" (pruned)" if args.prune else " (use --prune to remove)"))
        for k in gone:
            print(f"  - {fmt(k)}")
    if xenium_base_path:
        print(f"Xenium runs matched: {len(xenium_filled)}, ambiguous (alignment skipped): {len(xenium_ambiguous)}, unmatched: {len(xenium_unmatched)}")
        for k, run in xenium_filled:
            print(f"  {fmt(k)} -> {run}")
        for k, matches in xenium_ambiguous.items():
            print(f"  ambiguous, alignment skipped: {fmt(k)}: {', '.join(matches)}")
        for k in xenium_unmatched:
            print(f"  no Xenium run for donor {k[2]!r}, sample {k[3]!r}")
    # panels whose template metadata still has null values
    incomplete = [f"{t}/{p}" for t, tissue_node in layer_items(experiments) for p, panel_node in layer_items(tissue_node)
                  if any(panel_node.get(key) is None for key in PANEL_TEMPLATE)]
    if incomplete:
        print("panels with null metadata (fill in before running the pipeline): " + ", ".join(incomplete))

    # write, or show the diff
    new_text = dump(yaml, experiments)
    if new_text == old_text:
        print("no changes")
    elif args.dry_run:
        sys.stdout.writelines(difflib.unified_diff(
            old_text.splitlines(True), new_text.splitlines(True), "experiments.yml", "experiments.yml (updated)"))
    else:
        os.makedirs(os.path.dirname(experiments_path), exist_ok=True)
        with open(experiments_path, "w") as f:
            f.write(new_text)
        print(f"written: {experiments_path}")


if __name__ == "__main__":
    main()
