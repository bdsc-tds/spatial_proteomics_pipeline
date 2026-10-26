"""ESM-2 embedding of every protein of a VirTues marker table: the identity VirTues gives each channel.

Reads the panel's marker table (channel_name, uniprot_id; channels without an ID are not used), downloads
each protein sequence from UniProt (isoforms like P08575-8 included) and embeds it with ESM-2
(esm2_t30_150M_UR50D, last layer averaged over the residues, as VirTues' compute_esm_embeddings.py).
Writes one .pt file: {"uniprot_ids": [...], "embeddings": (n_proteins, 640) tensor}.
"""

import argparse
import os
import urllib.request

import pandas as pd
import torch

ESM_MODEL = "esm2_t30_150M_UR50D"
ESM_LAYER = 30


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--markers", required=True, help="marker table of the panel (channel_name, uniprot_id)")
    parser.add_argument("--model_dir", required=True, help="folder holding (or receiving) the ESM-2 weights")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--out_pt", required=True)
    return parser.parse_args()


def fetch_sequence(uniprot_id):
    """Protein sequence of a UniProt accession or isoform, from the UniProt REST API."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    with urllib.request.urlopen(url, timeout=60) as response:
        lines = response.read().decode().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"no FASTA for {uniprot_id} at {url}")
    return "".join(line.strip() for line in lines[1:])


def main():
    args = parse_args()
    torch.set_num_threads(args.threads)
    import esm

    # unique proteins of the panel, sorted (one embedding each, shared by channels of the same protein)
    markers = pd.read_csv(args.markers, dtype=str)
    uniprot_ids = sorted(markers["uniprot_id"].dropna().str.strip().unique())
    print(f"{len(uniprot_ids)} proteins in {args.markers}", flush=True)

    # ESM-2 on CPU; its weights are cached under model_dir (torch hub) and downloaded on first use
    torch.hub.set_dir(os.path.join(args.model_dir, "torch_hub"))
    model, alphabet = esm.pretrained.load_model_and_alphabet(ESM_MODEL)
    model.eval()
    batch_converter = alphabet.get_batch_converter()

    # last-layer representation averaged over the residues (BOS / EOS tokens excluded); no contact
    # prediction, which VirTues' script also runs but which does not change the representation and
    # needs (length²) memory for long proteins such as Ki-67
    embeddings = []
    for uniprot_id in uniprot_ids:
        sequence = fetch_sequence(uniprot_id)
        _, _, tokens = batch_converter([(uniprot_id, sequence)])
        with torch.no_grad():
            representation = model(tokens, repr_layers=[ESM_LAYER])["representations"][ESM_LAYER]
        embeddings.append(representation[0, 1 : len(sequence) + 1].mean(0))
        print(f"{uniprot_id}: {len(sequence)} residues", flush=True)

    torch.save({"uniprot_ids": uniprot_ids, "embeddings": torch.stack(embeddings)}, args.out_pt)
    print(f"written: {args.out_pt}")


if __name__ == "__main__":
    main()
