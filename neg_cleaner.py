#!/usr/bin/env python3
"""
Remove contaminação de controle negativo em leituras Nanopore (fastq).

Para cada placa, identifica as leituras associadas aos poços de controle
negativo (arquivos .fa cujo nome contém "Neg"), extrai uma subsequência de
assinatura (posições 50-100) de cada leitura negativa e remove, dos demais
arquivos fastq da mesma placa, qualquer leitura cuja sequência contenha essa
assinatura com distância de edição <= 1 (indício de contaminação/index
hopping a partir do controle negativo).
"""

import argparse
import gzip
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

import edlib
from Bio import SeqIO
from tqdm import tqdm

SUBSEQ_START = 50
SUBSEQ_END = 100
MAX_EDIT_DISTANCE = 1


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--demux-folder",
        required=True,
        help="Pasta com os .fa demultiplexados por placa (saída do ONTbarcoder).",
    )
    parser.add_argument(
        "--fastq-folder",
        required=True,
        help="Pasta com os arquivos .fastq.gz brutos a serem limpos.",
    )
    parser.add_argument(
        "--output-dict",
        default=None,
        help="Arquivo onde salvar o dicionário de IDs/subsequências por placa "
        "(padrão: dic_<timestamp>.txt dentro de --fastq-folder).",
    )
    parser.add_argument(
        "--output-suffix",
        default="_clean",
        help="Sufixo adicionado ao nome dos arquivos fastq limpos (padrão: _clean).",
    )
    return parser.parse_args()


def load_plate_ids(demux_folder):
    """Carrega, por placa, os IDs de todas as leituras e dos negativos."""
    placas_data = {}
    id_neg = []

    for fasta in os.listdir(demux_folder):
        if not fasta.endswith(".fa"):
            continue
        placa_clean = fasta.split("_")[2]
        if placa_clean not in placas_data:
            placas_data[placa_clean] = {"ids": [], "ids_neg": [], "subseqs": set()}

        for record in SeqIO.parse(os.path.join(demux_folder, fasta), "fasta"):
            placas_data[placa_clean]["ids"].append(record.id)
            if "Neg" in fasta:
                placas_data[placa_clean]["ids_neg"].append(record.id)
                id_neg.append(record.id)

    return placas_data, id_neg


def extract_negative_subseqs(fastq_folder, placas_data, id_neg):
    """Extrai a subsequência de assinatura de cada leitura de controle negativo."""
    for arqfastq in tqdm(
        os.listdir(fastq_folder), desc="Extraindo subsequências dos negativos", unit="arquivo"
    ):
        if not arqfastq.endswith(".fastq.gz"):
            continue
        file_path = os.path.join(fastq_folder, arqfastq)
        with gzip.open(file_path, "rt") as handle:
            for record in SeqIO.parse(handle, "fastq"):
                if record.id not in id_neg:
                    continue
                subseq = str(record.seq)[SUBSEQ_START:SUBSEQ_END]
                for data in placas_data.values():
                    if record.id in data["ids_neg"]:
                        data["subseqs"].add(subseq)
                        break


def save_dict(output_dict, placas_data):
    with open(output_dict, mode="w", encoding="utf-8") as f:
        for placa_clean, data in placas_data.items():
            f.write(f"Placa: {placa_clean}\n")
            f.write(f"  IDs: {data['ids']}\n")
            f.write(f"  IDs Negativos: {data['ids_neg']}\n")
            f.write(f"  Subsequências: {list(data['subseqs'])}\n\n")


def clean_file(fastq, fastq_folder, output_suffix, placas_data):
    file_path = os.path.join(fastq_folder, fastq)
    file_prefix = fastq.replace(".fastq.gz", "")
    output_file = os.path.join(fastq_folder, f"{file_prefix}{output_suffix}.fastq.gz")

    with gzip.open(file_path, "rt") as in_handle, gzip.open(output_file, "wt") as out_handle:
        for record in SeqIO.parse(in_handle, "fastq"):
            data = next(
                (d for d in placas_data.values() if record.id in d["ids"]), None
            )
            if data is None:
                SeqIO.write(record, out_handle, "fastq")
                continue

            seq_fastq = str(record.seq)
            match_found = any(
                edlib.align(s, seq_fastq, mode="HW", k=MAX_EDIT_DISTANCE)["editDistance"]
                in (0, MAX_EDIT_DISTANCE)
                for s in data["subseqs"]
            )
            if not match_found:
                SeqIO.write(record, out_handle, "fastq")

    print(f"[{datetime.now().strftime('%d/%m/%Y %H:%M:%S')}] Arquivo processado: {fastq}")


def main():
    args = parse_args()

    output_dict = args.output_dict or os.path.join(
        args.fastq_folder, f"dic_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    )

    print("Carregando IDs por placa...")
    placas_data, id_neg = load_plate_ids(args.demux_folder)

    extract_negative_subseqs(args.fastq_folder, placas_data, id_neg)

    save_dict(output_dict, placas_data)
    print(f"Dicionário salvo em {output_dict}")

    fastqs = [f for f in os.listdir(args.fastq_folder) if f.endswith(".fastq.gz")]
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(clean_file, fastq, args.fastq_folder, args.output_suffix, placas_data)
            for fastq in fastqs
        ]
        for _ in tqdm(as_completed(futures), total=len(futures), desc="Progresso"):
            pass

    print("Todos os arquivos foram processados com sucesso.")


if __name__ == "__main__":
    main()
