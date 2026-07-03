


import os
import gzip
from Bio import SeqIO
import edlib
from concurrent.futures import ProcessPoolExecutor, as_completed
from tqdm import tqdm
from datetime import datetime

# ================================
# CONFIGURAÇÕES GLOBAIS
# ================================
input_folder = 'C:/Users/biodo/OneDrive/Documentos/Outputs_de_cada_Sequenciamento/BIMSeq004_480_14jul25/basecalling/Demultiplexed'
input_folder_fastq = 'C:/Users/biodo/OneDrive/Documentos/Outputs_de_cada_Sequenciamento/BIMSeq004_480_14jul25/basecalling/pass/'
dic = 'C:/Users/biodo/OneDrive/Documentos/Vscode/dic.txt'


def clean_file(fastq, placas_data):
    file_path = os.path.join(input_folder_fastq, fastq)
    file_prefix = fastq.replace(".fastq.gz", "")
    output_file = os.path.join(input_folder_fastq, f"{file_prefix}_sem_query.fastq.gz")

    with gzip.open(file_path, "rt") as in_handle, gzip.open(output_file, "wt") as out_handle:
        for record in SeqIO.parse(in_handle, "fastq"):
            id_found = False
            for placa_clean, data in placas_data.items():
                if record.id in data["ids"]:
                    id_found = True
                    seq_fastq = str(record.seq)
                    match_found = any(
                        edlib.align(s, seq_fastq, mode="HW")["editDistance"] in (0, 1)
                        for s in data["subseqs"]
                    )
                    if not match_found:
                        SeqIO.write(record, out_handle, "fastq")
                    break
            if not id_found:
                SeqIO.write(record, out_handle, "fastq")

    print(f"[{datetime.now().strftime('%d/%m/%Y %H:%M:%S')}] ✅ Arquivo processado: {fastq}")


if __name__ == "__main__":
    placas_data = {}
    id_neg = []
    countnaoescritas = 0

    # ================================
    # ETAPA 1: CARREGAR IDS
    # ================================
    for fasta in os.listdir(input_folder):
        placa = fasta.split('_')
        placa_clean = placa[2]
        if placa_clean not in placas_data:
            placas_data[placa_clean] = {"ids": [], "ids_neg": [], "subseqs": set()}
        if fasta.endswith(".fa") and placa_clean in fasta:
            for record in SeqIO.parse(os.path.join(input_folder, fasta), "fasta"):
                placas_data[placa_clean]["ids"].append(record.id)
                if "Neg" in fasta:
                    placas_data[placa_clean]["ids_neg"].append(record.id)
                    id_neg.append(record.id)

    # ================================
    # ETAPA 2: EXTRAÇÃO DE SUBSEQUÊNCIAS
    # ================================
    for arqfastq in tqdm(os.listdir(input_folder_fastq), desc="Processando arquivos FASTQ", unit="arquivo"):
        if arqfastq.endswith(".fastq.gz"):
            file_path = os.path.join(input_folder_fastq, arqfastq)
            with gzip.open(file_path, "rt") as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    if record.id in id_neg:
                        seq_fastq = str(record.seq)
                        subseq = seq_fastq[50:100]
                        for placa_clean, data in placas_data.items():
                            if record.id in data["ids_neg"]:
                                data["subseqs"].add(subseq)
                                break

    # ================================
    # ETAPA 3: SALVAR DICIONÁRIO
    # ================================
    with open(dic, mode='a', encoding='utf-8') as f:
        for placa_clean, data in placas_data.items():
            f.write(f"Placa: {placa_clean}\n")
            f.write(f"  IDs: {data['ids']}\n")
            f.write(f"  IDs Negativos: {data['ids_neg']}\n")
            f.write(f"  Subsequências: {list(data['subseqs'])}\n\n")

    print("Processamento finalizado! Resultados salvos em", dic)

    # ================================
    # ETAPA 4: LIMPAR ARQUIVOS FASTQ
    # ================================
    fastqs = [f for f in os.listdir(input_folder_fastq) if f.endswith(".fastq.gz")]
    with ProcessPoolExecutor() as executor:
        futures = []
        for fastq in fastqs:
            futures.append(executor.submit(clean_file, fastq, placas_data))

        for _ in tqdm(as_completed(futures), total=len(futures), desc="🔄 Progresso"):
            pass

    print("✅ Todos os arquivos foram processados com sucesso.")
    print(f"nao escritas: {countnaoescritas}")
