# neg-cleaner

Ferramenta para remover contaminação de controle negativo em leituras de
sequenciamento Nanopore (metabarcoding), usada no fluxo de trabalho do
BIODOSSEL/BIOINSECTA.

## O problema

Em placas de sequenciamento, poços de controle negativo (sem amostra) às
vezes acabam recebendo leituras por contaminação cruzada ou *index hopping*.
Se essas sequências "vazarem" para os poços com amostra real, elas podem ser
confundidas com resultado verdadeiro.

## Como funciona

Para cada placa:

1. Carrega os IDs de todas as leituras (a partir dos `.fa` demultiplexados
   pelo ONTbarcoder) e identifica quais pertencem a poços de controle
   negativo (arquivos cujo nome contém `Neg`).
2. Extrai, de cada leitura de controle negativo, uma subsequência de
   assinatura (posições 50–100 da sequência).
3. Nos arquivos `.fastq.gz` brutos da mesma placa, remove qualquer leitura
   cuja sequência contenha essa assinatura com distância de edição ≤ 1
   (indicando que veio do controle negativo).
4. Salva um dicionário (`dic_<timestamp>.txt`) com os IDs e assinaturas por
   placa, e grava os arquivos fastq limpos com o sufixo configurado
   (padrão `_clean`).

O processamento dos arquivos fastq é paralelizado (`ProcessPoolExecutor`).

## Instalação

```bash
pip install -r requirements.txt
```

## Uso

```bash
python neg_cleaner.py \
  --demux-folder "caminho/para/Demultiplexed" \
  --fastq-folder "caminho/para/fastq_pass" \
  [--output-dict "caminho/para/dic.txt"] \
  [--output-suffix "_clean"]
```

- `--demux-folder`: pasta com os `.fa` demultiplexados por placa (saída do
  ONTbarcoder). Os nomes dos arquivos devem seguir o padrão
  `algo_algo_<placa>_...`, e os arquivos de controle negativo devem conter
  `Neg` no nome.
- `--fastq-folder`: pasta com os `.fastq.gz` brutos a serem limpos. Os
  arquivos de saída são gravados nessa mesma pasta.
- `--output-dict` (opcional): caminho do dicionário de IDs/subsequências.
  Se omitido, é salvo como `dic_<timestamp>.txt` dentro de `--fastq-folder`.
- `--output-suffix` (opcional, padrão `_clean`): sufixo dos arquivos fastq
  limpos gerados.

## Gerar executável (Windows)

Para rodar em máquinas sem Python instalado, dá para gerar um `.exe`
standalone com o [PyInstaller](https://pyinstaller.org/):

```bash
pip install -r requirements.txt pyinstaller
pyinstaller --onefile --name neg_cleaner --console neg_cleaner.py
```

O executável fica em `dist/neg_cleaner.exe` e aceita os mesmos argumentos do
script (`--demux-folder`, `--fastq-folder` etc.):

```bash
dist\neg_cleaner.exe --demux-folder "caminho\Demultiplexed" --fastq-folder "caminho\fastq_pass"
```

## Histórico

As três versões anteriores deste script (evoluídas entre setembro e
novembro de 2025) foram consolidadas num único script parametrizável. As
versões originais continuam disponíveis no histórico do git
(commit inicial deste repositório).
