# neg-cleaner — Biodossel/Bioinsecta

O neg-cleaner é um programa que limpa a contaminação de controle negativo dos dados de sequenciamento (Nanopore/metabarcoding) usados no fluxo de trabalho do projeto Biodossel/Bioinsecta.

## Qual problema ele resolve

Em cada placa de sequenciamento, alguns poços são reservados como "controle negativo" — não deveriam conter DNA de nenhuma amostra real. Às vezes, por contaminação cruzada durante o processo, sequências que na verdade vieram desses poços de controle acabam aparecendo, por engano, nos poços com amostra. Isso pode fazer parecer que uma amostra teve um resultado que, na verdade, veio do controle negativo.

O neg-cleaner identifica esse "vazamento" e remove essas leituras contaminadas dos arquivos brutos de sequenciamento, gerando uma versão limpa de cada arquivo.

## O que você precisa ter em mãos

Duas pastas:

1. **Pasta com os arquivos demultiplexados por placa**, gerados pelo ONTbarcoder (arquivos terminados em `.fa`). Os arquivos dos poços de controle negativo precisam ter **"Neg" no nome** — é assim que o programa reconhece quais são os controles. Além disso, o nome de todos os arquivos dessa pasta precisa seguir o padrão usado pelo ONTbarcoder (algo assim: `qualquercoisa_qualquercoisa_NOMEDAPLACA_...`), pois é dali que o programa descobre a qual placa cada arquivo pertence.
2. **Pasta com os arquivos brutos do sequenciamento** (`.fastq.gz`), os que serão limpos.

## Instalação

Requer [Python](https://www.python.org/downloads/) 3.8 ou mais recente instalado no computador.

1. Baixe/copie a pasta do projeto (com `neg_cleaner.py` e `requirements.txt`).
2. Abra o **Prompt de Comando** (cmd) ou PowerShell nessa pasta: segure Shift, clique com o botão direito num espaço vazio dentro da pasta e escolha "Abrir janela do PowerShell aqui" (ou "Abrir Prompt de Comando aqui").
3. Instale as dependências do programa (só precisa fazer isso uma vez):

   ```
   pip install -r requirements.txt
   ```

## Como usar

1. Na mesma janela de terminal, digite o comando abaixo, trocando os caminhos pelos das suas pastas (mantenha as aspas):

   ```
   python neg_cleaner.py --demux-folder "C:\caminho\para\Demultiplexed" --fastq-folder "C:\caminho\para\fastq_pass"
   ```

2. Aperte Enter e aguarde. O programa mostra o andamento na tela (quantos arquivos já foram processados).
3. Quando terminar, os arquivos limpos aparecem **na mesma pasta dos arquivos brutos**, com o mesmo nome seguido de `_clean` (ex: `barcode01.fastq.gz` vira `barcode01_clean.fastq.gz`). Os arquivos originais não são apagados nem alterados.
4. Também é salvo, nessa mesma pasta, um arquivo de texto (`dic_<data_hora>.txt`) com o registro de quais leituras foram identificadas como controle negativo em cada placa — serve como conferência, não é necessário abri-lo no uso normal.

### Ajustes opcionais

Só use se precisar mudar o comportamento padrão:

- `--output-dict "C:\caminho\dic.txt"` — escolhe onde salvar o arquivo de registro (`dic_...txt`), em vez do padrão.
- `--output-suffix "_limpo"` — troca o sufixo `_clean` dos arquivos gerados por outro texto de sua escolha.

---

## Para quem for mexer no código-fonte

### Estrutura do projeto

```
neg-cleaner/
├── neg_cleaner.py     # script único: parsing de args, pipeline completo
└── requirements.txt
```

O script é mantido como um arquivo único (não dividido em módulos) — as etapas do pipeline (carregar IDs, extrair assinaturas, salvar dicionário, limpar fastqs) já são funções bem separadas dentro dele.

### Como o pareamento entre negativo e leitura funciona

Para cada leitura de controle negativo encontrada nos fastqs brutos, o programa extrai a subsequência das posições 50–100 e associa essa "assinatura" à placa correspondente. Em seguida, para cada leitura de cada fastq, verifica se ela pertence a uma placa conhecida e, se sim, se sua sequência contém alguma assinatura de negativo daquela placa (distância de edição ≤ 1); em caso positivo, a leitura é descartada do arquivo limpo. O processamento dos arquivos fastq é paralelizado (`ProcessPoolExecutor`).

## Histórico

As três versões anteriores deste script (evoluídas entre setembro e novembro de 2025) foram consolidadas num único script parametrizável. As versões originais continuam disponíveis no histórico do git (commit inicial deste repositório).
