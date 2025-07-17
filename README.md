# GenomeAnonymizer

GenomeAnonymizer is a tool for anonymizing genomic data by processing BAM/CRAM files and applying various privacy-preserving transformations. It is designed for researchers and organizations who need to share or analyze genomic data while protecting sensitive information.

## Features

- Anonymizes genomic alignments in BAM/CRAM files
- Supports multiple anonymization strategies
- Flexible configuration via command-line arguments
- Output in standard formats for downstream analysis

## Requirements

- Java 21 or higher
- Gradle (for building from source)
- [htsjdk](https://samtools.github.io/htsjdk/) library (included as dependency)

## Installation

Clone the repository and build the project using Gradle:

```bash
git clone https://github.com/Computational-Genomics-BSC/GenomeAnonymizer.git
cd GenomeAnonymizer
./gradlew build
```

The executable JAR will be located in `build/libs/`.

## Usage

You can run GenomeAnonymizer using the following command:

```bash
java -jar build/libs/GenomeAnonymizer.jar [arguments]
```

### Arguments

| Argument                | Description                                                                                 | Required | Example                                 |
|-------------------------|---------------------------------------------------------------------------------------------|----------|-----------------------------------------|
| `-i`, `--input`         | Input BAM/CRAM file                                                                         | Yes      | `-i input.bam`                          |
| `-o`, `--output`        | Output BAM/CRAM file                                                                        | Yes      | `-o anonymized.bam`                     |
| `-r`, `--reference`     | Reference FASTA file (required for CRAM input/output)                                       | No       | `-r reference.fa`                       |
| `-m`, `--mode`          | Anonymization mode (`mask`, `shuffle`, `remove`, etc.)                                      | Yes      | `-m mask`                               |
| `-s`, `--seed`          | Random seed for reproducibility                                                             | No       | `-s 42`                                 |
| `-t`, `--threads`       | Number of threads to use                                                                    | No       | `-t 4`                                  |
| `--regions`             | BED file with regions to anonymize                                                          | No       | `--regions regions.bed`                 |
| `--exclude`             | BED file with regions to exclude from anonymization                                         | No       | `--exclude exclude.bed`                 |
| `--min-mapq`            | Minimum mapping quality for reads to be considered                                          | No       | `--min-mapq 20`                         |
| `--include-duplicates`  | Include duplicate reads in anonymization                                                    | No       | `--include-duplicates`                  |
| `--help`                | Show help message                                                                           | No       | `--help`                                |

> **Note:** For a full list of arguments and their descriptions, run:
> ```bash
> java -jar build/libs/GenomeAnonymizer.jar --help
> ```

### Example Commands

**Basic anonymization:**
```bash
java -jar build/libs/GenomeAnonymizer.jar -i input.bam -o anonymized.bam -m mask
```

**Anonymize with a specific region and random seed:**
```bash
java -jar build/libs/GenomeAnonymizer.jar -i input.bam -o anonymized.bam -m shuffle --regions regions.bed -s 1234
```

**Anonymize a CRAM file with reference:**
```bash
java -jar build/libs/GenomeAnonymizer.jar -i input.cram -o anonymized.cram -r reference.fa -m remove
```

## Output

The output file will be a BAM or CRAM file with the specified anonymization applied. The format matches the input file unless otherwise specified.

## Citation

If you use GenomeAnonymizer in your research, please cite the repository and the associated publication (if available).

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

---

For questions or support, please open an issue on the [GitHub repository](https://github.com/Computational-Genomics-BSC/GenomeAnonymizer).