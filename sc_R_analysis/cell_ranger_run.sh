#!/bin/sh

module purge
module load CellRanger/3.0.2-bcl2fastq-2.20.0

sample_id=$1

instruction="cellranger count --id=$sample_id --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-3.0.0 --fastqs=/camp/stp/babs/inputs/sequencing/fastq/190510_K00102_0338_BH77LFBBXY/fastq_tenx,/camp/stp/babs/inputs/sequencing/fastq/190513_K00102_0340_BH77MKBBXY/fastq_tenx --sample=$sample_id --project=SC19093"

echo $instruction

cellranger count --id=$sample_id --transcriptome=/camp/svc/reference/Genomics/10x/10x_transcriptomes/refdata-cellranger-GRCh38-3.0.0 --fastqs=/camp/stp/babs/inputs/sequencing/fastq/190510_K00102_0338_BH77LFBBXY/fastq_tenx,/camp/stp/babs/inputs/sequencing/fastq/190513_K00102_0340_BH77MKBBXY/fastq_tenx --sample=$sample_id --project=SC19093
