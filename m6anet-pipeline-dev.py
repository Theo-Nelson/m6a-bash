"""
python pipeline for m6a-net pipeline command generation

Author: Theo Nelson

Lab: Mason Lab
"""

input_path="/athena/masonlab/scratch/users/tmn2126/m6a_cortex/samples/"
output_path="/athena/masonlab/scratch/users/tmn2126/m6a_cortex/transcript_pipeline/"
ref_transcripts="/athena/masonlab/scratch/users/tmn2126/public_meta-analysis/GCF_000001635.27_GRCm39_rna.fna"
ref_gtf="/athena/masonlab/scratch/users/tmn2126/public_meta-analysis/GCF_000001635.27_GRCm39_genomic.gtf"


def minimap2(i):
    print("conda activate minimap2 ; ",end="")
    print("minimap2 -ax splice -t 16 -uf -k14 " + ref_transcripts + " " + input_path + i + "/" + i + ".fastq > " + output_path + i + "/" + i + ".sam ; ",end="")
    print("conda deactivate ; ",end ="")
    
def samtools(i):
    print("conda activate samtools ; ",end="")
    print("samtools view -S -b " + output_path + i + "/" + i + ".sam > " + output_path + i + "/" + i + ".bam ; ",end="")
    print("samtools sort " + output_path + i + "/" + i + ".bam -@ 8 -o " + output_path + i + "/" + i + ".sorted.bam ; ",end="")
    print("samtools index " + output_path + i + "/" + i + ".sorted.bam ; ",end="")
    print("samtools stats " + output_path + i + "/" + i + ".sorted.bam > " + output_path + i + "/" + i + "_alignment_stats.txt ; ",end="")
    print("conda deactivate ; ",end ="")

def featureCounts(i):
    print("conda activate featurecounts ; ",end="")
    print("featureCounts -O -L -a " + ref_gtf + " -t exon -g gene_id -o " + output_path + i + "/" + i + "_counts.txt " + output_path + i + "/" + i + ".sorted.bam ; ",end="")
    print("conda deactivate ; ",end ="")


def f5c(i):
    print("conda activate f5c ; ",end="")
    print("/athena/masonlab/scratch/users/tmn2126/m6a_cortex/samples/1338_KO/install_vbz.sh ; ",end="")
    print("export HDF5_PLUGIN_PATH=/home/tmn2126/.local/hdf5/lib/plugin ; ",end="")
    print("f5c eventalign --iop 600 --rna --scale-events --signal-index -b " + output_path + i + "/" + i + ".sorted.bam" + " -g " + ref_transcripts + " -r " + input_path + i + "/" + i + ".fastq" + " > " + output_path + i + "/" + i + "_eventalign.txt ; ",end="")
    print("conda deactivate ; ",end ="")

def m6anet(i):
    print("conda activate m6anet ; ",end="")
    print("m6anet dataprep --eventalign " + output_path + i + "/" + i + "_eventalign.txt " + "--out_dir " + output_path + i + "/dataprep --n_processes 600 ; ",end="")
    print("m6anet inference --input_dir " + output_path + i + "/dataprep --out_dir " + output_path + i + "/inference --n_processes 56 --num_iterations 1000 ; ",end="")
    print("conda deactivate ; ",end ="")

def generate_command(name):
    for i in name:
        #print("mkdir " + output_path + i + " ; ",end="")
        #minimap2(i)
        #samtools(i)
        #featureCounts(i)
        f5c(i)
        m6anet(i)

if __name__ == '__main__':
    generate_command(["1338_KO","1338_WT","9548_KO","9548_WT"])