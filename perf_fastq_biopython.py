from Bio import SeqIO
import time
from memory_profiler import profile
import pyfastx


# @profile
def main():
    tic = time.perf_counter()
    # biopython
    # with open("C:/Users/Maciej/PycharmProjects/fasta2png/fastqs/ERR6099228.fastq") as fq:
    #     all_rec = [record for record in SeqIO.parse(fq, 'fastq')]
    # toc = time.perf_counter()
    # elapsed_time = toc - tic
    # print(f"[biopython-parser] elapsed time: {elapsed_time:0.8f} seconds")
    # print(all_rec[0])
    # print(len(all_rec))

    # pyfastx
    all_rec = [r for r in pyfastx.Fasta("X:/edwards2016/models/wrapped/random_train-genus-dim_10-len_125.fasta", build_index=False)]
    toc = time.perf_counter()
    elapsed_time = toc - tic
    print(all_rec[0])
    print(f"[biopython-parser] elapsed time: {elapsed_time:0.8f} seconds")

if __name__ == "__main__":
    main()