from Bio import SeqIO
import fasta
import time
from memory_profiler import profile


@profile
def main():
    tic = time.perf_counter()
    # X:/edwards2016/GCA_005890655.1_Cmin_1.0_genomic.fna.gz
    # X:/edwards2016/models/stats_test/random_train-species-dim_10-len_125.fasta
    # /home/hyperscroll/edwards2016/Cmin_scaff.fna
    # X:/edwards2016/models/wrapped/random_train-genus-dim_10-len_125.fasta
    record1 = list(fasta.parse('X:/edwards2016/GCA_005890655.1_Cmin_1.0_genomic.fna.gz'))
    toc = time.perf_counter()
    elapsed_time = toc - tic
    print(f"[fasta-parser] elapsed time: {elapsed_time:0.8f} seconds")
    print(record1[0].id)
    print(len(record1))


if __name__ == "__main__":
    main()
