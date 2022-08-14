import mmap
import time
import re
import fasta
import gzip


def make_record(index, max_size):
    seqid = str(d_mmap[record_ind2_extended[index][0]:record_ind2_extended[index][1]], 'utf-8').split()[0][1:]
    desc = str(d_mmap[record_ind2_extended[index][0] + len(seqid) + 1:record_ind2_extended[index][1]], 'utf-8').strip()
    if index == max_size - 1:
        mmap_size = len(d_mmap)
        seq = str(d_mmap[record_ind2_extended[index][1]:mmap_size], 'utf-8').replace('\n', '')
        return fasta.Record(seqid, seq, desc)

    seq = str(d_mmap[record_ind2_extended[index][1]:record_ind2_extended[index + 1][0]], 'utf-8').replace('\n', '')
    return fasta.Record(seqid, seq, desc)


fasta_re2 = re.compile(rb">.*")
tic = time.perf_counter()
with open('X:/edwards2016/GCA_005890655.1_Cmin_1.0_genomic.fna.gz', mode="rt") as file_obj:
    with mmap.mmap(file_obj.fileno(), length=0, access=mmap.ACCESS_READ) as mmap_obj:
        # in case of compression:
        d_mmap = gzip.decompress(mmap_obj)
        record_ind2_extended = [m.span() for m in re.finditer(fasta_re2, d_mmap[::])]
        # print(str(mmap_obj[0:20]))
        # r = [fasta.Record(id=str(mmap_obj[record_ind2_extended[elem][0]:record_ind2_extended[elem][1]], 'utf-8'), seq=str(mmap_obj[record_ind2_extended[elem][1]:record_ind2_extended[elem + 1][0]], 'utf-8')) for elem in range(len(record_ind2_extended) - 1)]
        rec_size = len(record_ind2_extended)
        r = [make_record(elem, rec_size) for elem in range(rec_size)]
        # r = [make_record(elem, rec_size) for elem in range(rec_size)]
toc = time.perf_counter()
elapsed_time = toc - tic
print(f"[fasta-parser] elapsed time: {elapsed_time:0.8f} seconds")
print(r[0].id)
