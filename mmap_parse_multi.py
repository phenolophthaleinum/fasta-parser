import itertools
import mmap
import time
import re
import regex
import fasta
import gzip
import numpy as np
from joblib import Parallel, delayed
from multiprocessing import cpu_count, Pool
import threading


# def make_record(index, max_size):
#     seqid = str(d_mmap[record_ind2_extended[index][0]:record_ind2_extended[index][1]], 'utf-8').split()[0][1:]
#     desc = str(d_mmap[record_ind2_extended[index][0] + len(seqid) + 1:record_ind2_extended[index][1]], 'utf-8').strip()
#     if index == max_size - 1:
#         mmap_size = len(d_mmap)
#         seq = str(d_mmap[record_ind2_extended[index][1]:mmap_size], 'utf-8').replace('\n', '')
#         return fasta.Record(seqid, seq, desc)
#
#     seq = str(d_mmap[record_ind2_extended[index][1]:record_ind2_extended[index + 1][0]], 'utf-8').replace('\n', '')
#     return fasta.Record(seqid, seq, desc)


def make_record_fasta(index, max_size):
    seqid = str(mmap_obj[record_ind2_extended[index][0]:record_ind2_extended[index][1]], 'utf-8').split()[0][1:]
    desc = str(mmap_obj[record_ind2_extended[index][0] + len(seqid) + 1:record_ind2_extended[index][1]],
               'utf-8').strip()
    if index == max_size - 1:
        mmap_size = len(mmap_obj)
        seq = str(mmap_obj[record_ind2_extended[index][1]:mmap_size], 'utf-8').replace('\n', '')
        return fasta.Record(seqid, seq, desc)

    seq = str(mmap_obj[record_ind2_extended[index][1]:record_ind2_extended[index + 1][0]], 'utf-8').replace('\n', '')
    return fasta.Record(seqid, seq, desc)


def partition_queries(records, partitions: int = cpu_count() - 1) -> np.ndarray:
    partitions = np.array_split(np.array(records), partitions)
    # debug for partitions
    # print(partitions)
    # print([ar.shape for ar in partitions])
    return partitions


def batch_exec(part, max_size):
    # local_list = []
    # for elem in part:
    #     local_list.extend(make_record_fasta(elem, max_size))
    # return local_list
    return [make_record_fasta(elem, max_size) for elem in part]


# def regex_search(start, end):
#     return [m.span() for m in re.finditer(fasta_re2, mmap_obj[start:end])]
    # record_ind2_extended.extend([m.span() for m in re.finditer(fasta_re2, mmap_obj[start:end])])
    # return [m.span() for m in re.finditer(fasta_re2, mmap_obj[slice[0]:slice[1]])]


# record_ind2_extended = []
fasta_re2 = re.compile(rb">.*")
tic = time.perf_counter()
# X:/edwards2016/GCA_005890655.1_Cmin_1.0_genomic.fna.gz
# X:/edwards2016/models/stats_test/random_train-species-dim_10-len_125.fasta
# X:/edwards2016/models/wrapped/random_train-genus-dim_10-len_125.fasta
with open('X:/edwards2016/models/stats_test/random_train-species-dim_10-len_125.fasta', mode="rt") as file_obj:
    with mmap.mmap(file_obj.fileno(), length=0, access=mmap.ACCESS_READ) as mmap_obj:
        file_size = len(mmap_obj)
        file_part_sizes = [file_size // cpu_count() + (1 if x < file_size % cpu_count() else 0) for x in
                           range(cpu_count())]
        # file_parts = [0]
        # file_parts.extend(list(itertools.accumulate(file_part_sizes)))
        # in case of compression:
        # d_mmap = gzip.decompress(mmap_obj)

        # normal parallel
        record_ind2_extended = [m.span() for m in re.finditer(fasta_re2, mmap_obj[::])]

        # paralleled regex
        # threads = []
        # for i in range(len(file_parts) - 1):
        #     t = threading.Thread(target=regex_search, args=(file_parts[i], file_parts[i + 1]))
        #     threads.append(t)
        #     t.start()
        # for index, thread in enumerate(threads):
        #     thread.join()

        # r_ind = Parallel(prefer='threads', verbose=True, n_jobs=-1)(delayed(regex_search)(file_parts[i], file_parts[i + 1]) for i in range(len(file_parts) - 1))
        # record_ind2_extended = []
        # for e in r_ind:
        #     record_ind2_extended.extend(e)

        # pool
        # file_slices = [(file_parts[i], file_parts[i + 1]) for i in range(len(file_parts) - 1)]
        # print(file_slices[0:2])
        # with Pool() as pool:
        #     record_ind2_extended = pool.imap(regex_search, file_slices)
        #
        # print(list(record_ind2_extended)[0])
        # print(str(mmap_obj[0:20]))
        # r = [fasta.Record(id=str(mmap_obj[record_ind2_extended[elem][0]:record_ind2_extended[elem][1]], 'utf-8'), seq=str(mmap_obj[record_ind2_extended[elem][1]:record_ind2_extended[elem + 1][0]], 'utf-8')) for elem in range(len(record_ind2_extended) - 1)]
        rec_size = len(record_ind2_extended)
        # r = [make_record_fasta(elem, rec_size) for elem in range(rec_size)]
        # r = [make_record(elem, rec_size) for elem in range(rec_size)]
        r_temp = Parallel(prefer='threads', verbose=False, n_jobs=-1)(delayed(batch_exec)(batch, rec_size) for batch in
                                                                     partition_queries(range(rec_size)))
        # r_temp = Parallel(prefer='threads', n_jobs=-1)(delayed(make_record_fasta)(elem, rec_size) for elem in range(rec_size))
        r = []
        for e in r_temp:
            r.extend(e)

toc = time.perf_counter()
elapsed_time = toc - tic
print(f"[fasta-parser] elapsed time: {elapsed_time:0.8f} seconds")
print(r[0].id)
print(len(r))
