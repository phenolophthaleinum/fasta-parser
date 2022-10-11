import mmap
import time
import re
import fastq
from memory_profiler import profile
import array
import gzip
from io import BytesIO


# @profile
def main():
    # def make_record_fasta(index, max_size):
    #     seqid = str(mmap_obj[record_ind2_extended[index][0]:record_ind2_extended[index][1]], 'utf-8').split()[0][1:]
    #     desc = str(mmap_obj[record_ind2_extended[index][0] + len(seqid) + 1:record_ind2_extended[index][1]],
    #                'utf-8').strip()
    #     if index == max_size - 1:
    #         mmap_size = len(mmap_obj)
    #         seq = str(mmap_obj[record_ind2_extended[index][1]:mmap_size], 'utf-8').replace('\n', '')
    #         return fasta.Record(seqid, seq, desc)
    #
    #     seq = str(mmap_obj[record_ind2_extended[index][1]:record_ind2_extended[index + 1][0]], 'utf-8').replace('\n',
    #                                                                                                             '')
    #     return fasta.Record(seqid, seq, desc)
    def make_record_fastq(index):
        seqid = str(mmap_obj[offsets[index]:offsets[index + 1]], 'utf-8').split()[0][1:]
        desc = str(mmap_obj[offsets[index]:offsets[index + 1]], 'utf-8')[1:]
        seq = str(mmap_obj[offsets[index + 1]:offsets[index + 2]], 'utf-8')[:-1]
        quality = str(mmap_obj[offsets[index + 3]:offsets[index + 4]], 'utf-8')[:-1]

        return fastq.Record(seqid, seq, quality, desc)

    tic = time.perf_counter()
    # X:/edwards2016/GCA_005890655.1_Cmin_1.0_genomic.fna.gz
    # X:/edwards2016/models/stats_test/random_train-species-dim_10-len_125.fasta
    # X:/edwards2016/models/wrapped/random_train-genus-dim_10-len_125.fasta
    with gzip.open('C:/Users/Maciej/PycharmProjects/fasta2png/fastqs/ERR6099228.fastq', mode="r+b") as file_obj:
        # file_obj = gzip.decompress(file_obj)
    # with fi.input('X:/edwards2016/GCA_005890655.1_Cmin_1.0_genomic.fna.gz', openhook=fi.hook_compressed) as file_obj:
        # with mmap its easier to lookup specific parts of file
        with mmap.mmap(file_obj.fileno(), length=0, access=mmap.ACCESS_READ) as mmap_obj:
            # line_map = iter(mmap_obj.readline, b"")
            # print(type(mmap_obj))

            # if compressed
            # mmap_obj = gzip.decompress(mmap_obj)
            # mmap_obj_bio = BytesIO(mmap_obj)

            offsets = [0]

            # if compressed
            # lines = [mmap_obj_bio.tell() for _ in enumerate(iter(mmap_obj_bio.readline, b""))]

            lines = [mmap_obj.tell() for _ in enumerate(iter(mmap_obj.readline, b""))]
            offsets.extend(lines)
            offsets_size = len(offsets)
            # print(offsets[0], offsets[1])
            # print(mmap_obj[offsets[0]:offsets[1]])
            # print(mmap_obj[:4])
            records = [make_record_fastq(i) for i in range(0, offsets_size - 1, 4)]
            # in case of compression:
            # mmap_obj = gzip.decompress(mmap_obj)
            # record_ind2_extended = [m.span() for m in re.finditer(fasta_re2, mmap_obj[::])]
            # print(str(mmap_obj[0:20]))
            # r = [fasta.Record(id=str(mmap_obj[record_ind2_extended[elem][0]:record_ind2_extended[elem][1]], 'utf-8'), seq=str(mmap_obj[record_ind2_extended[elem][1]:record_ind2_extended[elem + 1][0]], 'utf-8')) for elem in range(len(record_ind2_extended) - 1)]
            # rec_size = len(record_ind2_extended)
            # r = [make_record_fasta(elem, rec_size) for elem in range(rec_size)]
            # r = [make_record(elem, rec_size) for elem in range(rec_size)]
    # interesting memory release shown by memory_profiler
    # del mmap_obj
    # gc.collect()
    toc = time.perf_counter()
    elapsed_time = toc - tic
    print(f"[fasta-parser] elapsed time: {elapsed_time:0.8f} seconds")
    print(records[0])
    # print(r[0].id)
    # print(len(r))
    # print(records[0].id)
    # print(' '.join(records[0].desc.split()[1:]))
    # print(records[0].seq)
    # print(records[0].phred_quality)
    # print(len(records))


if __name__ == "__main__":
    main()