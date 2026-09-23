# -*- coding: utf-8 -*-
import argparse
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from interval3 import Interval
from interval3 import IntervalSet
from collections import Counter
from multiprocessing import Process
from multiprocessing import Pool
from Bio.AlignIO import MultipleSeqAlignment
import gc
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Process MAF files and extract common alignment columns."
    )

    parser.add_argument(
        "--input", type=str, nargs='+', required=True,
        help="Input MAF files (space-separated)."
    )

    parser.add_argument(
        "--output", type=str, required=True,
        help="Path to the output MAF file."
    )

    parser.add_argument(
        "--global_num", type=int, required=True,
        help="Number of segments to divide the alignment (subtracted by 1)."
    )

    parser.add_argument(
        "--global_start", type=int, required=True,
        help="Start coordinate for segmentation."
    )

    parser.add_argument(
        "--global_end", type=int, required=True,
        help="End coordinate for segmentation."
    )

    parser.add_argument(
        "--separate_workers", type=int, default=4,
        help="Number of worker processes for the separate() stage (default: 4)."
    )

    parser.add_argument(
        "--common_workers", type=int, default=15,
        help="Number of worker processes for the get_common_colomn() stage (default: 15)."
    )

    return parser.parse_args()


def separate(start, gap, end, list0, num):
    # When only one segment is requested, return the whole interval directly.
    # This also avoids referring to an undefined loop variable `i`.
    if num == 0:
        real_interval = IntervalSet()
        real_interval.add(Interval(start, end, upper_closed=False))
        return [real_interval], [list0]

    lsts = [[] for i in range(num)]
    real_intervals = [IntervalSet() for i in range(num)]

    last_boundary = start
    for i, lst, real_interval in zip(range(start + gap, end, gap), lsts, real_intervals):
        real_interval.add(Interval(i - gap, i, upper_closed=False))
        last_boundary = i
        for aln in list0[:]:
            inter = Interval(
                aln[0].annotations["start"],
                aln[0].annotations["start"] + aln[0].annotations["size"],
                upper_closed=False
            )
            if i in inter:
                lst.append(aln)
            if aln[0].annotations["start"] + aln[0].annotations["size"] <= i:
                lst.append(aln)
                list0.remove(aln)
            if aln[0].annotations["start"] >= i:
                break

    lsts.append(list0)
    real_intervals.append(IntervalSet())
    real_intervals[-1].add(Interval(last_boundary, end, upper_closed=False))

    del list0
    gc.collect()

    return real_intervals, lsts


def get_common_colomn(lists, intervalsets_used):
    blocks = []
    con_lst = []
    ks = [0 for _ in range(len(lists))]

    for used_interval in list(intervalsets_used):
        for i in range(used_interval.lower_bound, used_interval.upper_bound):
            lst_infos = [[] for _ in range(len(lists))]
            for lst, k, lst_info in zip(lists, range(len(ks)), lst_infos):
                for aln in lst[:]:
                    aln_start = aln[0].annotations["start"]
                    aln_end = aln[0].annotations["start"] + aln[0].annotations["size"]
                    aln_interval = Interval(aln_start, aln_end, upper_closed=False)

                    if i in aln_interval:
                        if i == used_interval.lower_bound:
                            ks[k] = 0
                            pos1 = aln_start
                            while pos1 <= i:
                                if aln[0].seq[ks[k]] != '-':
                                    pos1 += 1
                                ks[k] += 1
                            ks[k] -= 1
                        else:
                            if i == aln_start:
                                ks[k] = 0
                                pos1 = aln_start
                                while pos1 <= i:
                                    if aln[0].seq[ks[k]] != '-':
                                        pos1 += 1
                                    ks[k] += 1
                                ks[k] -= 1
                            elif i > aln_start:
                                ks[k] += 1
                                nue = 0
                                while nue < 1:
                                    if aln[0].seq[ks[k]] != '-':
                                        nue += 1
                                    ks[k] += 1
                                ks[k] -= 1

                        for seqs in range(len(aln)):
                            lst_info.append([
                                aln[seqs].id,
                                aln[seqs].annotations["start"] + ks[k] - aln[seqs].seq[0:ks[k] + 1].count('-') + (1 if aln[seqs].seq[ks[k]] == '-' else 0),
                                0 if aln[seqs].seq[ks[k]] == '-' else 1,
                                aln[seqs].annotations["strand"],
                                aln[seqs].annotations["srcSize"],
                                aln[seqs].seq[ks[k]].upper()
                            ])
                    elif i >= aln_end:
                        lst.remove(aln)
                    elif i < aln_start:
                        break

            temp_con_lst = []
            for elem in lst_infos[0]:
                if all(elem in lst_info for lst_info in lst_infos[1:]):
                    temp_con_lst.append(elem)

            if len(temp_con_lst) == 1:
                temp_con_lst = []
            elif all(lst[-1] == '-' for lst in temp_con_lst[1:]):
                temp_con_lst = []

            if temp_con_lst != []:
                if con_lst == []:
                    con_lst = temp_con_lst
                    temp_con_lst = []
                elif len(temp_con_lst) == len(con_lst):
                    score = 0
                    for l in range(0, len(temp_con_lst)):
                        if (
                            temp_con_lst[l][0] == con_lst[l][0]
                            and temp_con_lst[l][1] == con_lst[l][1] + con_lst[l][2]
                            and temp_con_lst[l][3] == con_lst[l][3]
                        ):
                            score += 1
                    if score == len(temp_con_lst):
                        for l in range(len(temp_con_lst)):
                            con_lst[l][2] += temp_con_lst[l][2]
                            con_lst[l][5] += temp_con_lst[l][5]
                        temp_con_lst = []
                    else:
                        one_block = []
                        if len(con_lst) != 1:
                            for l in range(len(con_lst)):
                                one_block.append(
                                    SeqRecord(
                                        Seq(con_lst[l][5]),
                                        id=con_lst[l][0],
                                        annotations={
                                            'start': con_lst[l][1],
                                            'size': con_lst[l][2],
                                            'strand': con_lst[l][3],
                                            'srcSize': con_lst[l][4]
                                        }
                                    )
                                )
                            blocks.append(MultipleSeqAlignment(one_block))
                        con_lst = temp_con_lst
                else:
                    one_block = []
                    if len(con_lst) != 1:
                        for l in range(len(con_lst)):
                            one_block.append(
                                SeqRecord(
                                    Seq(con_lst[l][5]),
                                    id=con_lst[l][0],
                                    annotations={
                                        'start': con_lst[l][1],
                                        'size': con_lst[l][2],
                                        'strand': con_lst[l][3],
                                        'srcSize': con_lst[l][4]
                                    }
                                )
                            )
                        blocks.append(MultipleSeqAlignment(one_block))
                    con_lst = temp_con_lst
            else:
                if con_lst != []:
                    one_block = []
                    if len(con_lst) != 1:
                        for l in range(len(con_lst)):
                            one_block.append(
                                SeqRecord(
                                    Seq(con_lst[l][5]),
                                    id=con_lst[l][0],
                                    annotations={
                                        'start': con_lst[l][1],
                                        'size': con_lst[l][2],
                                        'strand': con_lst[l][3],
                                        'srcSize': con_lst[l][4]
                                    }
                                )
                            )
                        blocks.append(MultipleSeqAlignment(one_block))
                    con_lst = []

        if con_lst != []:
            one_block = []
            if len(con_lst) != 1:
                for l in range(len(con_lst)):
                    one_block.append(
                        SeqRecord(
                            Seq(con_lst[l][5]),
                            id=con_lst[l][0],
                            annotations={
                                'start': con_lst[l][1],
                                'size': con_lst[l][2],
                                'strand': con_lst[l][3],
                                'srcSize': con_lst[l][4]
                            }
                        )
                    )
                blocks.append(MultipleSeqAlignment(one_block))
            con_lst = []

    return blocks


def main():
    args = parse_args()

    input_files = args.input
    output_file = args.output
    global_num = args.global_num - 1
    global_start = args.global_start
    global_end = args.global_end
    separate_workers = args.separate_workers
    common_workers = args.common_workers

    lst_files = [list(AlignIO.parse(f, "maf")) for f in input_files]

    interval_length = global_end - global_start

    # Do not create more segments than there are reference positions.
    # This guarantees a positive step for range() even for very short intervals.
    requested_segments = global_num + 1
    actual_segments = min(requested_segments, interval_length)
    global_num = actual_segments - 1
    global_gap = max(1, round(interval_length / actual_segments))

    separate_workers = max(1, min(separate_workers, len(lst_files)))
    pool0 = Pool(processes=separate_workers)
    res0 = [
        pool0.apply_async(
            func=separate,
            args=(global_start, global_gap, global_end, lst, global_num)
        )
        for lst in lst_files
    ]
    pool0.close()
    pool0.join()

    result_lsts = []
    real_intervalss = []

    for resm in res0:
        result = resm.get()
        real_intervalss.append(result[0])
        result_lsts.append(result[1])

    num_common_tasks = len(real_intervalss[0])
    common_workers = max(1, min(common_workers, num_common_tasks))
    pool = Pool(processes=common_workers)
    res = [
        pool.apply_async(
            func=get_common_colomn,
            args=([result_lsts[j][i] for j in range(len(result_lsts))], list(real_intervalss[0][i]))
        )
        for i in range(num_common_tasks)
    ]
    pool.close()
    pool.join()

    record_blocks = []
    for resm in res:
        result = resm.get()
        record_blocks.extend(result)

    AlignIO.write(record_blocks, output_file, "maf")


if __name__ == "__main__":
    main()
