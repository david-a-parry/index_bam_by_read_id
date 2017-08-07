#!/usr/bin/env python
import pysam
import pickle
from collections import OrderedDict
import os

_by_qname = lambda x: x.query_name

class IndexByReadId(object):
    ''' 
        Class for sorting, indexing and retrieving reads from BAMs 
        by Read ID
    '''

    def __init__(self, bam, index=None):
        ''' 
            Initialization requires a bam and optionally a filename
            for the index (defaults to bam + '.ibbr').
        '''
        if index is None:
            self.index_file = bam + '.ibbr'
        else:
            self.index_file = index
        self.bam = bam
        self.idx = None
        self.k_idx = None
        self._cache = []
        self.bamfile = pysam.AlignmentFile(bam, 'rb')
        
    def sort_bam(self, outfile=None, batch_size=2000000):
        ''' 
            Samtools sort does not appear to guarantee ascibetical or
            numeric (including hexidecimal) order of field in read 
            names. This method sorts ascibetically on full read name.
            If outfile is not given, output will be the name of self.bam
            without the extension + '_rid_sorted.bam' (i.e. in.bam 
            becomes in_rid_sorted.bam).
        '''
        if outfile is None:
            (f, ext) = os.path.splitext(self.bam)
            outfile = f + '_rid_sorted.bam'
        header = self.bamfile.header
        if 'HD' not in header:
            header['HD'] = dict()
        header['HD']['SO'] = 'queryname'
        sink = pysam.AlignmentFile(outfile, 'wb', header=header)
        mergers = []
        source = pysam.AlignmentFile(self.bam, 'rb')
        n = 0
        m = 0 
        recs = []
        for r in source.fetch(until_eof=True):
            n += 1
            recs.append(r)
            if n % batch_size == 0:
                merge_fn = self._merge_write(recs, outfile, header, m)
                recs[:] = []
                mergers.append(pysam.AlignmentFile(merge_fn, 'rb'))
                m += 1
        if recs:
            merge_fn = self._merge_write(recs, outfile, header, m)
            recs[:] = []
            mergers.append(pysam.AlignmentFile(merge_fn, 'rb'))
        source.close()
        # merge onto sink
        stack_tops = [next(f) for f in mergers]
        #import pdb; pdb.set_trace() 
        while stack_tops:
            c = min(stack_tops, key=_by_qname)
            sink.write(c)
            i = stack_tops.index(c)
            try:
                t = next(mergers[i])
                stack_tops[i] = t
            except StopIteration:
                del stack_tops[i]
                mergers[i].close()
                os.remove(mergers[i].filename.decode())
                del mergers[i]  # __del__ method of file_opener should delete the file
        sink.close()

    def _merge_write(self, recs, outfile, header, n):
        recs.sort(key=_by_qname)
        merge_fn = outfile + str.format(".{:05d}.ibbr.tmp.bam", n)
        merge_bam = pysam.AlignmentFile(merge_fn, 'wb', header=header)
        for mrec in recs:
            merge_bam.write(mrec)
        merge_bam.close()
        return merge_fn
     
    def create_index(self, chunk_size=50000):
        ''' Create an index for a BAM sorted by Read ID '''
        idx = OrderedDict()
        n = 0
        prev_qname = ''
        with pysam.AlignmentFile( self.bam, "rb" ) as b_fh:
            p = b_fh.tell()
            for r in b_fh.fetch(until_eof=True):
                if prev_qname and r.qname < prev_qname:
                    raise UnsortedBamError("Input ({}) is not sorted by Read ID"
                                    .format(self.bam))
                if n % chunk_size == 0:
                    idx[r.qname] = p
                n += 1
                p = b_fh.tell()
                prev_qname = r.qname
        if prev_qname:
            idx[prev_qname] = p
        pfh = open(self.index_file, 'wb')
        pickle.dump(idx, pfh)
        pfh.close()

    def read_index(self):
        ''' Read an index file created using the create_index function.'''
        pfh = open(self.index_file, 'rb')
        self.idx = pickle.load(pfh)
        self.k_idx = list(self.idx.keys())
        pfh.close()

    def get_reads_by_id(self, rid):
        '''
            Retrieve reads with matching ID. Returns a list of 
            pysam.AlignedSegment objects (or an empty list if no reads
            are found). The index must exist.
        '''
        if self.idx is None:
            self.read_index()
        reads = self._get_matching_reads(rid)
        return reads

    def _get_matching_reads(self, rid):
        if rid < self.k_idx[0]:
            #first read in index is gt rid
            return []
        if rid > self.k_idx[-1]:
            #last read in index is lt rid
            return []
        (before, after) = self._get_nearest_indices(rid)
        return self._match_within_indices(rid, before, after)

    def _match_within_indices(self, rid, l, u):
        matches = []
        if l == u: #exact match
            #look at reads before matching index
            first = self.k_idx[l-1]
            last = self.k_idx[l]
            self._set_cache(first, last)
            for i in range(len(self._cache)-2, -1, -1):
                if self._cache[i].query_name == rid:
                    matches.append(self._cache[i])
                else:
                    break
            #look at reads including and after matching index
            first = self.k_idx[l]
            last = self.k_idx[l+1]
            self._set_cache(first, last)
            for i in range(0, len(self._cache)):
                if self._cache[i].query_name == rid:
                    matches.append(self._cache[i])
                else:
                    break
        else: #match is somewhere inbetween self.k_idx[l] and self.k_idx[u]
            matches = self._search_reads(rid, l, u)
        return matches

    def _search_reads(self, rid, l, u):
        first = self.k_idx[l]
        last = self.k_idx[u]
        self._set_cache(first, last)
        i = self._binsearch_cache(rid)
        matches = []
        if i < 0:
            return []
        for j in range(i-1, -1, -1):
            if self._cache[j].query_name == rid:
                matches.append(self._cache[j])
            else:
                break
        for j in range(i, len(self._cache)):
            if self._cache[j].query_name == rid:
                matches.append(self._cache[j])
            else:
                break
        return matches

    def _binsearch_cache(self, x):
        l = 0
        u = len(self._cache)
        while l <= u:
            i = int((u + l)/2)
            q = self._cache[i].query_name
            if x < q:
                u = i - 1
            elif x > q:
                l = i + 1
            else:
                return i
        return -1 #not found

    def _set_cache(self, id1, id2):
        if (self._cache and self._cache[0].query_name == id1
            and self._cache[-1].query_name == id2):
            return
        else:
            self._cache = self._reads_inbetween(id1, id2)

    def _reads_inbetween(self, id1, id2):
        reads = []
        start = self.idx[id1]
        stop  = self.idx[id2]
        self.bamfile.seek(start)  
        p = self.bamfile.tell()
        while p <= stop:
            try:
                reads.append(next(self.bamfile))
                p = self.bamfile.tell()
            except StopIteration:
               break 
        return reads
        
    def _read_from_index(self, i):
        pass

    def _get_nearest_indices(self, rid):
        '''
            Binary search of indexed read names (we have already checked
            that rid is within the range of self.k_idx)
        '''
        l = 0
        u = len(self.k_idx) - 1
        while l <= u:
            i = int((u + l)/2)
            kid = self.k_idx[i]
            if rid < kid: 
                pre = self.k_idx[i-1]
                if rid > pre: #rid is gt prev in index
                    return (i-1, i)
                u = i - 1
            elif rid > kid: 
                nxt = self.k_idx[i+1]
                if rid < nxt: #rid is lt next in index
                    return (i, i+1)
                l = i + 1
            else: #exact match
                return (i, i)
        #not found (shouldn't happen after check in _get_matching_reads)
        return (None, None) 


class UnsortedBamError(Exception):
    pass

