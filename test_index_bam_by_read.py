from nose.tools import *
import os
import pysam
from collections import OrderedDict
try:
    from index_bam_by_read_id.index_bam_by_read_id  import IndexByReadId,UnsortedBamError
except ImportError:
    from index_bam_by_read_id  import IndexByReadId,UnsortedBamError

_f = 'test_data/test.bam'
_s = "test_data/test_rid_sorted.bam"
_sidx = "test_data/test_rid_sorted.bam.ibbr"
_f2 = 'test_data/sim_test.bam'
_s2 = 'test_data/sim_test_rid_sorted.bam'
_s2idx = 'test_data/sim_test_rid_sorted.bam.ibbr'
_scram = "test_data/test_rid_sorted.cram"
_ssam = "test_data/test_rid_sorted.sam"

def test_sort_input():
    for f in [_s, _s2, _sidx, _s2idx]:
        if os.path.exists(f):
            os.remove(f)
    ibbr = IndexByReadId(_f)
    ibbr.sort_bam(outfile=_s, batch_size=100)
    ibbr2 = IndexByReadId(_f2)
    ibbr2.sort_bam()
    assert_true(os.path.exists(_s2))
    ibbr.sort_bam(outfile=_scram, batch_size=100)
    ibbr.sort_bam(outfile=_ssam, batch_size=100)

def test_unsorted_error():
    ibbr = IndexByReadId(_f)
    assert_raises(UnsortedBamError, ibbr.create_index)

def test_create_index():
    if os.path.exists(_sidx):
        os.remove(_sidx)
        assert_false(os.path.exists(_sidx))
    ibbr = IndexByReadId(_s)
    ibbr.create_index(chunk_size=100)
    assert_true(os.path.exists(_sidx))
    ibbr2 = IndexByReadId(_s2)
    ibbr2.create_index()
    assert_true(os.path.exists(_s2idx))
    sibbr = IndexByReadId(_ssam)
    sibbr.create_index(chunk_size=100)
    cribbr = IndexByReadId(_scram)
    assert_raises(NotImplementedError, cribbr.create_index)
    
def test_read_index():
    assert_true(os.path.exists(_sidx))
    ibbr = IndexByReadId(_s)
    idx = ibbr.read_index()
    #import pdb; pdb.set_trace() 
    assert_true(isinstance(ibbr.idx, OrderedDict))
    assert_equal(ibbr.n_records, 1010)
    
def test_find_targets():
    targets = [ "SRR043396.16790502", "SRR043348.7107241", "SRR043348.14284382", 
                "SRR043354.12939592", "SRR043366.2208972", "SRR043372.1247100", 
                "SRR043372.11318637", "SRR043378.2924881", "SRR043378.13722892", 
                "SRR043386.7232182", "SRR043386.16375485", "SRR043396.6476248",  
                "SRR043354.3858155", ]
    targets2 = ["16_81676054_81676426_0_1_0_0_6:0:0_12:0:0_2089aa",
                "NC_003063.2_591540_591498_1_0_0_0_2:0:0_3:0:0_557e",
                "NC_005945.1_2248952_2249021_0_1_0_0_1:0:0_3:0:0_9b39",
                "NC_009715.2_846767_846681_1_0_0_0_5:0:0_4:0:0_d1db",
                "NC_013192.1_482003_481932_1_0_0_0_3:0:0_2:0:0_8231",
                "NC_015311.1_1039211_1039135_1_0_0_0_1:0:0_2:0:0_3431e",
                "NZ_ACDZ02000008.1_259070_258949_1_0_0_0_2:0:0_0:0:0_1c06",
                "NZ_ATUH01000004.1_46254_46186_1_0_0_0_2:0:0_6:0:0_295",
                "NZ_CP007062.1_455537_455567_0_1_0_0_0:0:0_2:0:0_1e538",
                "NZ_FOJP01000008.1_159948_159985_0_1_0_0_2:0:0_5:0:0_5767",
                "NZ_GL397214.1_371347_371440_0_1_0_0_1:0:0_5:0:0_32114",
                "NZ_GL892076.1_394940_395033_0_1_0_0_4:0:0_1:0:0_3b18",
                "NZ_JH806632.1_998649_998668_0_1_0_0_6:0:0_1:0:0_189da",
                "NZ_KB900367.1_231568_231461_1_0_0_0_3:0:0_2:0:0_20bb",
                "NZ_KE993508.1_78376_78338_1_0_0_0_3:0:0_2:0:0_d543",]
    for f, t, e in zip((_s, _ssam, _s2),
                    (targets, targets, targets2),
                    (2, 2, 1)):
        _find_targets(f, t, e)
    
def _find_targets(f, targets, expect):
    ibbr = IndexByReadId(f)
    ibbr.read_index()
    for t in targets:
        hits = ibbr.get_reads_by_id(t)
        n = 0
        for h in hits:
            n += 1
            assert_equal(t, h.query_name)
        assert_equal(n, expect)

def test_non_existent():
    ibbr = IndexByReadId(_s2)
    hits = ibbr.get_reads_by_id( #search for non-existent read
                        "NZ_GG703879.1_967999_968061_0_1_0_0_4:0:0_1:0:0_62bb")
    assert_equal(len(hits), 0)
    #search first and last records

def test_first_last():
    ibbr = IndexByReadId(_s)
    first_last = ("SRR043348.1017172", "SRR043396.9722121")
    for t in first_last:
        hits = ibbr.get_reads_by_id(t)
        n = 0
        for h in hits:
            n += 1
            assert_equal(t, h.query_name)
        assert_equal(n, 2)

def test_sort_creates_index():
    if os.path.exists(_sidx):
        os.remove(_sidx)
    ibbr = IndexByReadId(_f)
    ibbr.sort_bam(outfile=_s, batch_size=100)
    assert_equal(ibbr.bam, _s)
    assert_equal(ibbr.index_file, _sidx)
    ibbr.create_index()
    t = "SRR043372.11318637"
    hits = ibbr.get_reads_by_id(t)
    n = 0
    for h in hits:
        n += 1
        assert_equal(t, h.query_name)
    assert_equal(n, 2)

if __name__ == '__main__':
    test_sort_input()
    test_unsorted_error()
    test_create_index()
    test_read_index()
    test_find_targets()
