from nose.tools import *
import os
import pysam
from collections import OrderedDict
from index_bam_by_read_id import *

_f = 'test_data/test.bam'
_s = "test_data/test_rid_sorted.bam"
_sidx = "test_data/test_rid_sorted.bam.ibbr"
_f2 = 'test_data/sim_test.bam'
_s2 = 'test_data/sim_test_rid_sorted.bam'
_s2idx = 'test_data/sim_test_rid_sorted.bam.ibbr'

def test_sort_input():
    for f in [_s, _s2, _sidx, _s2idx]:
        if os.path.exists(f):
            os.remove(f)
    ibbr = IndexByReadId(_f)
    ibbr.sort_bam(outfile=_s, batch_size=100)
    ibbr2 = IndexByReadId(_f2)
    ibbr2.sort_bam()
    assert_true(os.path.exists(_s2))

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
    #import pdb; pdb.set_trace() 
    ibbr2.create_index()
    assert_true(os.path.exists(_s2idx))
    
def test_read_index():
    assert_true(os.path.exists(_sidx))
    ibbr = IndexByReadId(_s)
    idx = ibbr.read_index()
    isinstance(idx, OrderedDict)
    
def test_find_targets():
    ibbr = IndexByReadId(_s)
    ibbr.read_index()
    ibbr2 = IndexByReadId(_s2)
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
    for t in targets:
        hits = ibbr.get_reads_by_id(t)
        n = 0
        for h in hits:
            n += 1
            assert_equal(t, h.query_name)
        assert_equal(n, 2)
    for t in targets2[::-1]:
        hits = ibbr2.get_reads_by_id(t)
        n = 0
        for h in hits:
            n += 1
            assert_equal(t, h.query_name)
        assert_equal(n, 1)
    hits = ibbr2.get_reads_by_id( #search for non-existent read
                        "NZ_GG703879.1_967999_968061_0_1_0_0_4:0:0_1:0:0_62bb")
    assert_equal(len(hits), 0)
