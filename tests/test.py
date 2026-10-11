from subprocess import check_call, check_output, PIPE, Popen

import os
import os.path as op
wd = op.dirname(op.realpath(__file__))
os.chdir(wd)

def rm(f):
    try:
        os.unlink(f)
    except OSError:
        pass

assert op.exists('../MethylDackel')
MPath = os.path.abspath('../MethylDackel')
check_call([MPath, '--version'])

rm('test1_CpG.bedGraph')
check_call([MPath, 'extract', 'ct100.fa', 'ct_aln.bam', '-q', '2', '-o', 'test1'])
assert op.exists('test1_CpG.bedGraph')
lines = sum(1 for _ in open('test1_CpG.bedGraph'))
assert lines == 1
rm('test1_CpG.bedGraph')

rm('test2_CpG.bedGraph')
check_call([MPath, 'extract', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test2'])
assert op.exists('test2_CpG.bedGraph')
for line in open('test2_CpG.bedGraph'):
    print(line)
lines = sum(1 for _ in open('test2_CpG.bedGraph'))
assert lines > 1
rm('test2_CpG.bedGraph')

# should be none with q > 10
rm('test3_CpG.bedGraph')
check_call([MPath, 'extract', 'cg100.fa', 'cg_aln.bam', '-q', '10', '-o', 'test3'])
assert op.exists('test3_CpG.bedGraph')
lines = sum(1 for _ in open('test3_CpG.bedGraph'))
assert lines == 1
rm('test3_CpG.bedGraph')

# Test the new methylKit option
rm('test4_CpG.methylKit')
rm('test4_CHG.methylKit')
rm('test4_CHH.methylKit')
check_call([MPath, 'extract', '--methylKit', '--CHH', '--CHG', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test4'])
assert op.exists('test4_CpG.methylKit')
lines = sum(1 for _ in open('test4_CpG.methylKit'))
assert lines > 1
rm('test4_CpG.methylKit')
lines = sum(1 for _ in open('test4_CHG.methylKit'))
assert lines == 1
rm('test4_CHG.methylKit')
assert op.exists('test4_CHH.methylKit')
lines = sum(1 for _ in open('test4_CHH.methylKit'))
assert lines == 2
rm('test4_CHH.methylKit')

# Check that --minDepth is working, which means there should be no called sites
rm('test5_CpG.bedGraph')
check_call([MPath, 'extract', '--minDepth', '2', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test5'])
assert op.exists('test5_CpG.bedGraph')
lines = sum(1 for _ in open('test5_CpG.bedGraph'))
assert lines == 1
rm('test5_CpG.bedGraph')

# Check that --ignoreFlags is working, which means that there are now called sites
rm('test6_CpG.bedGraph')
check_call([MPath, 'extract', '--ignoreFlags', '0xD00', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test6'])
assert op.exists('test6_CpG.bedGraph')
lines = sum(1 for _ in open('test6_CpG.bedGraph'))
assert lines == 49
rm('test6_CpG.bedGraph')

# Check that --requireFlags is working
rm('test7_CpG.bedGraph')
check_call([MPath, 'extract', '--requireFlags', '0xD00', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test7'])
assert op.exists('test7_CpG.bedGraph')
lines = sum(1 for _ in open('test7_CpG.bedGraph'))
assert lines == 49
rm('test7_CpG.bedGraph')

# Check absolute trimming bounds
rm('test8_CpG.bedGraph')
check_call([MPath, 'extract', '--nOT', '50,50,40,40', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test8'])
assert op.exists('test8_CpG.bedGraph')
lines = sum(1 for _ in open('test8_CpG.bedGraph'))
assert lines == 11
rm('test8_CpG.bedGraph')

# Excluding one fewer base from each 3' end keeps one more CpG
rm('test8_CpG.bedGraph')
check_call([MPath, 'extract', '--nOT', '50,49,40,39', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test8'])
assert op.exists('test8_CpG.bedGraph')
lines = sum(1 for _ in open('test8_CpG.bedGraph'))
assert lines == 12
rm('test8_CpG.bedGraph')

# Check variant filtering (there are 49 lines otherwise)
rm('test9_CpG.bedGraph')
check_call([MPath, 'extract', '-p', '1', '-q', '0', '-o', 'test9', '--minOppositeDepth', '3', '--maxVariantFrac', '0.25', 'cg100.fa', 'cg_with_variants.bam'])
assert op.exists('test9_CpG.bedGraph')
lines = sum(1 for _ in open('test9_CpG.bedGraph'))
assert lines == 48
rm('test9_CpG.bedGraph')

# Check conversion efficiency. 2 read pairs, one mostly converted
# By default, 1 read is MAPQ filtered, another is kept
rm('test10_CpG.bedGraph')
check_call([MPath, 'extract', '-o', 'test10', 'chgchh.fa', 'chgchh_aln.bam'])
assert op.exists('test10_CpG.bedGraph')
lines = sum(1 for _ in open('test10_CpG.bedGraph'))
assert lines == 2
rm('test10_CpG.bedGraph')

# Ensure 2 reads/positions are covered by changing MAPQ
rm('test11_CpG.bedGraph')
check_call([MPath, 'extract', '-o', 'test11', '-q', '5', 'chgchh.fa', 'chgchh_aln.bam'])
assert op.exists('test11_CpG.bedGraph')
lines = sum(1 for _ in open('test11_CpG.bedGraph'))
assert lines == 3
rm('test11_CpG.bedGraph')

# Only 1 read has a conversion efficiency >=0.9
rm('test12_CpG.bedGraph')
check_call([MPath, 'extract', '-o', 'test12', '-q', '5', '--minConversionEfficiency', '0.9', 'chgchh.fa', 'chgchh_aln.bam'])
assert op.exists('test12_CpG.bedGraph')
lines = sum(1 for _ in open('test12_CpG.bedGraph'))
assert lines == 2
rm('test12_CpG.bedGraph')

# No perfectly converted reads
rm('test13_CpG.bedGraph')
check_call([MPath, 'extract', '-o', 'test13', '-q', '5', '--minConversionEfficiency', '1.0', 'chgchh.fa', 'chgchh_aln.bam'])
assert op.exists('test13_CpG.bedGraph')
lines = sum(1 for _ in open('test13_CpG.bedGraph'))
assert lines == 1
rm('test13_CpG.bedGraph')

# Test ignoreNH
rm('test14_CpG.bedGraph')
check_call([MPath, 'extract', '-o', 'test14', '-q', '1', 'cg100.fa', 'NH.bam'])
assert op.exists('test14_CpG.bedGraph')
lines = sum(1 for _ in open('test14_CpG.bedGraph'))
assert lines == 1
rm('test14_CpG.bedGraph')

# Test ignoreNH
rm('test15_CpG.bedGraph')
check_call([MPath, 'extract', '-o', 'test15', '--ignoreNH', '-q', '1', 'cg100.fa', 'NH.bam'])
assert op.exists('test15_CpG.bedGraph')
lines = sum(1 for _ in open('test15_CpG.bedGraph'))
assert lines == 49
rm('test15_CpG.bedGraph')

# Bounds strings may contain 0 in any field, and the fields after a 0 are still applied. The original
# failure depended on a stale errno at parse time, so these reproduce it only where errno is non-zero by then.
rm('test16_CpG.bedGraph')
check_call([MPath, 'extract', '--OT', '0,50,0,40', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test16'])
assert op.exists('test16_CpG.bedGraph')
lines = sum(1 for _ in open('test16_CpG.bedGraph'))
assert lines == 25
rm('test16_CpG.bedGraph')

rm('test17_CpG.bedGraph')
check_call([MPath, 'extract', '--nOT', '0,50,0,40', 'cg100.fa', 'cg_aln.bam', '-q', '2', '-o', 'test17'])
assert op.exists('test17_CpG.bedGraph')
lines = sum(1 for _ in open('test17_CpG.bedGraph'))
assert lines == 30
rm('test17_CpG.bedGraph')

# Malformed bounds strings are an error (rather than being partially applied or crashing)
def bounds_error(args):
    p = Popen(args, stdout=PIPE, stderr=PIPE)
    _, err = p.communicate()
    return p.returncode == 1 and b'Invalid bounds string' in err

malformed = ['5', '5,5', '5,0,0,0,0', '5,,0,0', '5,,0,0,0', ',5,0,0,0', '5,0,0,0,', 'a,0,0,0', '-1,0,0,0',
             '5x,0,0,0', '99999999999,0,0,0']
for bounds in malformed:
    assert bounds_error([MPath, 'extract', '--nOT', bounds, '-o', 'test18', 'cg100.fa', 'cg_aln.bam'])
    assert bounds_error([MPath, 'mbias', '--nOT', bounds, 'cg100.fa', 'cg_aln.bam', 'test18'])
for opt in ['OT', 'OB', 'CTOT', 'CTOB', 'nOT', 'nOB', 'nCTOT', 'nCTOB']:
    assert bounds_error([MPath, 'extract', '--' + opt, '5,,0,0,0', '-o', 'test18', 'cg100.fa', 'cg_aln.bam'])
for opt in ['nOT', 'nOB', 'nCTOT', 'nCTOB']:
    assert bounds_error([MPath, 'mbias', '--' + opt, '5,,0,0,0', 'cg100.fa', 'cg_aln.bam', 'test18'])

# With --mergeContext, a CpG must be dropped whether the likely variant is at its C or its G.
# cg_variants_merged.bam has a C>T at the C of the CpG at 10 and a G>A at the G of the CpG
# at 20, each seen on 10 reads of the opposite strand (49 CpGs without the filter).
rm('test18_CpG.bedGraph')
check_call([MPath, 'extract', '-p', '1', '-q', '0', '-o', 'test18', '--mergeContext', '--minOppositeDepth', '3', '--maxVariantFrac', '0.25', 'cg100.fa', 'cg_variants_merged.bam'])
assert op.exists('test18_CpG.bedGraph')
starts = [l.split('\t')[1] for l in open('test18_CpG.bedGraph') if not l.startswith('track')]
assert len(starts) == 47
assert '10' not in starts and '20' not in starts
rm('test18_CpG.bedGraph')

# Test min/max InsertSize
rm('test19_CpG.bedGraph')
check_call([MPath, 'extract', '-o', 'test16', '--minInsertSize', '0', '--maxInsertSize', '200', '-q', '1', 'cg100.fa', 'cg_aln.bam'])
assert op.exists('test16_CpG.bedGraph')
lines = sum(1 for _ in open('test16_CpG.bedGraph'))
assert lines == 49
rm('test16_CpG.bedGraph')

# Test mbias functionality
print("Testing mbias functionality...")

# Test basic mbias functionality with relaxed quality settings
result = check_output([MPath, 'mbias', '-q', '0', '-p', '1', '--txt', '--noSVG', 'cg100.fa', 'cg_aln.bam'])
assert b'Strand' in result and b'Position' in result and b'OT' in result
print("Basic mbias test passed")

# Test mbias with insert size filtering (should include reads with ~100bp insert size)
print("Testing mbias with insert size filtering (inclusive range)...")
result = check_output([MPath, 'mbias', '-q', '0', '-p', '1', '--minInsertSize', '90', '--maxInsertSize', '110', '--txt', '--noSVG', 'cg100.fa', 'cg_aln.bam'])
# Should have same results as basic test since insert size ~100 is included
assert b'Strand' in result and b'Position' in result and b'OT' in result
print("mbias insert size filtering (inclusive) test passed")

# Test mbias with restrictive insert size filtering (should exclude all reads)
print("Testing mbias with restrictive insert size filtering...")
result = check_output([MPath, 'mbias', '-q', '0', '-p', '1', '--minInsertSize', '200', '--maxInsertSize', '300', '--txt', '--noSVG', 'cg100.fa', 'cg_aln.bam'])
# Should only have header since no reads match the insert size range
assert b'Strand' in result and b'Position' in result and b'OT' not in result
print("mbias restrictive insert size filtering test passed")

print("Finished correctly")