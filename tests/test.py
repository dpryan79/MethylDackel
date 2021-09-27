from subprocess import check_call
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
49 for --ignoreNH
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
print("Finished correctly")

