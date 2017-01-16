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

rm('ct_aln_CpG.bedGraph')
check_call('../MethylDackel extract ct100.fa ct_aln.bam -q 2', shell=True)
assert op.exists('ct_aln_CpG.bedGraph')
lines = sum(1 for _ in open('ct_aln_CpG.bedGraph'))
assert lines == 1
rm('ct_aln_CpG.bedGraph')

rm('cg_aln_CpG.bedGraph')
check_call('../MethylDackel extract cg100.fa cg_aln.bam -q 2', shell=True)
assert op.exists('cg_aln_CpG.bedGraph')
lines = sum(1 for _ in open('cg_aln_CpG.bedGraph'))
assert lines > 1
rm('cg_aln_CpG.bedGraph')

# should be none with q > 10
check_call('../MethylDackel extract cg100.fa cg_aln.bam -q 10', shell=True)
assert op.exists('cg_aln_CpG.bedGraph')
lines = sum(1 for _ in open('cg_aln_CpG.bedGraph'))
assert lines == 1
rm('cg_aln_CpG.bedGraph')

# Test the new methylKit option
check_call('../MethylDackel extract --methylKit --CHH --CHG cg100.fa cg_aln.bam -q 2', shell=True)
assert op.exists('cg_aln_CpG.methylKit')
lines = sum(1 for _ in open('cg_aln_CpG.methylKit'))
assert lines > 1
rm('cg_aln_CpG.methylKit')
lines = sum(1 for _ in open('cg_aln_CHG.methylKit'))
assert lines == 1
rm('cg_aln_CHG.methylKit')
assert op.exists('cg_aln_CHH.methylKit')
lines = sum(1 for _ in open('cg_aln_CHH.methylKit'))
assert lines == 1
rm('cg_aln_CHH.methylKit')

# Check that --minDepth is working, which means there should be no called sites
check_call('../MethylDackel extract --minDepth 2 cg100.fa cg_aln.bam -q 2', shell=True)
assert op.exists('cg_aln_CpG.bedGraph')
lines = sum(1 for _ in open('cg_aln_CpG.bedGraph'))
assert lines == 1
rm('cg_aln_CpG.bedGraph')

# Check that --ignoreFlags is working, which means that there are now called sites
check_call('../MethylDackel extract --ignoreFlags 0xD00 cg100.fa cg_aln.bam -q 2', shell=True)
assert op.exists('cg_aln_CpG.bedGraph')
lines = sum(1 for _ in open('cg_aln_CpG.bedGraph'))
assert lines == 49
rm('cg_aln_CpG.bedGraph')

<<<<<<< HEAD
# Check that --requireFlags is working
check_call('../PileOMeth extract --requireFlags 0xD00 cg100.fa cg_aln.bam -q 2', shell=True)
assert op.exists('cg_aln_CpG.bedGraph')
lines = sum(1 for _ in open('cg_aln_CpG.bedGraph'))
assert lines == 49
=======
# Check absolute trimming bounds
check_call('../MethylDackel extract --nOT 50,50,40,40 cg100.fa cg_aln.bam -q 2', shell=True)
assert op.exists('cg_aln_CpG.bedGraph')
lines = sum(1 for _ in open('cg_aln_CpG.bedGraph'))
assert lines == 12
>>>>>>> master
rm('cg_aln_CpG.bedGraph')
