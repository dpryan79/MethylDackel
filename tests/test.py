from subprocess import check_call
import os
import os.path as op


def rm(f):
    try:
        os.unlink(f)
    except OSError:
        pass

rm('po')
check_call('cd .. && make', shell=True)
assert op.exists('../PileOMeth')

rm('ct_aln_CpG.bedGraph')
check_call('../PileOMeth extract ct100.fa ct_aln.bam -q 2', shell=True)
assert op.exists('ct_aln_CpG.bedGraph')
lines = sum(1 for _ in open('ct_aln_CpG.bedGraph'))
assert lines == 1


rm('cg_aln_CpG.bedGraph')
check_call('../PileOMeth extract cg100.fa cg_aln.bam -q 2', shell=True)
assert op.exists('cg_aln_CpG.bedGraph')
lines = sum(1 for _ in open('cg_aln_CpG.bedGraph'))
assert lines > 1

# should be none with q > 10
check_call('../PileOMeth extract cg100.fa cg_aln.bam -q 10', shell=True)
assert op.exists('cg_aln_CpG.bedGraph')
lines = sum(1 for _ in open('cg_aln_CpG.bedGraph'))
assert lines == 1
