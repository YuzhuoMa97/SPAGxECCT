"""Create VCF and BGEN 1.2 fixtures independently from the R-generated dosage table.
Only Python's standard library is used. Run from cpp/.
"""
from pathlib import Path
import csv, struct, zlib

def main():
    folder=Path('examples/synthetic')
    with (folder/'synthetic.genotypes.tsv').open() as f:
        rows=list(csv.reader(f,delimiter='\t'))
    ids=[r[0] for r in rows[1:]]
    dosage=[[None if r[j+1]=='NA' else int(r[j+1]) for r in rows[1:]] for j in range(len(rows[0])-1)]
    bim=[r.split() for r in (folder/'synthetic.bim').read_text().splitlines()]
    with (folder/'synthetic.vcf').open('w',newline='\n') as f:
        f.write('##fileformat=VCFv4.2\n##contig=<ID=1>\n##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(ids)+'\n')
        for b,g in zip(bim,dosage):
            f.write('\t'.join([b[0],b[3],b[1],b[5],b[4],'.','PASS','.','GT']+['./.' if x is None else ('0/0','0/1','1/1')[x] for x in g])+'\n')
    def s16(s):
        v=s.encode(); return struct.pack('<H',len(v))+v
    def s32(s):
        v=s.encode(); return struct.pack('<I',len(v))+v
    samples=struct.pack('<I',len(ids))+b''.join(s16(i) for i in ids)
    samples=struct.pack('<I',4+len(samples))+samples
    for codec in (0,1):
        for bits in (8,16):
            path=folder/('synthetic'+('' if codec else '.raw')+('' if bits==8 else '.16bit')+'.bgen')
            with path.open('wb') as f:
                f.write(struct.pack('<IIII4sI',20+len(samples),20,len(bim),len(ids),b'bgen',(1<<31)|8|codec)); f.write(samples)
                for b,g in zip(bim,dosage):
                    f.write(s16(b[1])+s16(b[1])+s16(b[0])+struct.pack('<IH',int(b[3]),2)+s32(b[5])+s32(b[4]))
                    limit=(1<<bits)-1
                    probs=[]
                    for x in g:
                        probs.extend((limit if x==0 else 0,limit if x==1 else 0))
                    raw=struct.pack('<IHBB',len(ids),2,2,2)+bytes(130 if x is None else 2 for x in g)+bytes((0,bits))+struct.pack('<' + ('B' if bits==8 else 'H')*len(probs),*probs)
                    if codec:
                        z=zlib.compress(raw);f.write(struct.pack('<II',4+len(z),len(raw))+z)
                    else:f.write(struct.pack('<I',len(raw))+raw)
    # Additional independent encodings: external sample IDs, reversed alleles,
    # phased layout 2 and legacy layout 1.1 with/without compression.
    for mode in ('noids','reflast','phased','v11','v11raw'):
        layout=4 if mode.startswith('v11') else 8
        embedded=mode not in ('noids','v11','v11raw')
        codec=0 if mode=='v11raw' else 1
        path=folder/('synthetic.'+mode+'.bgen')
        with path.open('wb') as f:
            f.write(struct.pack('<IIII4sI',20+(len(samples) if embedded else 0),20,len(bim),len(ids),b'bgen',(1<<31 if embedded else 0)|layout|codec))
            if embedded:f.write(samples)
            for b,g in zip(bim,dosage):
                if layout==4:f.write(struct.pack('<I',len(ids)))
                f.write(s16(b[1])+s16(b[1])+s16(b[0])+struct.pack('<I',int(b[3])))
                if layout==8:f.write(struct.pack('<H',2))
                f.write(s32(b[4] if mode=='reflast' else b[5])+s32(b[5] if mode=='reflast' else b[4]))
                if layout==4:
                    probs=[]
                    for x in g:probs.extend(tuple(32768 if x==i else 0 for i in range(3)))
                    raw=struct.pack('<'+'H'*len(probs),*probs)
                else:
                    probs=[]
                    for x in g:
                        if mode=='phased':probs.extend((255 if x in (0,1) else 0,255 if x==0 else 0))
                        else:probs.extend((255 if x==(2 if mode=='reflast' else 0) else 0,255 if x==1 else 0))
                    raw=struct.pack('<IHBB',len(ids),2,2,2)+bytes(130 if x is None else 2 for x in g)+bytes((int(mode=='phased'),8))+bytes(probs)
                if codec:
                    z=zlib.compress(raw)
                    f.write(struct.pack('<I',len(z)+(4 if layout==8 else 0)))
                    if layout==8:f.write(struct.pack('<I',len(raw)))
                    f.write(z)
                else:f.write(raw)
        if not embedded:
            path.with_suffix('.sample').write_text('ID_1 ID_2 missing\n0 0 0\n'+''.join(i+' '+i+' 0\n' for i in ids))
    print('Wrote VCF and 9 independent BGEN encodings (layouts 1/2, phased/unphased, raw/zlib, 8/16-bit, embedded/external IDs, REF/ALT modes).')
if __name__=='__main__':main()

