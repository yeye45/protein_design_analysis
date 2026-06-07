#!/usr/bin/env python3
from __future__ import annotations
import argparse, logging, os, re, tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

LOGGER=logging.getLogger('generate_a3m')
SAFE_ID_RE=re.compile(r'[^A-Za-z0-9_.-]+')
AA_RE=re.compile(r'^[A-Za-z*.-]+$')
@dataclass(frozen=True)
class FastaRecord:
    seq_id:str; description:str; sequence:str
@dataclass(frozen=True)
class A3MRecord:
    header:str; sequence:str

def parse_args():
    p=argparse.ArgumentParser(description='Create mutant A3M files by replacing the WT A3M query sequence.')
    p.add_argument('--wt-a3m',required=True,type=Path); p.add_argument('--mutants-fasta',required=True,type=Path); p.add_argument('--out-dir',required=True,type=Path)
    p.add_argument('--num-batches',type=int,default=5); p.add_argument('--wrap',type=int,default=100); p.add_argument('--skip-existing',action='store_true')
    p.add_argument('--rename-query-header',action='store_true'); p.add_argument('--allow-length-mismatch',action='store_true'); p.add_argument('--allow-duplicate-ids',action='store_true')
    p.add_argument('--manifest-name',default='manifest.tsv'); p.add_argument('--log-level',default='INFO',choices=['DEBUG','INFO','WARNING','ERROR'])
    return p.parse_args()
def setup_logging(level): logging.basicConfig(level=getattr(logging,level),format='%(asctime)s %(levelname)s %(name)s: %(message)s')
def sanitize_id(seq_id):
    safe=SAFE_ID_RE.sub('_',seq_id.strip()).strip('._-')
    if not safe: raise ValueError(f'FASTA ID {seq_id!r} becomes empty after sanitization')
    return safe
def wrap_sequence(seq,width): return seq+'\n' if width<=0 else ''.join(seq[i:i+width]+'\n' for i in range(0,len(seq),width))
def parse_fasta(path:Path)->Iterable[FastaRecord]:
    seq_id=None; desc=''; chunks=[]
    with path.open('r',encoding='utf-8') as h:
        for line_no,raw in enumerate(h,1):
            line=raw.strip()
            if not line: continue
            if line.startswith('>'):
                if seq_id is not None: yield FastaRecord(seq_id,desc,''.join(chunks).replace(' ',''))
                desc=line[1:].strip()
                if not desc: raise ValueError(f'Empty FASTA header at {path}:{line_no}')
                seq_id=desc.split()[0]; chunks=[]
            else:
                if seq_id is None: raise ValueError(f'Sequence line before FASTA header at {path}:{line_no}')
                chunks.append(line)
    if seq_id is not None: yield FastaRecord(seq_id,desc,''.join(chunks).replace(' ',''))
def parse_a3m(path:Path):
    records=[]; header=None; chunks=[]
    with path.open('r',encoding='utf-8') as h:
        for line_no,raw in enumerate(h,1):
            line=raw.rstrip('\n\r')
            if not line: continue
            if line.startswith('>'):
                if header is not None: records.append(A3MRecord(header,''.join(chunks)))
                header=line; chunks=[]
            else:
                if header is None: raise ValueError(f'A3M sequence line before header at {path}:{line_no}')
                chunks.append(line)
    if header is not None: records.append(A3MRecord(header,''.join(chunks)))
    if not records: raise ValueError(f'No A3M records found in {path}')
    return records
def ungapped_upper_len(seq): return len([c for c in seq if c.isalpha() and c.isupper()])
def validate_mutant_sequence(rec):
    seq=rec.sequence.strip().replace(' ','').replace('\t','').upper()
    if not seq: raise ValueError(f'Mutant {rec.seq_id!r} has an empty sequence')
    if not AA_RE.match(seq): raise ValueError(f'Mutant {rec.seq_id!r} contains unsupported characters')
    return seq
def render_a3m(records,mutant_id,mutant_seq,wrap,rename_query):
    first_header=f'>{mutant_id}' if rename_query else records[0].header
    parts=[first_header+'\n',wrap_sequence(mutant_seq,wrap)]
    for rec in records[1:]: parts.extend([rec.header+'\n',wrap_sequence(rec.sequence,wrap)])
    return ''.join(parts)
def atomic_write(path,content):
    path.parent.mkdir(parents=True,exist_ok=True); fd,tmp=tempfile.mkstemp(prefix=path.name+'.',suffix='.tmp',dir=str(path.parent))
    try:
        with os.fdopen(fd,'w',encoding='utf-8') as h: h.write(content); h.flush(); os.fsync(h.fileno())
        os.replace(tmp,path)
    except Exception:
        try: os.unlink(tmp)
        except OSError: pass
        raise
def main():
    a=parse_args(); setup_logging(a.log_level)
    if a.num_batches<1: raise ValueError('--num-batches must be >= 1')
    if not a.wt_a3m.is_file(): raise FileNotFoundError(a.wt_a3m)
    if not a.mutants_fasta.is_file(): raise FileNotFoundError(a.mutants_fasta)
    records=parse_a3m(a.wt_a3m); wt_len=ungapped_upper_len(records[0].sequence); LOGGER.info('Loaded WT A3M records=%d query_length=%d',len(records),wt_len)
    a.out_dir.mkdir(parents=True,exist_ok=True); batch_dirs=[a.out_dir/f'batch_{i}' for i in range(a.num_batches)]
    for d in batch_dirs: d.mkdir(parents=True,exist_ok=True)
    seen={}; written=0; skipped=0
    with (a.out_dir/a.manifest_name).open('w',encoding='utf-8') as m:
        m.write('index\tsequence_id\tbatch\tsequence_length\ta3m_path\n')
        for idx,rec in enumerate(parse_fasta(a.mutants_fasta)):
            sid=sanitize_id(rec.seq_id)
            if sid in seen:
                if not a.allow_duplicate_ids: raise ValueError(f'Duplicate FASTA ID after sanitization: {sid!r}')
                seen[sid]+=1; sid=f'{sid}__dup{seen[sid]}'
            else: seen[sid]=0
            seq=validate_mutant_sequence(rec)
            if not a.allow_length_mismatch and len(seq)!=wt_len: raise ValueError(f'Length mismatch for {rec.seq_id!r}: mutant={len(seq)} WT_query={wt_len}')
            batch=idx%a.num_batches; out=batch_dirs[batch]/f'{sid}.a3m'
            if a.skip_existing and out.exists(): skipped+=1
            else: atomic_write(out,render_a3m(records,sid,seq,a.wrap,a.rename_query_header)); written+=1
            m.write(f'{idx}\t{sid}\t{batch}\t{len(seq)}\t{out.as_posix()}\n')
            if (idx+1)%10000==0: LOGGER.info('Processed %d mutants; written=%d skipped=%d',idx+1,written,skipped)
    LOGGER.info('Done. written=%d skipped=%d',written,skipped)
    return 0
if __name__=='__main__':
    try: raise SystemExit(main())
    except Exception as exc: LOGGER.exception('Failed: %s',exc); raise SystemExit(1)
